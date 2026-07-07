#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

//' Pocock-Simon Minimization Allocation Logic
//'
//' @param counts A matrix of dimensions (sum of levels) x (number of treatments).
//'               Each row corresponds to a specific level of a specific covariate.
//' @param subject_levels_idx An integer vector of indices indicating which rows 
//'                           of the counts matrix the current subject belongs to.
//' @param weights A numeric vector of weights for each covariate.
//' @param p_best The probability of assigning the treatment that minimizes the imbalance.
//' @param prob_T Target probability for treatment (usually 0.5).
//'
//' @return The assigned treatment (0 or 1).
//' @export
//' @keywords internal
// [[Rcpp::export]]
int pocock_simon_assign_cpp(NumericMatrix counts, IntegerVector subject_levels_idx, NumericVector weights, double p_best, double prob_T) {
  int num_trts = counts.cols(); // Should be 2 (Control=0, Treatment=1)
  int num_covs = subject_levels_idx.size();
  
  std::vector<double> G(num_trts, 0.0);
  
  for (int k = 0; k < num_trts; ++k) {
    double G_k = 0.0;
    for (int j = 0; j < num_covs; ++j) {
      int row_idx = subject_levels_idx[j] - 1; // 1-based to 0-based
      
      // Calculate imbalance for covariate j if assigned to treatment k
      // We use the variance of counts across treatments as the measure d_ik
      std::vector<double> counts_after(num_trts);
      double sum = 0.0;
      for (int t = 0; t < num_trts; ++t) {
        counts_after[t] = counts(row_idx, t) + (t == k ? 1 : 0);
        sum += counts_after[t];
      }
      
      double mean = sum / num_trts;
      double var = 0.0;
      for (int t = 0; t < num_trts; ++t) {
        var += (counts_after[t] - mean) * (counts_after[t] - mean);
      }
      var /= (num_trts - 1);
      
      G_k += weights[j] * var;
    }
    G[k] = G_k;
  }
  
  int best_trt = 0;
  if (G[1] < G[0]) {
    best_trt = 1;
  } else if (G[0] < G[1]) {
    best_trt = 0;
  } else {
    // Tie: use prob_T
    return (R::unif_rand() < prob_T) ? 1 : 0;
  }
  
  // Assign best treatment with probability p_best
  if (R::unif_rand() < p_best) {
    return best_trt;
  } else {
    return 1 - best_trt;
  }
}

//' Pocock-Simon Minimization Assignment and Update
//'
//' @param counts A matrix of dimensions (sum of levels) x (number of treatments).
//'               Modified in place.
//' @param subject_levels_idx An integer vector of indices indicating which rows 
//'                           of the counts matrix the current subject belongs to.
//' @param weights A numeric vector of weights for each covariate.
//' @param p_best The probability of assigning the treatment that minimizes the imbalance.
//' @param prob_T Target probability for treatment (usually 0.5).
//'
//' @return The assigned treatment (0 or 1).
//' @export
//' @keywords internal
// [[Rcpp::export]]
int pocock_simon_assign_and_update_cpp(NumericMatrix counts, IntegerVector subject_levels_idx, NumericVector weights, double p_best, double prob_T) {
  int w = pocock_simon_assign_cpp(counts, subject_levels_idx, weights, p_best, prob_T);
  
  // Update counts in place
  for (int j = 0; j < subject_levels_idx.size(); ++j) {
    int row_idx = subject_levels_idx[j] - 1;
    counts(row_idx, w)++;
  }
  
  return w;
}

//' Pocock-Simon Minimization Redraw Assignments
//'
//' @param x_levels_matrix A matrix where each row is a subject and each column is the row 
//'   index in counts for that covariate.
//' @param num_levels_total Total number of levels across all covariates.
//' @param weights A numeric vector of weights for each covariate.
//' @param p_best The probability of assigning the treatment that minimizes the imbalance.
//' @param prob_T Target probability for treatment (usually 0.5).
//'
//' @return An integer vector of treatment assignments.
//' @export
//' @keywords internal
// [[Rcpp::export]]
IntegerVector pocock_simon_redraw_w_cpp(IntegerMatrix x_levels_matrix, int num_levels_total, NumericVector weights, double p_best, double prob_T) {
  int n = x_levels_matrix.nrow();
  int num_covs = x_levels_matrix.ncol();
  if (num_levels_total <= 0) stop("num_levels_total must be positive");
  if (weights.size() != num_covs) stop("weights length must match the number of covariates");

  std::vector<int> level_rows(static_cast<std::size_t>(n) * num_covs);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < num_covs; ++j) {
      const int row_idx = x_levels_matrix(i, j) - 1;
      if (row_idx < 0 || row_idx >= num_levels_total) {
        stop("x_levels_matrix contains a level index outside 1..num_levels_total");
      }
      level_rows[static_cast<std::size_t>(i) * num_covs + j] = row_idx;
    }
  }

  IntegerVector w(n);
  int* w_ptr = w.begin();
  const double* weights_ptr = weights.begin();
  std::vector<int> counts(static_cast<std::size_t>(num_levels_total) * 2, 0);
  
  for (int i = 0; i < n; ++i) {
    const int* subject_levels = level_rows.data() + static_cast<std::size_t>(i) * num_covs;
    double G[2] = {0.0, 0.0};
    for (int k = 0; k < 2; ++k) {
      double G_k = 0.0;
      for (int j = 0; j < num_covs; ++j) {
        const int row_idx = subject_levels[j];
        const double c0 = counts[static_cast<std::size_t>(row_idx) * 2] + (k == 0 ? 1 : 0);
        const double c1 = counts[static_cast<std::size_t>(row_idx) * 2 + 1] + (k == 1 ? 1 : 0);
        const double mean = (c0 + c1) / 2.0;
        const double variance = (c0 - mean) * (c0 - mean) + (c1 - mean) * (c1 - mean);
        G_k += weights_ptr[j] * variance;
      }
      G[k] = G_k;
    }

    int assigned_w;
    if (G[0] == G[1]) {
      assigned_w = (R::unif_rand() < prob_T) ? 1 : 0;
    } else {
      const int best_treatment = (G[1] < G[0]) ? 1 : 0;
      assigned_w = (R::unif_rand() < p_best) ? best_treatment : 1 - best_treatment;
    }
    w_ptr[i] = assigned_w;
    
    // Update counts
    for (int j = 0; j < num_covs; ++j) {
      counts[static_cast<std::size_t>(subject_levels[j]) * 2 + assigned_w]++;
    }
  }
  
  return w;
}
