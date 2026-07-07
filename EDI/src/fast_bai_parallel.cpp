#include <RcppEigen.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Fast Bai Adjusted T Statistic for Multiple Permutations
//'
//' @param w_mat Integer matrix of permuted treatment assignments (n x r).
//' @param m_mat Integer matrix of match indicators (n x r).
//' @param y Numeric response vector.
//' @param delta Null treatment effect shift.
//' @param halves_idx Integer matrix of half-sample indices.
//' @param convex_flag Logical flag for convex combination.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of Bai adjusted T statistics.
// [[Rcpp::export]]
NumericVector compute_bai_distr_parallel_cpp(
    SEXP w_mat_sexp,
    SEXP m_mat_sexp,
    SEXP y_sexp,
    double delta,
    SEXP halves_idx_sexp,
    bool convex_flag,
    int num_cores) {

  Eigen::Map<const Eigen::MatrixXi> w_mat(INTEGER(w_mat_sexp), Rf_nrows(w_mat_sexp), Rf_ncols(w_mat_sexp));
  Eigen::Map<const Eigen::MatrixXi> m_mat(INTEGER(m_mat_sexp), Rf_nrows(m_mat_sexp), Rf_ncols(m_mat_sexp));
  Eigen::Map<const Eigen::VectorXd> y(REAL(y_sexp), Rf_length(y_sexp));
  Eigen::Map<const Eigen::MatrixXi> halves_idx(INTEGER(halves_idx_sexp), Rf_nrows(halves_idx_sexp), Rf_ncols(halves_idx_sexp));

  int n = y.size();
  int nsim = w_mat.cols();
  int n_halves = halves_idx.rows();
  std::vector<double> results_vec(nsim);

  const double* y_ptr = y.data();
  const int* w_ptr = w_mat.data();
  const int* m_ptr = m_mat.data();
  const int* h_ptr = halves_idx.data();
  double* res_ptr = results_vec.data();

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel
  {
    std::vector<double> d_i;
    std::vector<double> y_r;
    std::vector<int> w_r;
    std::vector<double> match_T, match_C;
    std::vector<char> has_T, has_C;
    std::vector<double> yT, yC;

    d_i.reserve(n / 2);
    y_r.reserve(n);
    w_r.reserve(n);
    yT.reserve(n);
    yC.reserve(n);

#pragma omp for schedule(static)
    for (int b = 0; b < nsim; ++b) {
      const int* w_col = w_ptr + (size_t)b * n;
      const int* m_col = m_ptr + (size_t)b * n;

      d_i.clear();
      y_r.clear();
      w_r.clear();
      yT.clear();
      yC.clear();

      int max_match = 0;
      for (int i = 0; i < n; ++i) if (m_col[i] > max_match) max_match = m_col[i];

      if (max_match > 0) {
        match_T.assign(max_match, 0.0);
        match_C.assign(max_match, 0.0);
        has_T.assign(max_match, 0);
        has_C.assign(max_match, 0);
        for (int i = 0; i < n; ++i) {
          double y_val = y_ptr[i] + (w_col[i] == 1 ? delta : 0.0);
          int m = m_col[i];
          if (m > 0) {
            if (w_col[i] == 1) { match_T[m-1] = y_val; has_T[m-1] = 1; }
            else { match_C[m-1] = y_val; has_C[m-1] = 1; }
          } else {
            y_r.push_back(y_val);
            w_r.push_back(w_col[i]);
          }
        }
        for (int m = 0; m < max_match; ++m) if (has_T[m] && has_C[m]) d_i.push_back(match_T[m] - match_C[m]);
      } else {
        for (int i = 0; i < n; ++i) {
          y_r.push_back(y_ptr[i] + (w_col[i] == 1 ? delta : 0.0));
          w_r.push_back(w_col[i]);
        }
      }

      int m_size = d_i.size();
      if (m_size == 0 && y_r.empty()) { res_ptr[b] = NA_REAL; continue; }

      double d_bar = 0;
      if (m_size > 0) {
        for (double d : d_i) d_bar += d;
        d_bar /= m_size;
      }

      double r_bar = 0;
      double ssqR = 0;
      int nRT = 0, nRC = 0;
      if (!y_r.empty()) {
        double sumT = 0, sumC = 0;
        for (size_t i = 0; i < y_r.size(); ++i) {
          if (w_r[i] == 1) { sumT += y_r[i]; nRT++; yT.push_back(y_r[i]); }
          else { sumC += y_r[i]; nRC++; yC.push_back(y_r[i]); }
        }
        if (nRT > 0 && nRC > 0) {
          r_bar = (sumT/nRT) - (sumC/nRC);
          if (nRT > 1 && nRC > 1) {
            double varT = 0, varC = 0;
            double meanT = sumT/nRT, meanC = sumC/nRC;
            for (double val : yT) varT += (val - meanT)*(val - meanT);
            for (double val : yC) varC += (val - meanC)*(val - meanC);
            ssqR = (varT/(nRT-1))/nRT + (varC/(nRC-1))/nRC;
          }
        }
      }

      double bai_var_d_bar = 0;
      if (m_size > 0) {
        double delta_sq = d_bar * d_bar;
        double tau_sq = 0;
        for (double d : d_i) tau_sq += d * d;
        tau_sq /= m_size;
        double lambda_squ = 0;
        if (n_halves > 0) {
          for (int i = 0; i < n_halves; ++i) {
            int id1 = h_ptr[i]; // Pair ID (1-indexed)
            int id2 = h_ptr[i + n_halves]; // Pair ID (1-indexed)

            if (id1 > 0 && id1 <= max_match && id2 > 0 && id2 <= max_match) {
              if (has_T[id1-1] && has_C[id1-1] && has_T[id2-1] && has_C[id2-1]) {
                lambda_squ += (match_T[id1-1] - match_C[id1-1]) * (match_T[id2-1] - match_C[id2-1]);
              }
            }
          }
          lambda_squ /= n_halves;
        }
        bai_var_d_bar = std::max(1e-8, tau_sq - (lambda_squ + delta_sq) / 2.0) / m_size;
      }

      if (convex_flag && nRT > 1 && nRC > 1 && m_size > 0 && ssqR > 0) {
        double w_star = ssqR / (ssqR + bai_var_d_bar);
        res_ptr[b] = w_star * d_bar + (1 - w_star) * r_bar;
      } else if (m_size > 0) {
        res_ptr[b] = d_bar;
      } else {
        res_ptr[b] = r_bar;
      }
    }
  }
  return wrap(results_vec);
}
