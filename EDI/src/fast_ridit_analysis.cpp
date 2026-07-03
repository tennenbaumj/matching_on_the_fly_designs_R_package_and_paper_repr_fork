#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace Eigen;

namespace {

struct RiditScoreData {
    NumericVector scores;
    std::vector<int> levels;
    NumericVector ref_p;
};

RiditScoreData fast_ridit_scores_from_ref_indices(const int* y_ptr,
                                                  int n,
                                                  const std::vector<int>& ref_idx_zero_based) {
    const int n_ref = static_cast<int>(ref_idx_zero_based.size());

    // 1. Get unique levels and their counts in the reference group
    std::vector<int> levels;
    std::vector<int> counts;
    levels.reserve(static_cast<std::size_t>(std::min(n_ref, 16)));
    counts.reserve(levels.capacity());
    for (int idx : ref_idx_zero_based) {
        const int yi = y_ptr[idx];
        auto level_it = std::lower_bound(levels.begin(), levels.end(), yi);
        const std::size_t pos = static_cast<std::size_t>(level_it - levels.begin());
        if (level_it != levels.end() && *level_it == yi) {
            counts[pos]++;
        } else {
            levels.insert(level_it, yi);
            counts.insert(counts.begin() + static_cast<std::ptrdiff_t>(pos), 1);
        }
    }

    const int K = static_cast<int>(levels.size());
    NumericVector ref_p(K);
    std::vector<double> ridit_values(static_cast<std::size_t>(K));

    // 2. Calculate Ridit scores for each category
    // R_k = sum_{j < k} p_j + 0.5 * p_k
    double cumulative_p = 0.0;
    for (int k = 0; k < K; ++k) {
        const double pk = static_cast<double>(counts[static_cast<std::size_t>(k)]) / n_ref;
        ref_p[k] = pk;
        const double ridit = cumulative_p + 0.5 * pk;
        ridit_values[static_cast<std::size_t>(k)] = ridit;
        cumulative_p += pk;
    }

    // 3. Assign scores to all subjects
    NumericVector scores(n);
    for (int i = 0; i < n; ++i) {
        const int yi = y_ptr[i];
        auto level_it = std::lower_bound(levels.begin(), levels.end(), yi);
        if (level_it != levels.end() && *level_it == yi) {
            scores[i] = ridit_values[static_cast<std::size_t>(level_it - levels.begin())];
        } else {
            // If a level wasn't in the reference group, find its place
            if (level_it == levels.begin()) {
                scores[i] = 0.0; // Extremely low
            } else if (level_it == levels.end()) {
                scores[i] = 1.0; // Extremely high
            } else {
                // Average of the ridits of the categories it falls between
                const int idx = static_cast<int>(std::distance(levels.begin(), level_it));
                scores[i] = 0.5 * (
                    ridit_values[static_cast<std::size_t>(idx - 1)] +
                    ridit_values[static_cast<std::size_t>(idx)]
                );
            }
        }
    }

    return RiditScoreData{scores, levels, ref_p};
}

}  // namespace

// [[Rcpp::export]]
List fast_ridit_scores_cpp(SEXP y_sexp, SEXP ref_idx_sexp) {
    IntegerVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> y(y_vec.begin(), y_vec.size());
    IntegerVector ref_idx_vec(ref_idx_sexp);
    Eigen::Map<const Eigen::VectorXi> ref_idx(ref_idx_vec.begin(), ref_idx_vec.size());
    int n = y.size();
    int n_ref = ref_idx.size();

    std::vector<int> ref_idx_zero_based;
    ref_idx_zero_based.reserve(static_cast<std::size_t>(n_ref));
    for (int i = 0; i < n_ref; ++i) {
        ref_idx_zero_based.push_back(ref_idx[i] - 1);
    }

    RiditScoreData ridit_data = fast_ridit_scores_from_ref_indices(y.data(), n, ref_idx_zero_based);
    return List::create(
        Named("scores") = ridit_data.scores,
        Named("levels") = wrap(ridit_data.levels),
        Named("ref_p") = ridit_data.ref_p
    );
}

// [[Rcpp::export]]
List fast_ridit_analysis_cpp(SEXP w_sexp, SEXP y_sexp, const std::string& reference = "control") {
    IntegerVector w_vec(w_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_vec.begin(), w_vec.size());
    IntegerVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> y(y_vec.begin(), y_vec.size());
    int n = y.size();
    std::vector<int> ref_idx;
    ref_idx.reserve(static_cast<std::size_t>(n));
    
    if (reference == "control") {
        for (int i = 0; i < n; ++i) if (w[i] == 0) ref_idx.push_back(i);
    } else if (reference == "treatment") {
        for (int i = 0; i < n; ++i) if (w[i] == 1) ref_idx.push_back(i);
    } else { // pooled
        for (int i = 0; i < n; ++i) ref_idx.push_back(i);
    }
    
    if (ref_idx.empty()) return List::create();
    
    RiditScoreData ridit_data = fast_ridit_scores_from_ref_indices(y.data(), n, ref_idx);
    NumericVector scores = ridit_data.scores;
    
    double sum_t = 0.0, sum_c = 0.0;
    int n_t = 0, n_c = 0;
    
    for (int i = 0; i < n; ++i) {
        if (w[i] == 1) {
            sum_t += scores[i];
            n_t++;
        } else {
            sum_c += scores[i];
            n_c++;
        }
    }
    
    double mean_ridit_t = (n_t > 0) ? sum_t / n_t : NA_REAL;
    double mean_ridit_c = (n_c > 0) ? sum_c / n_c : NA_REAL;
    
    // Variance of mean ridit (Bross 1958 approximation)
    // Var(W) = 1 / (12 * n_t) if reference is control and n_c is large
    // More generally, we can use the sample variance of the scores
    double var_t = 0.0;
    if (n_t > 1) {
        for (int i = 0; i < n; ++i) {
            if (w[i] == 1) {
                const double diff = scores[i] - mean_ridit_t;
                var_t += diff * diff;
            }
        }
        var_t /= (n_t - 1);
    }
    
    return List::create(
        Named("mean_ridit_t") = mean_ridit_t,
        Named("mean_ridit_c") = mean_ridit_c,
        Named("estimate") = mean_ridit_t - 0.5, // Centered at 0
        Named("se") = std::sqrt(var_t / n_t),
        Named("scores") = scores,
        Named("levels") = wrap(ridit_data.levels),
        Named("ref_p") = ridit_data.ref_p
    );
}
