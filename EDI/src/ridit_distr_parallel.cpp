#include "_helper_functions.h"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

// Helper to compute ridit scores and levels from a reference set
void get_ridit_map_cpp(const std::vector<int>& y_ref,
                  std::vector<int>& levels,
                  std::vector<double>& ridit_scores,
                  std::vector<int>& counts) {
    int n_ref = y_ref.size();
    if (n_ref == 0) return;

    levels = y_ref;
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

    ridit_scores.assign(levels.size(), 0.0);
    counts.assign(levels.size(), 0);
    for (int val : y_ref) {
        auto it = std::lower_bound(levels.begin(), levels.end(), val);
        counts[static_cast<std::size_t>(it - levels.begin())]++;
    }

    double cumulative_p = 0.0;
    for (std::size_t i = 0; i < levels.size(); ++i) {
        double p_k = static_cast<double>(counts[i]) / n_ref;
        ridit_scores[i] = cumulative_p + 0.5 * p_k;
        cumulative_p += p_k;
    }
}

double compute_mean_ridit_with_map_cpp(const std::vector<int>& y_target,
                                  const std::vector<double>& ridit_scores,
                                  const std::vector<int>& levels) {
    if (y_target.empty() || ridit_scores.empty()) return NA_REAL;

    double sum_t = 0.0;
    for (int val : y_target) {
        auto it = std::lower_bound(levels.begin(), levels.end(), val);
        if (it != levels.end() && *it == val) {
            sum_t += ridit_scores[static_cast<std::size_t>(it - levels.begin())];
        } else {
            if (it == levels.begin()) {
                // sum_t += 0.0
            } else if (it == levels.end()) {
                sum_t += 1.0;
            } else {
                int idx = static_cast<int>(std::distance(levels.begin(), it));
                sum_t += (ridit_scores[idx - 1] + ridit_scores[idx]) / 2.0;
            }
        }
    }
    return sum_t / y_target.size();
}

double compute_single_ridit_estimate_cpp(const int* y_b,
                                   const int* w_b,
                                   int n,
                                   const std::string& reference,
                                   std::vector<int>& y_ref,
                                   std::vector<int>& y_t,
                                   std::vector<int>& levels,
                                   std::vector<double>& ridit_scores,
                                   std::vector<int>& counts) {
    y_ref.clear();
    y_t.clear();

    for (int i = 0; i < n; ++i) {
        if (w_b[i] == 1) y_t.push_back(y_b[i]);
        
        if (reference == "control") {
            if (w_b[i] == 0) y_ref.push_back(y_b[i]);
        } else if (reference == "treatment") {
            if (w_b[i] == 1) y_ref.push_back(y_b[i]);
        } else { // pooled
            y_ref.push_back(y_b[i]);
        }
    }
    
    if (y_ref.empty() || y_t.empty()) return NA_REAL;
    
    get_ridit_map_cpp(y_ref, levels, ridit_scores, counts);
    
    return compute_mean_ridit_with_map_cpp(y_t, ridit_scores, levels) - 0.5;
}

} // namespace

// [[Rcpp::export]]
NumericVector compute_ridit_distr_parallel_cpp(const IntegerVector& y,
                                             const IntegerMatrix& w_mat, 
                                             std::string reference, 
                                             int num_cores) {
    int nsim = w_mat.cols();
    int n = y.size();
    std::vector<double> results_vec(nsim, NA_REAL);
    double* res_ptr = results_vec.data();
    
    const int* y_ptr = y.begin();
    const int* w_mat_ptr = w_mat.begin();

    std::vector<int> y_std(n);
    for(int i=0; i<n; ++i) y_std[i] = y_ptr[i];

    bool is_pooled = (reference == "pooled");
    std::vector<int> global_levels;
    std::vector<double> global_ridit_scores;
    std::vector<int> global_counts;
    if (is_pooled) {
        get_ridit_map_cpp(y_std, global_levels, global_ridit_scores, global_counts);
    }

    const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);
#ifdef _OPENMP
    if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
    {
        std::vector<int> y_ref;
        std::vector<int> y_t;
        std::vector<int> levels;
        std::vector<double> ridit_scores;
        std::vector<int> counts;
        y_ref.reserve(n);
        y_t.reserve(n);
        levels.reserve(n);
        ridit_scores.reserve(n);
        counts.reserve(n);

        #pragma omp for schedule(static)
        for (int b = 0; b < nsim; ++b) {
            const int* w_col = w_mat_ptr + (size_t)b * n;

            if (is_pooled) {
                y_t.clear();
                y_t.reserve(n);
                for (int i = 0; i < n; ++i) {
                    if (w_col[i] == 1) y_t.push_back(y_std[i]);
                }
                res_ptr[b] = y_t.empty() ? NA_REAL :
                    compute_mean_ridit_with_map_cpp(y_t, global_ridit_scores, global_levels) - 0.5;
            } else {
                res_ptr[b] = compute_single_ridit_estimate_cpp(
                    y_std.data(), w_col, n, reference, y_ref, y_t, levels, ridit_scores, counts
                );
            }
        }
    }
    
    return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_ridit_bootstrap_parallel_cpp(const IntegerVector& w,
                                                 const IntegerVector& y, 
                                                 const IntegerMatrix& indices_mat, 
                                                 std::string reference, 
                                                 int num_cores) {
    int B = indices_mat.ncol();
    int n = y.size();
    std::vector<double> results_vec(B, NA_REAL);
    double* res_ptr = results_vec.data();

    const int* y_ptr = y.begin();
    const int* w_ptr = w.begin();
    const int* idx_mat_ptr = indices_mat.begin();

    const bool use_parallel = should_parallelize_replicates(B, n, num_cores);
#ifdef _OPENMP
    if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
    {
        std::vector<int> y_b(n);
        std::vector<int> w_b(n);
        std::vector<int> y_ref;
        std::vector<int> y_t;
        std::vector<int> levels;
        std::vector<double> ridit_scores;
        std::vector<int> counts;
        y_ref.reserve(n);
        y_t.reserve(n);
        levels.reserve(n);
        ridit_scores.reserve(n);
        counts.reserve(n);

        #pragma omp for schedule(static)
        for (int b = 0; b < B; ++b) {
            const int* idx_col = idx_mat_ptr + (size_t)b * n;
            if (idx_col[0] == -1) {
                res_ptr[b] = NA_REAL;
                continue;
            }

            for (int i = 0; i < n; ++i) {
                int idx = idx_col[i] - 1; // 1-indexed to 0-indexed
                if (idx < 0 || idx >= n) {
                    y_b[i] = 0;
                    w_b[i] = 0;
                } else {
                    y_b[i] = y_ptr[idx];
                    w_b[i] = w_ptr[idx];
                }
            }
            res_ptr[b] = compute_single_ridit_estimate_cpp(
                y_b.data(), w_b.data(), n, reference, y_ref, y_t, levels, ridit_scores, counts
            );
        }
    }
    
    return wrap(results_vec);
}

// BRT variant of the ridit estimate: each replicate b resamples rows i_mat(., b) and
// pairs them with the fresh assignment w_mat(., b); the ridit reference map is rebuilt
// from the resampled responses per replicate (matching what a sub-inference on the
// bootstrap sample computes). No sharp-null shift (ordinal; the R hook declines
// delta != 0).
// [[Rcpp::export]]
NumericVector compute_ridit_rand_bootstrap_parallel_cpp(const IntegerVector& y0,
                                                        const IntegerMatrix& i_mat,
                                                        const IntegerMatrix& w_mat,
                                                        std::string reference,
                                                        int num_cores) {
    const int n = i_mat.nrow();
    const int nsim = i_mat.ncol();
    std::vector<double> results_vec(nsim, NA_REAL);
    double* res_ptr = results_vec.data();
    const int* y0_ptr = y0.begin();
    const int* i_ptr = i_mat.begin();
    const int* w_ptr = w_mat.begin();
    const bool is_pooled = (reference == "pooled");

    const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);
#ifdef _OPENMP
    if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
    {
        std::vector<int> y_b(n);
        std::vector<int> y_ref;
        std::vector<int> y_t;
        std::vector<int> levels;
        std::vector<double> ridit_scores;
        std::vector<int> counts;
        y_ref.reserve(n);
        y_t.reserve(n);
        levels.reserve(n);
        ridit_scores.reserve(n);
        counts.reserve(n);

        #pragma omp for schedule(static)
        for (int b = 0; b < nsim; ++b) {
            const int* i_col = i_ptr + (size_t)b * n;
            const int* w_col = w_ptr + (size_t)b * n;
            for (int i = 0; i < n; ++i) y_b[i] = y0_ptr[i_col[i] - 1]; // i_mat is 1-based

            if (is_pooled) {
                get_ridit_map_cpp(y_b, levels, ridit_scores, counts);
                y_t.clear();
                for (int i = 0; i < n; ++i) {
                    if (w_col[i] == 1) y_t.push_back(y_b[i]);
                }
                res_ptr[b] = y_t.empty() ? NA_REAL :
                    compute_mean_ridit_with_map_cpp(y_t, ridit_scores, levels) - 0.5;
            } else {
                res_ptr[b] = compute_single_ridit_estimate_cpp(
                    y_b.data(), w_col, n, reference, y_ref, y_t, levels, ridit_scores, counts
                );
            }
        }
    }

    return wrap(results_vec);
}
