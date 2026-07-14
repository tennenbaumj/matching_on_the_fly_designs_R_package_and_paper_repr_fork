#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <map>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

constexpr size_t kExactMedianMaterializeLimit = 4096;

double median_in_place(std::vector<double>& values) {
    const size_t n = values.size();
    if (n == 0) {
        return NA_REAL;
    }

    const size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    const double upper = values[mid];

    if (n % 2 == 1) {
        return upper;
    }

    std::nth_element(values.begin(), values.begin() + mid - 1, values.begin() + mid);
    const double lower = values[mid - 1];
    return 0.5 * (lower + upper);
}

inline double logit_cpp(double x, double clamp) {
    if (x < clamp) x = clamp;
    if (x > 1.0 - clamp) x = 1.0 - clamp;
    return std::log(x / (1.0 - x));
}

inline double inv_logit_cpp(double x, double clamp) {
    double p;
    if (x >= 0.0) {
        const double z = std::exp(-x);
        p = 1.0 / (1.0 + z);
    } else {
        const double z = std::exp(x);
        p = z / (1.0 + z);
    }
    if (p < clamp) p = clamp;
    if (p > 1.0 - clamp) p = 1.0 - clamp;
    return p;
}

size_t count_pairwise_diffs_leq(const std::vector<double>& y_t,
                                const std::vector<double>& y_c,
                                double x) {
    size_t count = 0;
    size_t idx_t = 0;
    const size_t n_t = y_t.size();
    for (double yc : y_c) {
        const double limit = x + yc;
        while (idx_t < n_t && y_t[idx_t] <= limit) ++idx_t;
        count += idx_t;
    }
    return count;
}

double max_pairwise_diff_leq(const std::vector<double>& y_t,
                             const std::vector<double>& y_c,
                             double x) {
    double best = -std::numeric_limits<double>::infinity();
    bool found = false;
    size_t idx_t = 0;
    const size_t n_t = y_t.size();
    for (double yc : y_c) {
        const double limit = x + yc;
        while (idx_t < n_t && y_t[idx_t] <= limit) ++idx_t;
        if (idx_t > 0) {
            const double candidate = y_t[idx_t - 1] - yc;
            if (!found || candidate > best) best = candidate;
            found = true;
        }
    }
    return found ? best : NA_REAL;
}

double select_pairwise_diff_sorted(const std::vector<double>& y_t,
                                   const std::vector<double>& y_c,
                                   size_t rank) {
    const size_t total = y_t.size() * y_c.size();
    if (total == 0 || rank >= total) return NA_REAL;

    double lo = y_t.front() - y_c.back();
    double hi = y_t.back() - y_c.front();
    if (rank == 0) return lo;
    if (rank + 1 == total) return hi;

    const size_t target_count = rank + 1;
    for (int iter = 0; iter < 96; ++iter) {
        const double mid = lo + 0.5 * (hi - lo);
        if (mid == lo || mid == hi) break;
        if (count_pairwise_diffs_leq(y_t, y_c, mid) >= target_count) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double snapped = max_pairwise_diff_leq(y_t, y_c, hi);
    return std::isfinite(snapped) ? snapped : hi;
}

size_t count_walsh_avgs_leq(const std::vector<double>& d, double x) {
    const int m = static_cast<int>(d.size());
    int j = m - 1;
    size_t count = 0;
    const double limit = 2.0 * x;
    for (int i = 0; i < m; ++i) {
        while (j >= i && d[i] + d[j] > limit) --j;
        if (j < i) break;
        count += static_cast<size_t>(j - i + 1);
    }
    return count;
}

double max_walsh_avg_leq(const std::vector<double>& d, double x) {
    const int m = static_cast<int>(d.size());
    int j = m - 1;
    double best = -std::numeric_limits<double>::infinity();
    bool found = false;
    const double limit = 2.0 * x;
    for (int i = 0; i < m; ++i) {
        while (j >= i && d[i] + d[j] > limit) --j;
        if (j < i) break;
        const double candidate = 0.5 * (d[i] + d[j]);
        if (!found || candidate > best) best = candidate;
        found = true;
    }
    return found ? best : NA_REAL;
}

double select_walsh_avg_sorted(const std::vector<double>& d, size_t rank) {
    const size_t m = d.size();
    const size_t total = m * (m + 1) / 2;
    if (total == 0 || rank >= total) return NA_REAL;

    double lo = d.front();
    double hi = d.back();
    if (rank == 0) return lo;
    if (rank + 1 == total) return hi;

    const size_t target_count = rank + 1;
    for (int iter = 0; iter < 96; ++iter) {
        const double mid = lo + 0.5 * (hi - lo);
        if (mid == lo || mid == hi) break;
        if (count_walsh_avgs_leq(d, mid) >= target_count) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double snapped = max_walsh_avg_leq(d, hi);
    return std::isfinite(snapped) ? snapped : hi;
}

double hl_from_groups(std::vector<double> y_t, std::vector<double> y_c) {
    if (y_t.empty() || y_c.empty()) {
        return NA_REAL;
    }

    const size_t n_t = y_t.size();
    const size_t n_c = y_c.size();
    const size_t total = n_t * n_c;
    if (total <= kExactMedianMaterializeLimit) {
        std::vector<double> diffs(total);
        for (size_t i = 0; i < n_t; ++i) {
            double* p = diffs.data() + i * n_c;
            const double yt = y_t[i];
            for (size_t j = 0; j < n_c; ++j) p[j] = yt - y_c[j];
        }
        return median_in_place(diffs);
    }

    std::sort(y_t.begin(), y_t.end());
    std::sort(y_c.begin(), y_c.end());
    const size_t mid = total / 2;
    const double upper = select_pairwise_diff_sorted(y_t, y_c, mid);
    if (total % 2 == 1) return upper;
    const double lower = select_pairwise_diff_sorted(y_t, y_c, mid - 1);
    return 0.5 * (lower + upper);
}

double hl_signed_rank(std::vector<double> pair_diffs) {
    if (pair_diffs.empty()) {
        return NA_REAL;
    }

    const size_t m = pair_diffs.size();
    const size_t total = m * (m + 1) / 2;
    if (total <= kExactMedianMaterializeLimit) {
        std::vector<double> walsh_avgs(total);
        size_t k = 0;
        for (size_t i = 0; i < m; ++i) {
            const double half_di = 0.5 * pair_diffs[i];
            double* p = walsh_avgs.data() + k;
            const size_t len = m - i;
            for (size_t j = 0; j < len; ++j) p[j] = half_di + 0.5 * pair_diffs[i + j];
            k += len;
        }

        return median_in_place(walsh_avgs);
    }

    std::sort(pair_diffs.begin(), pair_diffs.end());
    const size_t mid = total / 2;
    const double upper = select_walsh_avg_sorted(pair_diffs, mid);
    if (total % 2 == 1) return upper;
    const double lower = select_walsh_avg_sorted(pair_diffs, mid - 1);
    return 0.5 * (lower + upper);
}

double estimate_hl_ssq_rank_sum(const std::vector<double>& y_t, const std::vector<double>& y_c) {
    if (y_t.size() < 2 || y_c.size() < 2) return NA_REAL;
    double sum = 0;
    double sum_sq = 0;
    const size_t n_c = y_c.size();
    for (size_t i = 0; i < y_t.size(); ++i) {
        const double yt = y_t[i];
        for (size_t j = 0; j < n_c; ++j) {
            const double d = yt - y_c[j];
            sum += d;
            sum_sq += d * d;
        }
    }
    const int count = static_cast<int>(y_t.size() * n_c);
    double var_diffs = (sum_sq - (sum * sum) / count) / (count - 1);
    return var_diffs / (y_t.size() + y_c.size());
}

double estimate_hl_ssq_signed_rank(const std::vector<double>& pair_diffs) {
    if (pair_diffs.size() < 2) return NA_REAL;
    double sum = 0;
    double sum_sq = 0;
    int count = 0;
    const size_t m = pair_diffs.size();
    for (size_t i = 0; i < m; ++i) {
        const double half_di = 0.5 * pair_diffs[i];
        for (size_t j = i; j < m; ++j) {
            const double a = half_di + 0.5 * pair_diffs[j];
            sum += a;
            sum_sq += a * a;
        }
        count += static_cast<int>(m - i);
    }
    double var_walsh = (sum_sq - (sum * sum) / count) / (count - 1);
    return var_walsh / m;
}

double apply_shift(double y_val, double delta, int transform_code, double zero_one_logit_clamp) {
    if (transform_code == 1) {
        return y_val * std::exp(delta);
    }
    if (transform_code == 2) {
        return inv_logit_cpp(logit_cpp(y_val, zero_one_logit_clamp) + delta, zero_one_logit_clamp);
    }
    if (transform_code == 3) {
        return (y_val + 1.0) * std::exp(delta) - 1.0;
    }
    if (transform_code == 4) {
        // count response: multiplicative shift with rounding, matching
        // shift_randomization_responses' as.integer(round(y * exp(delta)))
        return std::round(y_val * std::exp(delta));
    }
    return y_val + delta;
}

} // namespace

// [[Rcpp::export]]
double wilcox_hl_signed_rank_point_estimate_cpp(SEXP dy_sexp) {
	NumericVector dy_vec(dy_sexp);
    Eigen::Map<const Eigen::VectorXd> dy(dy_vec.begin(), dy_vec.size());
    std::vector<double> pair_diffs;
    pair_diffs.reserve(dy.size());
    const double* dy_ptr = dy.data();
    for (int i = 0; i < dy.size(); ++i) {
        if (std::isfinite(dy_ptr[i])) pair_diffs.push_back(dy_ptr[i]);
    }
    return hl_signed_rank(std::move(pair_diffs));
}

// [[Rcpp::export]]
double wilcox_hl_point_estimate_cpp(SEXP w_sexp, SEXP y_sexp) {
	IntegerVector w_int(w_sexp);
	NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_int.begin(), w_int.size());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    std::vector<double> y_t;
    std::vector<double> y_c;
    y_t.reserve(y.size());
    y_c.reserve(y.size());

    const double* y_ptr = y.data();
    const int* w_ptr = w.data();
    int n = y.size();

    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(y_ptr[i])) continue;
        if (w_ptr[i] == 1) y_t.push_back(y_ptr[i]);
        else if (w_ptr[i] == 0) y_c.push_back(y_ptr[i]);
    }

    return hl_from_groups(std::move(y_t), std::move(y_c));
}

// [[Rcpp::export]]
NumericVector compute_wilcox_hl_bootstrap_parallel_cpp(
    SEXP w_sexp,
    SEXP y_sexp,
    SEXP indices_mat_sexp,
    int num_cores) {

	IntegerVector w_int(w_sexp);
	NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_int.begin(), w_int.size());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	IntegerMatrix indices_int_mat(indices_mat_sexp);
    Eigen::Map<const Eigen::MatrixXi> indices_mat(indices_int_mat.begin(), indices_int_mat.nrow(), indices_int_mat.ncol());

    const int n = y.size();
    const int B = indices_mat.cols();
    
    std::vector<double> results_vec(B, NA_REAL);
    const double* y_ptr = y.data();
    const int* w_ptr = w.data();
    const int* idx_ptr = indices_mat.data();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        const int* idx_col = idx_ptr + (size_t)b * n;
        if (idx_col[0] < 0) {
            res_ptr[b] = NA_REAL;
            continue;
        }

        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const int idx = idx_col[i];
            if (idx < 0 || idx >= n) continue;
            
            const double y_val = y_ptr[idx];
            if (!std::isfinite(y_val)) continue;
            
            if (w_ptr[idx] == 1) y_t.push_back(y_val);
            else if (w_ptr[idx] == 0) y_c.push_back(y_val);
        }

        res_ptr[b] = hl_from_groups(std::move(y_t), std::move(y_c));
    }

    return wrap(results_vec);
}

//' Fast Wilcoxon HL Statistic for Multiple Permutations
//'
//' @param w_mat_sexp Integer matrix of permuted treatment assignments (n x r).
//' @param y_sexp Numeric response vector.
//' @param delta Null treatment effect shift.
//' @param transform_code Integer code for response transformation.
//' @param zero_one_logit_clamp Clamp value for logit transformation.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of HL statistics.
// [[Rcpp::export]]
NumericVector compute_wilcox_hl_distr_parallel_cpp(
    SEXP w_mat_sexp,
    SEXP y_sexp,
    double delta,
    int transform_code,
    double zero_one_logit_clamp,
    int num_cores) {

	IntegerMatrix w_int_mat(w_mat_sexp);
	NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::MatrixXi> w_mat(w_int_mat.begin(), w_int_mat.nrow(), w_int_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());

    const int n = y.size();
    const int nsim = w_mat.cols();
    
    std::vector<double> results_vec(nsim, NA_REAL);
    const double* y_ptr = y.data();
    const int* w_ptr = w_mat.data();
    double* res_ptr = results_vec.data();

    std::vector<double> y_shifted(n);
    if (delta != 0.0) {
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(y_ptr[i])) y_shifted[i] = apply_shift(y_ptr[i], delta, transform_code, zero_one_logit_clamp);
            else y_shifted[i] = NA_REAL;
        }
    }

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(y_ptr[i])) continue;

            if (w_col[i] == 1) {
                y_t.push_back(delta != 0.0 ? y_shifted[i] : y_ptr[i]);
            } else if (w_col[i] == 0) {
                y_c.push_back(y_ptr[i]);
            }
        }

        res_ptr[b] = hl_from_groups(std::move(y_t), std::move(y_c));
    }

    return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_wilcox_matching_ivwc_bootstrap_parallel_cpp(
    SEXP w_sexp,
    SEXP y_sexp,
    SEXP m_vec_sexp,
    SEXP indices_mat_sexp,
    SEXP m_mat_sexp,
    int num_cores) {

	IntegerVector w_int(w_sexp);
	NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_int.begin(), w_int.size());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	IntegerVector m_int_vec(m_vec_sexp);
    Eigen::Map<const Eigen::VectorXi> m_vec(m_int_vec.begin(), m_int_vec.size());
	IntegerMatrix indices_int_mat(indices_mat_sexp);
    Eigen::Map<const Eigen::MatrixXi> indices_mat(indices_int_mat.begin(), indices_int_mat.nrow(), indices_int_mat.ncol());
	IntegerMatrix m_int_mat(m_mat_sexp);
    Eigen::Map<const Eigen::MatrixXi> m_mat(m_int_mat.begin(), m_int_mat.nrow(), m_int_mat.ncol());

    int B = indices_mat.cols();
    int n = y.size();
    
    std::vector<double> results_vec(B, NA_REAL);
    const double* y_ptr = y.data();
    const int* w_ptr = w.data();
    const int* idx_ptr = indices_mat.data();
    const int* m_ptr = m_mat.data();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        const int* idx_col = idx_ptr + (size_t)b * n;
        const int* m_col = m_ptr + (size_t)b * n;
        
        std::vector<double> pair_diffs;
        std::vector<double> res_y_t;
        std::vector<double> res_y_c;

        std::map<int, std::vector<size_t>> pairs;
        for (int i = 0; i < n; ++i) {
            int mid = m_col[i];
            int idx = idx_col[i] - 1;
            if (idx < 0) continue;
            
            if (mid == 0) {
                if (w_ptr[idx] == 1) res_y_t.push_back(y_ptr[idx]);
                else res_y_c.push_back(y_ptr[idx]);
            } else {
                pairs[mid].push_back(idx);
            }
        }

        for (auto const& [mid, idxs] : pairs) {
            if (idxs.size() == 2) {
                int idx1 = idxs[0];
                int idx2 = idxs[1];
                if (w_ptr[idx1] == 1 && w_ptr[idx2] == 0) pair_diffs.push_back(y_ptr[idx1] - y_ptr[idx2]);
                else if (w_ptr[idx1] == 0 && w_ptr[idx2] == 1) pair_diffs.push_back(y_ptr[idx2] - y_ptr[idx1]);
            }
        }

        double beta_m = hl_signed_rank(pair_diffs);
        double ssq_m = estimate_hl_ssq_signed_rank(pair_diffs);
        double beta_r = hl_from_groups(res_y_t, res_y_c);
        double ssq_r = estimate_hl_ssq_rank_sum(res_y_t, res_y_c);

        bool m_ok = std::isfinite(beta_m) && std::isfinite(ssq_m) && ssq_m > 0;
        bool r_ok = std::isfinite(beta_r) && std::isfinite(ssq_r) && ssq_r > 0;

        if (m_ok && r_ok) {
            double w_star = ssq_r / (ssq_r + ssq_m);
            res_ptr[b] = w_star * beta_m + (1.0 - w_star) * beta_r;
        } else if (m_ok) {
            res_ptr[b] = beta_m;
        } else if (r_ok) {
            res_ptr[b] = beta_r;
        }
    }

    return wrap(results_vec);
}

// BRT variant: each replicate b resamples rows i_mat(., b) of the sharp-null outcomes y0
// and pairs them with the fresh assignment w_mat(., b); delta is the sharp-null shift
// applied to the freshly treated on the transform_code scale (see apply_shift; code 4 =
// count response, multiplicative with rounding).
// [[Rcpp::export]]
NumericVector compute_wilcox_hl_rand_bootstrap_parallel_cpp(
    SEXP y0_sexp,
    SEXP i_mat_sexp,
    SEXP w_mat_sexp,
    double delta,
    int transform_code,
    double zero_one_logit_clamp,
    int num_cores) {

	NumericVector y0_vec(y0_sexp);
	IntegerMatrix i_int_mat(i_mat_sexp);
	IntegerMatrix w_int_mat(w_mat_sexp);
    Eigen::Map<const Eigen::VectorXd> y0(y0_vec.begin(), y0_vec.size());
    Eigen::Map<const Eigen::MatrixXi> i_mat(i_int_mat.begin(), i_int_mat.nrow(), i_int_mat.ncol());
    Eigen::Map<const Eigen::MatrixXi> w_mat(w_int_mat.begin(), w_int_mat.nrow(), w_int_mat.ncol());

    const int n = i_mat.rows();
    const int nsim = i_mat.cols();
    std::vector<double> results_vec(nsim, NA_REAL);
    const double* y0_ptr = y0.data();
    const int* i_ptr = i_mat.data();
    const int* w_ptr = w_mat.data();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        const int* i_col = i_ptr + (size_t)b * n;
        const int* w_col = w_ptr + (size_t)b * n;
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
            if (!std::isfinite(yv)) continue;
            if (w_col[i] == 1) {
                y_t.push_back(delta != 0.0 ? apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv);
            } else if (w_col[i] == 0) {
                y_c.push_back(yv);
            }
        }

        res_ptr[b] = hl_from_groups(std::move(y_t), std::move(y_c));
    }

    return wrap(results_vec);
}
