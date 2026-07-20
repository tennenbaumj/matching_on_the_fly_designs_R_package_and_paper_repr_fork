#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

inline double brt_logit(double x, double clamp) {
	if (x < clamp) x = clamp;
	if (x > 1.0 - clamp) x = 1.0 - clamp;
	return std::log(x / (1.0 - x));
}

inline double brt_inv_logit(double x, double clamp) {
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

// Forward sharp-null shift matching R's shift_randomization_responses:
// 0 = additive ("none"); 1 = multiplicative ("log", survival/continuous);
// 2 = logit ("logit", proportion); 4 = multiplicative with count rounding ("log", count).
inline double brt_apply_shift(double y_val, double delta, int transform_code, double clamp) {
	if (transform_code == 1) return y_val * std::exp(delta);
	if (transform_code == 2) return brt_inv_logit(brt_logit(y_val, clamp) + delta, clamp);
	if (transform_code == 4) return std::round(y_val * std::exp(delta));
	return y_val + delta;
}

} // namespace

// Batch kernel for the bootstrap randomization test (BRT) null distribution of the
// simple mean difference. Each replicate b resamples rows i_mat(., b) of the sharp-null
// control outcomes y0 and pairs them with the fresh design assignment w_mat(., b);
// delta is the sharp-null shift applied to the treated, on the scale given by
// transform_code (see brt_apply_shift).
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int transform_code,
	double zero_one_logit_clamp,
	Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = w_mat.rows();

	std::vector<double> results_vec(nsim);

	const double* y0_ptr = y0.begin();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

	const bool has_noise = noise_mat.isNotNull();
	NumericMatrix noise_m;
	const double* noise_ptr = nullptr;
	if (has_noise) {
		noise_m = NumericMatrix(noise_mat);
		noise_ptr = noise_m.begin();
	}

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;
		double sum_T = 0, sum_C = 0;
		int n_T = 0;

		for (int i = 0; i < n; ++i) {
			double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
			if (has_noise) yv += noise_ptr[(size_t)b * n + i];
			const int is_t = (w_col[i] == 1);
			if (is_t) {
				sum_T += (delta != 0.0) ? brt_apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv;
				++n_T;
			} else {
				sum_C += yv;
			}
		}

		const int n_C = n - n_T;
		if (n_T == 0 || n_C == 0) {
			res_ptr[b] = NA_REAL;
		} else {
			res_ptr[b] = (sum_T / n_T) - (sum_C / n_C);
		}
	}

	return wrap(results_vec);
}

// Compatibility entry point for the generated Rcpp wrapper.  The R API still
// uses the original five-argument additive-scale interface.
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {
	return compute_rand_bootstrap_mean_diff_parallel_cpp(
		y0, i_mat, w_mat, delta, 0, NA_REAL, R_NilValue, num_cores
	);
}

// Batch kernel returning a 2-row matrix: row 0 = t0_b (mean diff), row 1 = se0_b (Welch SE).
// Used by the studentized/symmetric-percentile-t BRT to get both statistics in one C++ pass.
// [[Rcpp::export]]
NumericMatrix compute_rand_bootstrap_mean_diff_se_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int transform_code,
	double zero_one_logit_clamp,
	int num_cores) {

	int nsim = w_mat.cols();
	int n    = w_mat.rows();

	NumericMatrix result(2, nsim);

	const double* y0_ptr  = y0.begin();
	const int*    i_ptr   = i_mat.begin();
	const int*    w_ptr   = w_mat.begin();
	double* t0_ptr  = result.begin();          // row 0
	double* se_ptr  = result.begin() + nsim;   // row 1

	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;
		double sum_T = 0, sum_T2 = 0;
		double sum_C = 0, sum_C2 = 0;
		int n_T = 0;

		for (int i = 0; i < n; ++i) {
			const double yv_raw = y0_ptr[i_col[i] - 1];  // i_mat is 1-based
			const int is_t = (w_col[i] == 1);
			if (is_t) {
				const double yv = (delta != 0.0) ? brt_apply_shift(yv_raw, delta, transform_code, zero_one_logit_clamp) : yv_raw;
				sum_T  += yv;
				sum_T2 += yv * yv;
				++n_T;
			} else {
				sum_C  += yv_raw;
				sum_C2 += yv_raw * yv_raw;
			}
		}

		const int n_C = n - n_T;
		if (n_T < 2 || n_C < 2) {
			t0_ptr[b] = NA_REAL;
			se_ptr[b] = NA_REAL;
		} else {
			const double mean_T = sum_T / n_T;
			const double mean_C = sum_C / n_C;
			t0_ptr[b] = mean_T - mean_C;
			// Welch SE: sqrt(var_T/n_T + var_C/n_C), each variance via one-pass formula
			const double var_T = (sum_T2 - sum_T * mean_T) / (n_T - 1);
			const double var_C = (sum_C2 - sum_C * mean_C) / (n_C - 1);
			const double se_sq = var_T / n_T + var_C / n_C;
			se_ptr[b] = (se_sq > 0.0) ? std::sqrt(se_sq) : 0.0;
		}
	}

	return result;
}
