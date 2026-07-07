#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_matching_compound_distr_parallel_cpp(
	SEXP y_sexp,
	SEXP w_mat_sexp,
	SEXP m_mat_sexp,
	int num_cores) {
	NumericVector y_r(y_sexp);
	IntegerMatrix w_mat_r(w_mat_sexp);
	IntegerMatrix m_mat_r(m_mat_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::MatrixXi> w_mat(w_mat_r.begin(), w_mat_r.nrow(), w_mat_r.ncol());
	Eigen::Map<const Eigen::MatrixXi> m_mat(m_mat_r.begin(), m_mat_r.nrow(), m_mat_r.ncol());

	int nsim = w_mat.cols();
	int n = y.size();
	std::vector<double> results_vec(nsim, NA_REAL);
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
	{
		std::vector<double> diffs;
		std::vector<int> treated_idx;
		std::vector<int> control_idx;

		#pragma omp for schedule(static)
		for (int b = 0; b < nsim; ++b) {
			const int* w_col = w_mat.data() + static_cast<size_t>(b) * n;
			const int* m_col = m_mat.data() + static_cast<size_t>(b) * n;

			// Single merged pass: max_m + scatter fill + unmatched sum+sumsq.
			// Reads m_col, w_col, y once each instead of 3-4 times.
			treated_idx.assign(n, -1);
			control_idx.assign(n, -1);
			int m = 0;
			int nRT = 0, nRC = 0;
			double sum_T = 0, sum_C = 0, sum_T2 = 0, sum_C2 = 0;
			for (int i = 0; i < n; ++i) {
				const int match_id = m_col[i];
				if (match_id > m) m = match_id;
				if (match_id > 0) {
					const int slot = match_id - 1;
					if (w_col[i] == 1) treated_idx[slot] = i;
					else control_idx[slot] = i;
				} else {
					const double yi = y[i];
					if (w_col[i] == 1) { nRT++; sum_T += yi; sum_T2 += yi * yi; }
					else { nRC++; sum_C += yi; sum_C2 += yi * yi; }
				}
			}

			// Compute matched-pair stats in one O(m) pass.
			double d_bar = NA_REAL;
			double ssqD_bar = NA_REAL;
			if (m > 0) {
				double sum_d = 0.0, sum_d2 = 0.0;
				diffs.resize(m);
				for (int slot = 0; slot < m; ++slot) {
					const double diff = y[treated_idx[slot]] - y[control_idx[slot]];
					diffs[slot] = diff;
					sum_d  += diff;
					sum_d2 += diff * diff;
				}
				d_bar = sum_d / m;
				if (m > 1)
					ssqD_bar = (sum_d2 - m * d_bar * d_bar) / (m - 1) / m;
			}

			// Unmatched stats: variance via sum-of-squares (no second pass).
			double r_bar = NA_REAL;
			double ssqR = NA_REAL;
			if (nRT > 0 && nRC > 0) {
				double mean_T = sum_T / nRT;
				double mean_C = sum_C / nRC;
				r_bar = mean_T - mean_C;
				if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
					double var_T = (sum_T2 - nRT * mean_T * mean_T) / (nRT - 1);
					double var_C = (sum_C2 - nRC * mean_C * mean_C) / (nRC - 1);
					int nR = nRT + nRC;
					ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
				}
			}

			double beta_hat_T = NA_REAL;
			if (nRT <= 1 || nRC <= 1) {
				beta_hat_T = d_bar;
			} else if (m == 0) {
				beta_hat_T = r_bar;
			} else {
				if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) beta_hat_T = r_bar;
				else if (!std::isfinite(ssqR) || ssqR <= 0) beta_hat_T = d_bar;
				else {
					double w_star = ssqR / (ssqR + ssqD_bar);
					beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
				}
			}
			res_ptr[b] = beta_hat_T;
		}
	}
	return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_matching_compound_bootstrap_parallel_cpp(
	SEXP w_mat_sexp,
	SEXP m_mat_sexp,
	SEXP y_mat_sexp,
	int num_cores) {
	IntegerMatrix w_mat_r(w_mat_sexp);
	IntegerMatrix m_mat_r(m_mat_sexp);
	NumericMatrix y_mat_r(y_mat_sexp);
	Eigen::Map<const Eigen::MatrixXi> w_mat(w_mat_r.begin(), w_mat_r.nrow(), w_mat_r.ncol());
	Eigen::Map<const Eigen::MatrixXi> m_mat(m_mat_r.begin(), m_mat_r.nrow(), m_mat_r.ncol());
	Eigen::Map<const Eigen::MatrixXd> y_mat(y_mat_r.begin(), y_mat_r.nrow(), y_mat_r.ncol());

	int nsim = w_mat.cols();
	int n = w_mat.rows();
	std::vector<double> results_vec(nsim, NA_REAL);
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
	{
		std::vector<double> diffs;
		std::vector<int> treated_idx;
		std::vector<int> control_idx;

		#pragma omp for schedule(static)
		for (int b = 0; b < nsim; ++b) {
			const double* y_col = y_mat.data() + static_cast<size_t>(b) * n;
			const int* w_col = w_mat.data() + static_cast<size_t>(b) * n;
			const int* m_col = m_mat.data() + static_cast<size_t>(b) * n;

			// Single merged pass: max_m + scatter fill + unmatched sum+sumsq.
			treated_idx.assign(n, -1);
			control_idx.assign(n, -1);
			int m = 0;
			int nRT = 0, nRC = 0;
			double sum_T = 0, sum_C = 0, sum_T2 = 0, sum_C2 = 0;
			for (int i = 0; i < n; ++i) {
				const int match_id = m_col[i];
				if (match_id > m) m = match_id;
				if (match_id > 0) {
					const int slot = match_id - 1;
					if (w_col[i] == 1) treated_idx[slot] = i;
					else control_idx[slot] = i;
				} else {
					const double yi = y_col[i];
					if (w_col[i] == 1) { nRT++; sum_T += yi; sum_T2 += yi * yi; }
					else { nRC++; sum_C += yi; sum_C2 += yi * yi; }
				}
			}

			double d_bar = NA_REAL;
			double ssqD_bar = NA_REAL;
			if (m > 0) {
				double sum_d = 0.0, sum_d2 = 0.0;
				diffs.resize(m);
				for (int slot = 0; slot < m; ++slot) {
					const double diff = y_col[treated_idx[slot]] - y_col[control_idx[slot]];
					diffs[slot] = diff;
					sum_d  += diff;
					sum_d2 += diff * diff;
				}
				d_bar = sum_d / m;
				if (m > 1)
					ssqD_bar = (sum_d2 - m * d_bar * d_bar) / (m - 1) / m;
			}

			double r_bar = NA_REAL;
			double ssqR = NA_REAL;
			if (nRT > 0 && nRC > 0) {
				double mean_T = sum_T / nRT;
				double mean_C = sum_C / nRC;
				r_bar = mean_T - mean_C;
				if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
					double var_T = (sum_T2 - nRT * mean_T * mean_T) / (nRT - 1);
					double var_C = (sum_C2 - nRC * mean_C * mean_C) / (nRC - 1);
					int nR = nRT + nRC;
					ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
				}
			}

			double beta_hat_T = NA_REAL;
			if (nRT <= 1 || nRC <= 1) {
				beta_hat_T = d_bar;
			} else if (m == 0) {
				beta_hat_T = r_bar;
			} else {
				if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) beta_hat_T = r_bar;
				else if (!std::isfinite(ssqR) || ssqR <= 0) beta_hat_T = d_bar;
				else {
					double w_star = ssqR / (ssqR + ssqD_bar);
					beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
				}
			}
			res_ptr[b] = beta_hat_T;
		}
	}
	return wrap(results_vec);
}
