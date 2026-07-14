#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Batch kernel for the BRT null distribution of the OLS treatment coefficient.
// Each replicate b builds the design matrix [1, w_fresh, Xc[i_b, ]] on the resampled
// rows, adds the additive sharp-null shift delta to the freshly treated responses,
// and returns the treatment coefficient (column 2) from the least-squares fit.
// Only the additive shift (transform "none") is supported — the OLS classes are
// continuous-response.
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_ols_parallel_cpp(
	const NumericVector& y0,
	const NumericMatrix& Xc,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {

	const int nsim = w_mat.cols();
	const int n = w_mat.rows();
	const int p_cov = Xc.ncol();
	const int p = 2 + p_cov;

	std::vector<double> results_vec(nsim, NA_REAL);

	const double* y0_ptr = y0.begin();
	const double* xc_ptr = (p_cov > 0) ? Xc.begin() : nullptr;
	const int n_full = y0.size();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n * p * p, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;

		Eigen::MatrixXd M(n, p);
		Eigen::VectorXd yb(n);
		int n_T = 0;
		for (int i = 0; i < n; ++i) {
			const int row0 = i_col[i] - 1; // i_mat is 1-based
			const int is_t = (w_col[i] == 1);
			M(i, 0) = 1.0;
			M(i, 1) = static_cast<double>(is_t);
			for (int j = 0; j < p_cov; ++j) {
				M(i, 2 + j) = xc_ptr[(size_t)j * n_full + row0];
			}
			yb(i) = y0_ptr[row0] + (is_t ? delta : 0.0);
			n_T += is_t;
		}
		if (n_T == 0 || n_T == n || n <= p) continue; // leaves NA_REAL

		const Eigen::MatrixXd MtM = M.transpose() * M;
		const Eigen::VectorXd Mty = M.transpose() * yb;
		Eigen::LDLT<Eigen::MatrixXd> ldlt(MtM);
		if (ldlt.info() != Eigen::Success) continue;
		const Eigen::VectorXd beta = ldlt.solve(Mty);
		// guard against a numerically singular design (LDLT does not fail on rank deficiency)
		const double resid_norm = (MtM * beta - Mty).norm();
		if (!std::isfinite(beta(1)) || resid_norm > 1e-6 * (1.0 + Mty.norm())) continue;
		res_ptr[b] = beta(1);
	}

	return wrap(results_vec);
}
