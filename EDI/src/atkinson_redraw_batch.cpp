// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// Helper: Find which columns have variation in rows 0..(n-1)
static void which_cols_vary_subset(
	const MatrixXd& X,
	int n_rows,
	std::vector<int>& varying_cols
) {
	varying_cols.clear();
	int p = X.cols();

	for (int j = 0; j < p; ++j) {
		double first_val = X(0, j);
		for (int i = 1; i < n_rows; ++i) {
			if (X(i, j) != first_val) {
				varying_cols.push_back(j);
				break;
			}
		}
	}
}

// Helper: Find linearly independent columns via QR with column pivoting
template <typename Derived>
static void find_independent_cols(
	const Eigen::MatrixBase<Derived>& M,
	ColPivHouseholderQR<MatrixXd>& qr,
	std::vector<int>& indep_cols,
	double tol = 1e-12
) {
	indep_cols.clear();
	if (M.cols() == 0 || M.rows() == 0) {
		return;
	}

	qr.compute(M);
	qr.setThreshold(tol);
	int rank = qr.rank();

	const auto& pivot = qr.colsPermutation().indices();
	for (int i = 0; i < rank; ++i) {
		indep_cols.push_back(pivot(i));
	}
	std::sort(indep_cols.begin(), indep_cols.end());
}

// Helper: Compute single Atkinson weight for subject t (0-indexed)
// Returns the probability used for assignment (not the actual 0/1 assignment)
static double compute_atkinson_weight(
	const double* w_prev,
	int rows,
	const MatrixXd& X_prev_workspace,
	const VectorXd& xt_prev_workspace,
	int rank_prev,
	int t,
	FullPivLU<MatrixXd>& lu,
	MatrixXd& design_workspace,
	MatrixXd& crossprod_workspace,
	MatrixXd& inverse_workspace,
	VectorXd& xt_workspace
) {
	int cols = rank_prev + 2;

	if (rows == 0 || cols < 2) {
		return 0.5;  // Will use Bernoulli
	}

	// Build design matrix [w, 1, X_prev]
	auto XprevWT = design_workspace.topLeftCorner(rows, cols);
	for (int i = 0; i < rows; ++i) {
		XprevWT(i, 0) = w_prev[i];
		XprevWT(i, 1) = 1.0;
		for (int j = 0; j < rank_prev; ++j) {
			XprevWT(i, j + 2) = X_prev_workspace(i, j);
		}
	}

	// Compute (X'X)^-1
	auto XwtXw = crossprod_workspace.topLeftCorner(cols, cols);
	XwtXw.noalias() = XprevWT.transpose() * XprevWT;
	lu.compute(XwtXw);
	if (!lu.isInvertible()) {
		return 0.5;  // Will use Bernoulli
	}

	auto M = inverse_workspace.topLeftCorner(cols, cols);
	M.noalias() = static_cast<double>(t - 1) * lu.inverse();

	// Build xt vector with intercept
	auto xt_full = xt_workspace.head(rank_prev + 1);
	xt_full(0) = 1.0;
	for (int j = 0; j < rank_prev; ++j) {
		xt_full(j + 1) = xt_prev_workspace(j);
	}

	double A = M.row(0).segment(1, rank_prev + 1).dot(xt_full);
	if (A == 0 || !std::isfinite(A)) {
		return 0.5;  // Will use Bernoulli
	}

	double val = M(0, 0) / A + 1.0;
	double s_over_A_plus_one_sq = val * val;
	double prob = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1.0);
	prob = std::max(0.0, std::min(1.0, prob));

	return prob;
}

// [[Rcpp::export]]
NumericVector atkinson_redraw_batch_cpp(
	SEXP X_sexp,                   // Full covariate matrix (n x p), already numeric
	int n,                         // Number of subjects
	int p_raw,                     // Number of raw covariates (for early-subject Bernoulli threshold)
	double prob_T = 0.5           // Treatment probability for Bernoulli
) {
	NumericMatrix X_r(X_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	std::vector<double> results_vec(n);
    double* w_ptr = results_vec.data();
	const int p = X.cols();
	const int max_cols = p + 2;

	std::vector<int> var_cols;
	std::vector<int> indep_cols;
	var_cols.reserve(p);
	indep_cols.reserve(p);
	MatrixXd X_var_workspace(n, p);
	MatrixXd X_prev_workspace(n, p);
	VectorXd xt_prev_workspace(p);
	MatrixXd design_workspace(n, max_cols);
	MatrixXd crossprod_workspace(max_cols, max_cols);
	MatrixXd inverse_workspace(max_cols, max_cols);
	VectorXd xt_workspace(p + 1);
	ColPivHouseholderQR<MatrixXd> qr_workspace(n, p);
	FullPivLU<MatrixXd> lu_workspace(max_cols, max_cols);

	// Threshold for using Bernoulli (early subjects)
	int bernoulli_threshold = p_raw + 2 + 1;

	for (int t = 1; t <= n; ++t) {
		// For early subjects, use Bernoulli
		if (t <= bernoulli_threshold) {
			w_ptr[t - 1] = (R::unif_rand() < prob_T) ? 1.0 : 0.0;
			continue;
		}

		// Find which columns vary in rows 0..(t-1)
		which_cols_vary_subset(X, t, var_cols);

		if (var_cols.empty()) {
			w_ptr[t - 1] = (R::unif_rand() < prob_T) ? 1.0 : 0.0;
			continue;
		}

		// Extract submatrix for rows 0..(t-2) with varying columns
		auto X_var = X_var_workspace.topLeftCorner(t - 1, static_cast<int>(var_cols.size()));
		for (int i = 0; i < t - 1; ++i) {
			for (std::size_t j = 0; j < var_cols.size(); ++j) {
				X_var(i, static_cast<int>(j)) = X(i, var_cols[j]);
			}
		}

		// Find linearly independent columns
		find_independent_cols(X_var, qr_workspace, indep_cols);

		if (indep_cols.empty()) {
			w_ptr[t - 1] = (R::unif_rand() < prob_T) ? 1.0 : 0.0;
			continue;
		}

		// Extract final processed matrix and vector
		auto X_prev = X_prev_workspace.topLeftCorner(t - 1, static_cast<int>(indep_cols.size()));
		for (int i = 0; i < t - 1; ++i) {
			for (size_t j = 0; j < indep_cols.size(); ++j) {
				X_prev(i, static_cast<int>(j)) = X_var(i, indep_cols[j]);
			}
		}

		// Get xt_prev (current subject's row with same column processing)
		for (size_t j = 0; j < indep_cols.size(); ++j) {
			xt_prev_workspace(static_cast<int>(j)) = X(t - 1, var_cols[indep_cols[j]]);
		}

		int rank_prev = static_cast<int>(indep_cols.size());

		// Compute Atkinson probability
		double prob = compute_atkinson_weight(
			w_ptr, t - 1, X_prev_workspace, xt_prev_workspace, rank_prev, t,
			lu_workspace, design_workspace, crossprod_workspace,
			inverse_workspace, xt_workspace);

		// Draw assignment
		w_ptr[t - 1] = (R::unif_rand() < prob) ? 1.0 : 0.0;
	}

	return wrap(results_vec);
}
