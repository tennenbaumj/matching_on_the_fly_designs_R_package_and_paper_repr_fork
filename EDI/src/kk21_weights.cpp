#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace Eigen;

inline Eigen::ArrayXd clamp_weibull_z_kk(const Eigen::ArrayXd& z) {
	return z.max(-20.0).min(20.0);
}

inline Eigen::ArrayXd plogis_array_kk(const Eigen::ArrayXd& eta) {
	const Eigen::Array<bool, Eigen::Dynamic, 1> nonnegative = (eta >= 0.0);
	const Eigen::ArrayXd pos = 1.0 / (1.0 + (-eta).exp());
	const Eigen::ArrayXd neg_exp = eta.exp();
	const Eigen::ArrayXd neg = neg_exp / (1.0 + neg_exp);
	return nonnegative.select(pos, neg);
}

// [[Rcpp::export]]
NumericVector kk21_continuous_weights_cpp(const NumericMatrix& X,
											const NumericVector& y) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p);

	if (n < 2 || p == 0) {
	std::fill(weights.begin(), weights.end(), std::numeric_limits<double>::epsilon());
	return weights;
	}

	double sumy = 0.0;
	double sumy2 = 0.0;
	for (int i = 0; i < n; ++i) {
	double yi = y[i];
	sumy += yi;
	sumy2 += yi * yi;
	}
	double ybar = sumy / static_cast<double>(n);
	double eps = std::numeric_limits<double>::epsilon();

	for (int j = 0; j < p; ++j) {
	double sumx = 0.0;
	double sumx2 = 0.0;
	double sumxy = 0.0;
	for (int i = 0; i < n; ++i) {
		double x = X(i, j);
		sumx += x;
		sumx2 += x * x;
		sumxy += x * y[i];
	}
	double xbar = sumx / static_cast<double>(n);
	double varx = sumx2 - static_cast<double>(n) * xbar * xbar;

	if (varx <= eps || n <= 2) {
		weights[j] = eps;
		continue;
	}

	double b1 = (sumxy - static_cast<double>(n) * xbar * ybar) / varx;
	double b0 = ybar - b1 * xbar;
	double sse = sumy2 + static_cast<double>(n) * b0 * b0 + b1 * b1 * sumx2 +
		2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;
	double denom = static_cast<double>(n - 2);
	if (denom <= 0.0 || sse <= 0.0) {
		weights[j] = eps;
		continue;
	}
	double sigma2 = sse / denom;
	double se_b1 = std::sqrt(sigma2 / varx);
	if (!R_finite(se_b1) || se_b1 <= 0.0) {
		weights[j] = eps;
	} else {
		weights[j] = std::fabs(b1 / se_b1);
	}
	}

	return weights;
}

// [[Rcpp::export]]
NumericVector kk21_logistic_weights_cpp(const NumericMatrix& X,
										const NumericVector& y,
										const int maxit = 100,
										const double tol = 1e-8) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p);
	double eps = std::numeric_limits<double>::epsilon();

	if (n == 0 || p == 0) {
	std::fill(weights.begin(), weights.end(), eps);
	return weights;
	}

	const Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	const Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

	for (int j = 0; j < p; ++j) {
	double b0 = 0.0;
	double b1 = 0.0;
	bool singular = false;
	const Eigen::ArrayXd x = X_map.col(j).array();

	for (int iter = 0; iter < maxit; ++iter) {
		const Eigen::ArrayXd eta = b0 + b1 * x;
		const Eigen::ArrayXd p_i = plogis_array_kk(eta);
		const Eigen::ArrayXd w = (p_i * (1.0 - p_i)).max(1e-10);
		const Eigen::ArrayXd z = eta + (y_vec.array() - p_i) / w;

		double S0 = w.sum();
		double S1 = (w * x).sum();
		double S2 = (w * x.square()).sum();
		double Sz0 = (w * z).sum();
		double Sz1 = (w * x * z).sum();

		double det = S0 * S2 - S1 * S1;
		if (det <= eps) {
		singular = true;
		break;
		}

		double b0_new = (Sz0 * S2 - Sz1 * S1) / det;
		double b1_new = (S0 * Sz1 - S1 * Sz0) / det;

		if (std::fabs(b0_new - b0) < tol && std::fabs(b1_new - b1) < tol) {
		b0 = b0_new;
		b1 = b1_new;
		break;
		}
		b0 = b0_new;
		b1 = b1_new;
	}

	if (singular) {
		weights[j] = eps;
		continue;
	}

	const Eigen::ArrayXd eta = b0 + b1 * x;
	const Eigen::ArrayXd p_i = plogis_array_kk(eta);
	const Eigen::ArrayXd w = (p_i * (1.0 - p_i)).max(1e-10);
	double S0 = w.sum();
	double S1 = (w * x).sum();
	double S2 = (w * x.square()).sum();

	double det = S0 * S2 - S1 * S1;
	if (det <= eps || S0 <= 0.0) {
		weights[j] = eps;
		continue;
	}
	double var_b1 = S0 / det;
	double se_b1 = std::sqrt(var_b1);
	if (!R_finite(se_b1) || se_b1 <= 0.0) {
		weights[j] = eps;
	} else {
		weights[j] = std::fabs(b1 / se_b1);
	}
	}

	return weights;
}

// Optimized using Frisch-Waugh-Lovell (FWL) theorem and incremental residual updates.
// We maintain E_x (residuals of candidates) and e_y (residuals of outcome) orthogonal
// to the selected basis Q. At each step, selecting best_j allows updating all
// remaining E_x columns in O(n*p) total, leading to O(n*p^2) overall.
// [[Rcpp::export]]
NumericVector kk21_stepwise_continuous_weights_cpp(const NumericMatrix& X,
													 const NumericVector& y,
													 const NumericVector& w) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p, NA_REAL);

	if (n == 0 || p == 0) {
	return weights;
	}

	const double eps = std::numeric_limits<double>::epsilon();

	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
	Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

	// Basis Q starts with intercept and treatment indicator w
	Eigen::MatrixXd Q(n, p + 2);
	int m = 0;

	// Intercept
	Q.col(0).setConstant(1.0 / std::sqrt(static_cast<double>(n)));
	m = 1;

	// Add w, orthogonalizing against intercept
	{
	Eigen::VectorXd u = w_vec - Q.col(0).dot(w_vec) * Q.col(0);
	double nm = u.norm();
	if (nm > 1e-8) {
		Q.col(1) = u / nm;
		m = 2;
	}
	}

	// Initial residuals: outcome e_y and ALL candidates E_x
	Eigen::VectorXd e_y = y_vec - Q.leftCols(m) * (Q.leftCols(m).transpose() * y_vec);
	double c = e_y.squaredNorm();

	Eigen::MatrixXd E_x = X_map;
	for (int j = 0; j < p; ++j) {
		E_x.col(j) -= Q.leftCols(m) * (Q.leftCols(m).transpose() * E_x.col(j));
	}

	std::vector<bool> used(p, false);

	for (int step = 0; step < p; ++step) {
	int df = n - m - 1;
	if (df <= 0) break;

	double best_stat = -1.0;
	int best_j = -1;
	double best_b = 0.0;
	double best_a = 0.0;

	for (int j = 0; j < p; ++j) {
		if (used[j]) continue;

		double b = E_x.col(j).squaredNorm();
		if (b < 1e-10) continue; // x_j collinear

		double a = E_x.col(j).dot(e_y);
		double sse = c - a * a / b;
		if (sse < 0.0) sse = 0.0;
		double sigma2 = sse / static_cast<double>(df);
		if (sigma2 <= 0.0 || !std::isfinite(sigma2)) continue;

		double stat = std::fabs(a) / std::sqrt(b * sigma2);
		if (!R_finite(stat)) stat = 0.0;

		if (stat > best_stat) {
			best_stat = stat;
			best_j = j;
			best_b = b;
			best_a = a;
		}
	}

	if (best_j < 0) break;

	weights[best_j] = best_stat;
	used[best_j] = true;

	// Update Q with best candidate: q_new = e_j / ||e_j||
	if (m < p + 2) {
		Q.col(m) = E_x.col(best_j) / std::sqrt(best_b);
		const Eigen::VectorXd q_new = Q.col(m);

		// Update e_y: e_y -= (q_new' e_y) * q_new
		double proj_y = q_new.dot(e_y);
		e_y -= proj_y * q_new;
		c -= proj_y * proj_y;
		if (c < 0.0) c = 0.0;

		// Update ALL remaining candidate residuals E_x: E_x_j -= (q_new' E_x_j) * q_new
		for (int j = 0; j < p; ++j) {
			if (!used[j]) {
				E_x.col(j) -= q_new.dot(E_x.col(j)) * q_new;
			}
		}
		m++;
	}
	}

	return weights;
}

// Helper: fit logistic regression on reduced design matrix XS via IRLS.
// beta is used as warm warm_start_params and updated in place. Returns false on failure.
static bool logistic_reduced_fit_for_score_test(
	const Eigen::MatrixXd& XS,
	const Eigen::VectorXd& y,
	Eigen::VectorXd& beta,
	Eigen::VectorXd& p_hat,
	Eigen::VectorXd& W_diag,
	Eigen::VectorXd& resid,
	Eigen::LDLT<Eigen::MatrixXd>& XtWX_ldlt,
	int maxit = 100,
	double tol = 1e-8) {
	int n = XS.rows();
	int m = XS.cols();
	if (n <= m) return false;

	p_hat.resize(n);
	W_diag.resize(n);

	for (int iter = 0; iter < maxit; ++iter) {
	Eigen::VectorXd eta = XS * beta;
	p_hat = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
	W_diag = (p_hat.array() * (1.0 - p_hat.array())).max(1e-10).matrix();
	Eigen::VectorXd z = eta + (y - p_hat).cwiseQuotient(W_diag);
	Eigen::MatrixXd XtWX = weighted_crossprod(XS, W_diag);
	Eigen::VectorXd beta_new = XtWX.ldlt().solve(weighted_crossprod_rhs(XS, W_diag, z));
	if ((beta - beta_new).norm() < tol) {
		beta = beta_new;
		break;
	}
	beta = beta_new;
	}

	// Final quantities for score test
	Eigen::VectorXd eta = XS * beta;
	p_hat = (1.0 / (1.0 + (-eta.array()).exp())).matrix();
	W_diag = (p_hat.array() * (1.0 - p_hat.array())).max(1e-10).matrix();
	resid = y - p_hat;
	Eigen::MatrixXd XtWX = weighted_crossprod(XS, W_diag);
	XtWX_ldlt = XtWX.ldlt();
	return XtWX_ldlt.info() == Eigen::Success;
}

// Optimized using score test (LM test) from the reduced model.
// At each step k, instead of fitting (p-k) logistic regressions of size n x (k+3), we:
//   1. Fit ONE logistic regression on the reduced model [1, X_sel, w] (with warm warm_start_params)
//   2. For each candidate x_j compute the score test statistic:
//      stat_j = |x_j' (y - p_hat)| / sqrt(I_jj.S)
//      where I_jj.S = x_j' W x_j - (X_S' W x_j)' (X_S' W X_S)^{-1} (X_S' W x_j)
// Asymptotically equivalent to the Wald t-statistic; gives the same ordering in practice.
// [[Rcpp::export]]
NumericVector kk21_stepwise_logistic_weights_cpp(const NumericMatrix& X,
												 const NumericVector& y,
												 const NumericVector& w) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p, NA_REAL);

	if (n == 0 || p == 0) {
	return weights;
	}

	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
	Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

	// Base design matrix X_S = [1, w, ...selected...], starting with [1, w]
	int m = 2;
	Eigen::MatrixXd XS(n, m);
	XS.col(0).setOnes();
	XS.col(1) = w_vec;

	Eigen::VectorXd beta = Eigen::VectorXd::Zero(m);
	Eigen::VectorXd p_hat, W_diag, resid;
	Eigen::LDLT<Eigen::MatrixXd> XtWX_ldlt;
	bool fit_ok = logistic_reduced_fit_for_score_test(XS, y_vec, beta, p_hat, W_diag, resid, XtWX_ldlt);

	std::vector<bool> used(p, false);

	for (int step = 0; step < p; ++step) {
	if (!fit_ok) break;

	// Precompute W * X_S once per step (used for each candidate's information)
	Eigen::MatrixXd WXS = (XS.array().colwise() * W_diag.array()).matrix(); // n x m

	double best_stat = -1.0;
	int best_j = -1;

	for (int j = 0; j < p; ++j) {
		if (used[j]) continue;

		Eigen::VectorXd xj = X_map.col(j);
		// Score: U_j = x_j' (y - p_hat)
		double Uj = xj.dot(resid);
		// Adjusted information: I_jj.S = x_j' W x_j - (X_S' W x_j)' (X_S' W X_S)^{-1} (X_S' W x_j)
		Eigen::VectorXd v_j = WXS.transpose() * xj;          // X_S' W x_j  (m-vector)
		Eigen::VectorXd h_j = XtWX_ldlt.solve(v_j);          // (X_S' W X_S)^{-1} v_j
		double xjWxj = (W_diag.array() * xj.array().square()).sum();
		double Ijj = xjWxj - v_j.dot(h_j);
		if (Ijj <= 0.0) continue;

		// |z_score| = |U_j| / sqrt(I_jj.S)  (asymptotically equivalent to |t| from full Wald)
		double stat = std::fabs(Uj) / std::sqrt(Ijj);
		if (!R_finite(stat)) stat = 0.0;

		if (stat > best_stat) {
		best_stat = stat;
		best_j = j;
		}
	}

	if (best_j < 0) break;

	weights[best_j] = best_stat;
	used[best_j] = true;

	// Update X_S: append x_{best_j}, refit with warm warm_start_params
	Eigen::MatrixXd XS_new(n, m + 1);
	XS_new.leftCols(m) = XS;
	XS_new.col(m) = X_map.col(best_j);
	XS = XS_new;
	m++;

	// Warm warm_start_params: extend beta with 0 for the newly added covariate
	Eigen::VectorXd beta_new(m);
	beta_new.head(m - 1) = beta;
	beta_new[m - 1] = 0.0;
	beta = beta_new;

	fit_ok = logistic_reduced_fit_for_score_test(XS, y_vec, beta, p_hat, W_diag, resid, XtWX_ldlt);
	}

	return weights;
}

// Helper: Univariate beta regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_beta_tstat(
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& x,
	int maxit = 100,
	double tol = 1e-8
) {
	int n = y.size();
	if (n < 3) return -1.0;

	// Build design matrix [1, x]
	Eigen::MatrixXd X(n, 2);
	X.col(0).setOnes();
	X.col(1) = x;

	// Fit logistic regression via IRLS (for mean model with logit link)
	Eigen::VectorXd beta = Eigen::VectorXd::Zero(2);
	Eigen::VectorXd prob(n);
	Eigen::VectorXd w(n);

	for (int iter = 0; iter < maxit; ++iter) {
		Eigen::VectorXd eta = X * beta;
		prob = (1.0 / (1.0 + (-eta.array()).exp())).min(1.0 - 1e-10).max(1e-10).matrix();
		w = (prob.array() * (1.0 - prob.array())).max(1e-10).matrix();

		Eigen::VectorXd z = eta + (y - prob).cwiseQuotient(w);
		Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
		Eigen::VectorXd XtWz = weighted_crossprod_rhs(X, w, z);

		Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
		if (ldlt.info() != Eigen::Success) return -1.0;

		Eigen::VectorXd beta_new = ldlt.solve(XtWz);
		if ((beta - beta_new).norm() < tol) {
			beta = beta_new;
			break;
		}
		beta = beta_new;
	}

	// Compute final mu and weights
	Eigen::VectorXd eta = X * beta;
	prob = (1.0 / (1.0 + (-eta.array()).exp())).min(1.0 - 1e-12).max(1e-12).matrix();
	w = (prob.array() * (1.0 - prob.array())).max(1e-10).matrix();

	// Profile likelihood optimization for phi
	// Grid search over log(phi) from log(0.01) to log(1000)
	double best_phi = 10.0;
	double best_ll = -std::numeric_limits<double>::infinity();

	for (double log_phi = -2.0; log_phi <= 7.0; log_phi += 0.5) {
		double phi = std::exp(log_phi);
		const double lgamma_phi = std::lgamma(phi);
		double ll = 0.0;
		for (int i = 0; i < n; ++i) {
			double mu_i = prob[i];
			double yi = y[i];
			double mu_phi = mu_i * phi;
			double one_minus_mu_phi = (1.0 - mu_i) * phi;
			mu_phi = std::max(1e-12, mu_phi);
			one_minus_mu_phi = std::max(1e-12, one_minus_mu_phi);

			ll += lgamma_phi - std::lgamma(mu_phi) - std::lgamma(one_minus_mu_phi)
				+ (mu_phi - 1.0) * std::log(yi)
				+ (one_minus_mu_phi - 1.0) * std::log1p(-yi);
		}
		if (ll > best_ll) {
			best_ll = ll;
			best_phi = phi;
		}
	}

	// Final weights for variance computation: w_final = (1+phi) * mu * (1-mu)
	Eigen::VectorXd w_final = ((1.0 + best_phi) * w.array()).matrix();

	// Compute XtWX and get variance of beta[1]
	Eigen::MatrixXd XtWX = weighted_crossprod(X, w_final);
	const double var_1 = compute_diagonal_inverse_entry(XtWX, 2);
	if (!R_finite(var_1) || var_1 <= 0) return -1.0;
	double se = std::sqrt(var_1);
	if (!std::isfinite(se) || se <= 0) return -1.0;

	return std::fabs(beta(1) / se);
}

// Helper: Univariate negative binomial regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_negbin_tstat(
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& x,
	int maxit = 50,
	double tol = 1e-6
) {
	int n = y.size();
	if (n < 3) return -1.0;

	// Build design matrix [1, x]
	Eigen::MatrixXd X(n, 2);
	X.col(0).setConstant(1.0);
	X.col(1) = x;

	// Initialize with Poisson regression coefficients (theta -> infinity)
	Eigen::VectorXd beta = Eigen::VectorXd::Zero(2);

	// Initialize theta using method of moments
	double y_mean = y.mean();
	double y_var = (y.array() - y_mean).square().sum() / (n - 1);
	double theta = (y_var > y_mean) ? (y_mean * y_mean) / (y_var - y_mean) : 10.0;
	theta = std::max(0.01, std::min(1000.0, theta));

	Eigen::VectorXd mu(n);
	Eigen::VectorXd w(n);

	for (int outer = 0; outer < maxit; ++outer) {
		// E-step / IRLS for beta given theta
		for (int iter = 0; iter < 25; ++iter) {
			Eigen::VectorXd eta = X * beta;
			mu = eta.array().min(20.0).exp().matrix();
			w = (mu.array() / (1.0 + mu.array() / theta)).max(1e-10).matrix();

			Eigen::VectorXd z = eta + (y - mu).cwiseQuotient(mu);
			Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
			Eigen::VectorXd XtWz = weighted_crossprod_rhs(X, w, z);

			Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
			if (ldlt.info() != Eigen::Success) return -1.0;

			Eigen::VectorXd beta_new = ldlt.solve(XtWz);
			if ((beta - beta_new).norm() < tol) {
				beta = beta_new;
				break;
			}
			beta = beta_new;
		}

		// Update mu with final beta
		Eigen::VectorXd eta = X * beta;
		mu = eta.array().min(20.0).exp().matrix();

		// M-step: update theta using one-step Newton-Raphson
		double score = 0.0;
		double info = 0.0;
		const double digamma_theta = fast_digamma(theta);
		const double trigamma_theta = R::trigamma(theta);
		const double log_theta = std::log(theta);
		const double inv_theta = 1.0 / theta;
		for (int i = 0; i < n; ++i) {
			double yi = y[i];
			double mui = mu[i];
			score += fast_digamma(yi + theta) - digamma_theta
					 + log_theta - std::log(theta + mui) + 1.0
					 - (yi + theta) / (theta + mui);
			info += -R::trigamma(yi + theta) + trigamma_theta
					- inv_theta + 2.0 / (theta + mui)
					- (yi + theta) / ((theta + mui) * (theta + mui));
		}

		if (std::fabs(info) < 1e-10) break;
		double theta_new = theta - score / info;
		theta_new = std::max(0.01, std::min(1000.0, theta_new));

		if (std::fabs(theta_new - theta) < tol) {
			theta = theta_new;
			break;
		}
		theta = theta_new;
	}

	// Compute final weights and information matrix
	Eigen::VectorXd eta = X * beta;
	mu = eta.array().min(20.0).exp().matrix();
	w = (mu.array() / (1.0 + mu.array() / theta)).max(1e-10).matrix();

	Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
	const double var_1 = compute_diagonal_inverse_entry(XtWX, 2);
	if (!R_finite(var_1) || var_1 <= 0) return -1.0;
	double se = std::sqrt(var_1);
	if (!std::isfinite(se) || se <= 0) return -1.0;

	return std::fabs(beta(1) / se);
}

// [[Rcpp::export]]
NumericVector kk21_beta_weights_cpp(const NumericMatrix& X,
									const NumericVector& y) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p);
	double eps = std::numeric_limits<double>::epsilon();

	if (n < 3 || p == 0) {
		std::fill(weights.begin(), weights.end(), eps);
		return weights;
	}

	// Convert to Eigen
	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

	for (int j = 0; j < p; ++j) {
		double tstat = univariate_beta_tstat(y_vec, X_map.col(j));

		if (tstat > 0 && std::isfinite(tstat)) {
			weights[j] = tstat;
		} else {
			// Fall back to continuous weights on logit transform
			Eigen::VectorXd y_logit(n);
			for (int i = 0; i < n; ++i) {
				double yi = std::max(1e-10, std::min(1.0 - 1e-10, y_vec[i]));
				y_logit[i] = std::log(yi / (1.0 - yi));
			}

			// Simple OLS
			double sumy = y_logit.sum();
			double sumy2 = y_logit.squaredNorm();
			double ybar = sumy / n;
			double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
			for (int i = 0; i < n; ++i) {
				double xi = X_map(i, j);
				sumx += xi;
				sumx2 += xi * xi;
				sumxy += xi * y_logit[i];
			}
			double xbar = sumx / n;
			double varx = sumx2 - n * xbar * xbar;

			if (varx <= eps || n <= 2) {
				weights[j] = eps;
				continue;
			}

			double b1 = (sumxy - n * xbar * ybar) / varx;
			double b0 = ybar - b1 * xbar;
			double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2
						 + 2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

			if (n <= 2 || sse <= 0) {
				weights[j] = eps;
				continue;
			}

			double sigma2 = sse / (n - 2);
			double se_b1 = std::sqrt(sigma2 / varx);

			if (!std::isfinite(se_b1) || se_b1 <= 0) {
				weights[j] = eps;
			} else {
				weights[j] = std::fabs(b1 / se_b1);
			}
		}
	}

	return weights;
}

// [[Rcpp::export]]
NumericVector kk21_negbin_weights_cpp(const NumericMatrix& X,
										const NumericVector& y) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p);
	double eps = std::numeric_limits<double>::epsilon();

	if (n < 3 || p == 0) {
		std::fill(weights.begin(), weights.end(), eps);
		return weights;
	}

	// Convert to Eigen
	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

	// Pre-compute OLS quantities for fallback (using log(y+1))
	Eigen::VectorXd log_y1(n);
	for (int i = 0; i < n; ++i) {
		log_y1[i] = std::log1p(y_vec[i]);
	}
	double sumy = log_y1.sum();
	double sumy2 = log_y1.squaredNorm();
	double ybar = sumy / n;

	for (int j = 0; j < p; ++j) {
		double tstat = univariate_negbin_tstat(y_vec, X_map.col(j));

		if (tstat > 0 && std::isfinite(tstat)) {
			weights[j] = tstat;
			continue;
		}

		// Fall back to OLS on log(y+1)
		double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
		for (int i = 0; i < n; ++i) {
			double x = X_map(i, j);
			sumx += x;
			sumx2 += x * x;
			sumxy += x * log_y1[i];
		}
		double xbar = sumx / n;
		double varx = sumx2 - n * xbar * xbar;

		if (varx <= eps || n <= 2) {
			weights[j] = eps;
			continue;
		}

		double b1 = (sumxy - n * xbar * ybar) / varx;
		double b0 = ybar - b1 * xbar;
		double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2
					 + 2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

		if (n <= 2 || sse <= 0) {
			weights[j] = eps;
			continue;
		}

		double sigma2 = sse / (n - 2);
		double se_b1 = std::sqrt(sigma2 / varx);

		if (!std::isfinite(se_b1) || se_b1 <= 0) {
			weights[j] = eps;
		} else {
			weights[j] = std::fabs(b1 / se_b1);
		}
	}

	return weights;
}

// Helper: Univariate Weibull AFT regression for a single covariate
// Returns t-statistic for the covariate, or -1 if fitting fails
static double univariate_weibull_tstat(
	const Eigen::VectorXd& log_y,
	const Eigen::VectorXd& delta,
	const Eigen::VectorXd& x,
	int maxit = 30,
	double tol = 1e-5
) {
	int n = log_y.size();
	if (n < 4) return -1.0;

	double n_events = delta.sum();
	if (n_events < 2) return -1.0;

	// Build design matrix [1, x]
	Eigen::MatrixXd X(n, 2);
	X.col(0).setConstant(1.0);
	X.col(1) = x;

	// Initialize with OLS on log(y)
	Eigen::MatrixXd XtX = X.transpose() * X;
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(XtX);
	if (qr.rank() < 2) return -1.0;

	Eigen::VectorXd beta = qr.solve(X.transpose() * log_y);
	Eigen::VectorXd resid = log_y - X * beta;

	// Initial scale estimate
	double scale = 0.0;
	for (int i = 0; i < n; ++i) {
		if (delta(i) > 0.5) {
			scale += resid(i) * resid(i);
		}
	}
	scale = std::sqrt(scale / std::max(1.0, n_events - 2));
	scale = std::max(0.01, std::min(10.0, scale));

	// Newton-Raphson for Weibull AFT
	for (int iter = 0; iter < maxit; ++iter) {
		Eigen::ArrayXd z = clamp_weibull_z_kk((resid / scale).array());
		Eigen::VectorXd w = z.exp().matrix();
		Eigen::ArrayXd adj = (delta.array() > 0.5).select(z - 1.0 + w.array(), w.array());

		// Weighted least squares update
		Eigen::VectorXd pseudo_response = (log_y.array() - scale * adj / w.array()).matrix();
		double sum_w = 0.0;
		double sum_wx = 0.0;
		double sum_wxx = 0.0;
		double sum_wz = 0.0;
		double sum_wxz = 0.0;
		for (int i = 0; i < n; ++i) {
			const double wi = w[i];
			const double xi = x[i];
			const double zi = pseudo_response[i];
			sum_w += wi;
			sum_wx += wi * xi;
			sum_wxx += wi * xi * xi;
			sum_wz += wi * zi;
			sum_wxz += wi * xi * zi;
		}
		Eigen::Matrix2d XwXw;
		XwXw << sum_w, sum_wx,
		          sum_wx, sum_wxx;
		Eigen::Vector2d Xwz(sum_wz, sum_wxz);
		Eigen::ColPivHouseholderQR<Eigen::Matrix2d> qr_w(XwXw);
		if (qr_w.rank() < 2) return -1.0;

		Eigen::VectorXd beta_new = qr_w.solve(Xwz);
		Eigen::VectorXd resid_new = log_y - X * beta_new;

		// Scale update
		double sum_d = 0.0;
		double scale_new = scale;
		for (int i = 0; i < n; ++i) {
			if (delta(i) > 0.5) sum_d += 1.0;
		}
		if (sum_d > 0) {
			const Eigen::ArrayXd z_new = clamp_weibull_z_kk((resid_new / scale).array());
			const Eigen::ArrayXd exp_z_new = z_new.exp();
			double score =
				(delta.array() * (z_new - 1.0) - z_new * exp_z_new).sum() / scale;
			double info = sum_d / (scale * scale);
			if (std::fabs(info) > 1e-10) {
				scale_new = scale - score / info;
				scale_new = std::max(0.01, std::min(10.0, scale_new));
			}
		}

		double diff = (beta_new - beta).norm() + std::fabs(scale_new - scale);
		beta = beta_new;
		scale = scale_new;
		resid = resid_new;

		if (diff < tol) break;
	}

	// Compute standard error of beta[1]
	Eigen::VectorXd w = (clamp_weibull_z_kk((resid / scale).array()).exp() / (scale * scale)).matrix();
	Eigen::MatrixXd info_mat = weighted_crossprod(X, w);

	const double var_1 = compute_diagonal_inverse_entry(info_mat, 2);
	if (!R_finite(var_1) || var_1 <= 0) return -1.0;
	double se = std::sqrt(var_1);
	if (!std::isfinite(se) || se <= 0) return -1.0;

	return std::fabs(beta(1) / se);
}

// [[Rcpp::export]]
NumericVector kk21_survival_weights_cpp(const NumericMatrix& X,
										const NumericVector& y,
										const NumericVector& delta) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p);
	double eps = std::numeric_limits<double>::epsilon();

	if (n < 2 || p == 0) {
		std::fill(weights.begin(), weights.end(), eps);
		return weights;
	}

	// Convert to Eigen
	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
	Eigen::VectorXd delta_vec = as<Eigen::VectorXd>(delta);

	// Compute log(y)
	Eigen::VectorXd log_y(n);
	for (int i = 0; i < n; ++i) {
		log_y(i) = std::log(std::max(y_vec(i), 1e-10));
	}

	// Pre-compute OLS quantities for fallback
	double sumy = log_y.sum();
	double sumy2 = log_y.squaredNorm();
	double ybar = sumy / n;

	for (int j = 0; j < p; ++j) {
		// Try Weibull first
		double tstat = univariate_weibull_tstat(log_y, delta_vec, X_map.col(j));

		if (tstat > 0 && std::isfinite(tstat)) {
			weights[j] = tstat;
			continue;
		}

		// Fall back to OLS on log(y)
		double sumx = 0.0, sumx2 = 0.0, sumxy = 0.0;
		for (int i = 0; i < n; ++i) {
			double x = X_map(i, j);
			sumx += x;
			sumx2 += x * x;
			sumxy += x * log_y(i);
		}
		double xbar = sumx / n;
		double varx = sumx2 - n * xbar * xbar;

		if (varx <= eps || n <= 2) {
			weights[j] = eps;
			continue;
		}

		double b1 = (sumxy - n * xbar * ybar) / varx;
		double b0 = ybar - b1 * xbar;
		double sse = sumy2 + n * b0 * b0 + b1 * b1 * sumx2 +
					 2.0 * b0 * b1 * sumx - 2.0 * b0 * sumy - 2.0 * b1 * sumxy;

		if (n <= 2 || sse <= 0) {
			weights[j] = eps;
			continue;
		}

		double sigma2 = sse / (n - 2);
		double se_b1 = std::sqrt(sigma2 / varx);

		if (!std::isfinite(se_b1) || se_b1 <= 0) {
			weights[j] = eps;
		} else {
			weights[j] = std::fabs(b1 / se_b1);
		}
	}

	return weights;
}

// Helper: Multivariate beta regression returning t-statistic for coefficient at index coef_idx
// Design matrix X should include intercept as first column
static double multivariate_beta_tstat(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	int coef_idx = 1,
	int maxit = 100,
	double tol = 1e-8
) {
	int n = X.rows();
	int p = X.cols();
	if (n < p + 1 || p < 2) return -1.0;

	// Fit logistic regression via IRLS (for mean model with logit link)
	Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
	Eigen::VectorXd prob(n);
	Eigen::VectorXd w(n);

	for (int iter = 0; iter < maxit; ++iter) {
		Eigen::VectorXd eta = X * beta;
		prob = (1.0 / (1.0 + (-eta.array()).exp())).min(1.0 - 1e-10).max(1e-10).matrix();
		w = (prob.array() * (1.0 - prob.array())).max(1e-10).matrix();

		Eigen::VectorXd z = eta + (y - prob).cwiseQuotient(w);
		Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
		Eigen::VectorXd XtWz = weighted_crossprod_rhs(X, w, z);

		Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
		if (ldlt.info() != Eigen::Success) return -1.0;

		Eigen::VectorXd beta_new = ldlt.solve(XtWz);
		if ((beta - beta_new).norm() < tol) {
			beta = beta_new;
			break;
		}
		beta = beta_new;
	}

	// Compute final mu and weights
	Eigen::VectorXd eta = X * beta;
	prob = (1.0 / (1.0 + (-eta.array()).exp())).min(1.0 - 1e-12).max(1e-12).matrix();
	w = (prob.array() * (1.0 - prob.array())).max(1e-10).matrix();

	// Profile likelihood optimization for phi (grid search)
	double best_phi = 10.0;
	double best_ll = -std::numeric_limits<double>::infinity();

	for (double log_phi = -2.0; log_phi <= 7.0; log_phi += 0.5) {
		double phi = std::exp(log_phi);
		const double lgamma_phi = std::lgamma(phi);
		double ll = 0.0;
		for (int i = 0; i < n; ++i) {
			double mu_i = prob[i];
			double yi = y[i];
			double mu_phi = mu_i * phi;
			double one_minus_mu_phi = (1.0 - mu_i) * phi;
			mu_phi = std::max(1e-12, mu_phi);
			one_minus_mu_phi = std::max(1e-12, one_minus_mu_phi);

			ll += lgamma_phi - std::lgamma(mu_phi) - std::lgamma(one_minus_mu_phi)
				+ (mu_phi - 1.0) * std::log(std::max(yi, 1e-12))
				+ (one_minus_mu_phi - 1.0) * std::log1p(-std::min(yi, 1.0 - 1e-12));
		}
		if (ll > best_ll) {
			best_ll = ll;
			best_phi = phi;
		}
	}

	// Final weights for variance computation
	Eigen::VectorXd w_final = ((1.0 + best_phi) * w.array()).matrix();

	// Compute XtWX and get variance of beta[coef_idx]
	Eigen::MatrixXd XtWX = weighted_crossprod(X, w_final);
	if (coef_idx >= p) return -1.0;
	const double var_k = compute_diagonal_inverse_entry(XtWX, coef_idx + 1);
	if (!R_finite(var_k) || var_k <= 0) return -1.0;
	double se = std::sqrt(var_k);
	if (!std::isfinite(se) || se <= 0) return -1.0;

	return std::fabs(beta(coef_idx) / se);
}

// Helper: Multivariate negative binomial regression returning t-statistic for coefficient at index coef_idx
// Design matrix X should include intercept as first column
static double multivariate_negbin_tstat(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	int coef_idx = 1,
	int maxit = 50,
	double tol = 1e-6
) {
	int n = X.rows();
	int p = X.cols();
	if (n < p + 1 || p < 2) return -1.0;

	Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);

	// Initialize theta using method of moments
	double y_mean = y.mean();
	double y_var = (y.array() - y_mean).square().sum() / (n - 1);
	double theta = (y_var > y_mean) ? (y_mean * y_mean) / (y_var - y_mean) : 10.0;
	theta = std::max(0.01, std::min(1000.0, theta));

	Eigen::VectorXd mu(n);
	Eigen::VectorXd w(n);

	for (int outer = 0; outer < maxit; ++outer) {
		// IRLS for beta given theta
		for (int iter = 0; iter < 25; ++iter) {
			Eigen::VectorXd eta = X * beta;
			mu = eta.array().min(20.0).exp().matrix();
			w = (mu.array() / (1.0 + mu.array() / theta)).max(1e-10).matrix();

			Eigen::VectorXd z = eta + (y - mu).cwiseQuotient(mu);
			Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
			Eigen::VectorXd XtWz = weighted_crossprod_rhs(X, w, z);

			Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
			if (ldlt.info() != Eigen::Success) return -1.0;

			Eigen::VectorXd beta_new = ldlt.solve(XtWz);
			if ((beta - beta_new).norm() < tol) {
				beta = beta_new;
				break;
			}
			beta = beta_new;
		}

		// Update mu
		Eigen::VectorXd eta = X * beta;
		mu = eta.array().min(20.0).exp().matrix();

		// Update theta using Newton-Raphson
		double score = 0.0;
		double info = 0.0;
		const double digamma_theta = fast_digamma(theta);
		const double trigamma_theta = R::trigamma(theta);
		const double log_theta = std::log(theta);
		const double inv_theta = 1.0 / theta;
		for (int i = 0; i < n; ++i) {
			double yi = y[i];
			double mui = mu[i];
			score += fast_digamma(yi + theta) - digamma_theta
					 + log_theta - std::log(theta + mui) + 1.0
					 - (yi + theta) / (theta + mui);
			info += -R::trigamma(yi + theta) + trigamma_theta
					- inv_theta + 2.0 / (theta + mui)
					- (yi + theta) / ((theta + mui) * (theta + mui));
		}

		if (std::fabs(info) < 1e-10) break;
		double theta_new = theta - score / info;
		theta_new = std::max(0.01, std::min(1000.0, theta_new));

		if (std::fabs(theta_new - theta) < tol) {
			theta = theta_new;
			break;
		}
		theta = theta_new;
	}

	// Compute final weights
	Eigen::VectorXd eta = X * beta;
	mu = eta.array().min(20.0).exp().matrix();
	w = (mu.array() / (1.0 + mu.array() / theta)).max(1e-10).matrix();

	Eigen::MatrixXd XtWX = weighted_crossprod(X, w);
	if (coef_idx >= p) return -1.0;
	const double var_k = compute_diagonal_inverse_entry(XtWX, coef_idx + 1);
	if (!R_finite(var_k) || var_k <= 0) return -1.0;
	double se = std::sqrt(var_k);
	if (!std::isfinite(se) || se <= 0) return -1.0;

	return std::fabs(beta(coef_idx) / se);
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_beta_weights_cpp(const NumericMatrix& X,
											 const NumericVector& y,
											 const NumericVector& w) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p, NA_REAL);

	if (n == 0 || p == 0) {
		return weights;
	}

	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
	Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

	// Clamp y to (0,1)
	for (int i = 0; i < n; ++i) {
		y_vec[i] = std::max(1e-10, std::min(1.0 - 1e-10, y_vec[i]));
	}

	std::vector<int> selected;
	selected.reserve(p);
	std::vector<bool> used(p, false);

	for (int step = 0; step < p; ++step) {
		double best_stat = -1.0;
		int best_j = -1;
		int k = static_cast<int>(selected.size());

		for (int j = 0; j < p; ++j) {
			if (used[j]) {
				continue;
			}

			// Build design matrix: [1, x_j, selected_covs..., w]
			Eigen::MatrixXd X(n, k + 3);
			X.col(0).setConstant(1.0);
			X.col(1) = X_map.col(j);
			for (int idx = 0; idx < k; ++idx) {
			        X.col(2 + idx) = X_map.col(selected[idx]);
			}
			X.col(k + 2) = w_vec;
			double stat = multivariate_beta_tstat(X, y_vec, 1);
			if (!R_finite(stat) || stat < 0) {
				stat = 0.0;
			}
			if (stat > best_stat) {
				best_stat = stat;
				best_j = j;
			}
		}

		if (best_j < 0) {
			break;
		}

		weights[best_j] = best_stat;
		used[best_j] = true;
		selected.push_back(best_j);
	}

	return weights;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_negbin_weights_cpp(const NumericMatrix& X,
												 const NumericVector& y,
												 const NumericVector& w) {
	int n = X.nrow();
	int p = X.ncol();
	NumericVector weights(p, NA_REAL);

	if (n == 0 || p == 0) {
		return weights;
	}

	Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
	Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
	Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

	std::vector<int> selected;
	selected.reserve(p);
	std::vector<bool> used(p, false);

	for (int step = 0; step < p; ++step) {
		double best_stat = -1.0;
		int best_j = -1;
		int k = static_cast<int>(selected.size());

		for (int j = 0; j < p; ++j) {
			if (used[j]) {
				continue;
			}

			// Build design matrix: [1, x_j, selected_covs..., w]
			Eigen::MatrixXd X(n, k + 3);
			X.col(0).setConstant(1.0);
			X.col(1) = X_map.col(j);
			for (int idx = 0; idx < k; ++idx) {
			        X.col(2 + idx) = X_map.col(selected[idx]);
			}
			X.col(k + 2) = w_vec;
			double stat = multivariate_negbin_tstat(X, y_vec, 1);
			if (!R_finite(stat) || stat < 0) {
				stat = 0.0;
			}
			if (stat > best_stat) {
				best_stat = stat;
				best_j = j;
			}
		}

		if (best_j < 0) {
			break;
		}

		weights[best_j] = best_stat;
		used[best_j] = true;
		selected.push_back(best_j);
	}

	return weights;
}

// Helper for ordinal t-stat
static double multivariate_ordinal_tstat(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    int coef_idx = 1
) {
    Function f("fast_ordinal_regression_with_var_cpp");
    List res = f(wrap(X), wrap(y));
    NumericVector b = res["b"];
    double ssq = res["ssq_b_2"];
    return std::fabs(b[coef_idx - 1] / std::sqrt(ssq));
}

// [[Rcpp::export]]
NumericVector kk21_ordinal_weights_cpp(const NumericMatrix& X,
                                         const NumericVector& y) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p);
    double eps = std::numeric_limits<double>::epsilon();

    if (n == 0 || p == 0) {
        std::fill(weights.begin(), weights.end(), eps);
        return weights;
    }

    Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);

    for (int j = 0; j < p; ++j) {
        Eigen::MatrixXd X(n, 1);
        X.col(0) = X_map.col(j);
        
        try {
            weights[j] = multivariate_ordinal_tstat(X, y_vec, 1);
        } catch (...) {
            weights[j] = eps;
        }
    }

    return weights;
}

// [[Rcpp::export]]
NumericVector kk21_stepwise_ordinal_weights_cpp(const NumericMatrix& X,
                                                 const NumericVector& y,
                                                 const NumericVector& w) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector weights(p, NA_REAL);

    if (n == 0 || p == 0) {
        return weights;
    }

    Eigen::MatrixXd X_map = as<Eigen::MatrixXd>(X);
    Eigen::VectorXd y_vec = as<Eigen::VectorXd>(y);
    Eigen::VectorXd w_vec = as<Eigen::VectorXd>(w);

    std::vector<int> selected;
    selected.reserve(p);
    std::vector<bool> used(p, false);

    for (int step = 0; step < p; ++step) {
        double best_stat = -1.0;
        int best_j = -1;
        int k = static_cast<int>(selected.size());

        for (int j = 0; j < p; ++j) {
            if (used[j]) {
                continue;
            }

            // Build design matrix: [x_j, selected_covs..., w]
            Eigen::MatrixXd X(n, k + 2);
            X.col(0) = X_map.col(j);
            for (int idx = 0; idx < k; ++idx) {
                X.col(1 + idx) = X_map.col(selected[idx]);
            }
            X.col(k + 1) = w_vec;

            try {
                double stat = multivariate_ordinal_tstat(X, y_vec, 1);
                if (!R_finite(stat) || stat < 0) {
                    stat = 0.0;
                }
                if (stat > best_stat) {
                    best_stat = stat;
                    best_j = j;
                }
            } catch (...) {
                continue;
            }
        }

        if (best_j < 0) {
            break;
        }

        weights[best_j] = best_stat;
        used[best_j] = true;
        selected.push_back(best_j);
    }

	return weights;
}
