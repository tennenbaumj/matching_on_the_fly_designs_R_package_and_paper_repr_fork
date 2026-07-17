#include "_helper_functions.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

inline Eigen::ArrayXd plogis_array_clamped(const Eigen::ArrayXd& eta) {
    const Eigen::ArrayXd eta_clamped = eta.max(-20.0).min(20.0);
    return 1.0 / (1.0 + (-eta_clamped).exp());
}

struct ContinuationRatioObjective {
    const Eigen::Ref<const MatrixXd> X_aug;
    const Eigen::Ref<const VectorXd> z;
    VectorXd eta;
    VectorXd mu;
    VectorXd work;
    ArrayXd log_mu;
    ArrayXd log_one_minus_mu;

    ContinuationRatioObjective(const Eigen::Ref<const MatrixXd>& X_aug, const Eigen::Ref<const VectorXd>& z) :
        X_aug(X_aug), z(z), eta(X_aug.rows()), mu(X_aug.rows()), work(X_aug.rows()),
        log_mu(X_aug.rows()), log_one_minus_mu(X_aug.rows()) {}

    double operator()(const VectorXd& beta, VectorXd& grad) {
        eta.noalias() = X_aug * beta;
        mu = plogis_array_clamped(eta.array()).matrix();
        work.noalias() = mu - z;
        grad.noalias() = X_aug.transpose() * work; // Negative log-likelihood gradient

        log_mu = mu.array().max(1e-12).log();
        log_one_minus_mu = (1.0 - mu.array()).max(1e-12).log();
        return -(z.array() * log_mu + (1.0 - z.array()) * log_one_minus_mu).sum();
    }

    MatrixXd hessian(const VectorXd& beta) {
        eta.noalias() = X_aug * beta;
        mu = plogis_array_clamped(eta.array()).matrix();
        work = (mu.array() * (1.0 - mu.array())).matrix();
        return weighted_crossprod(X_aug, work);
    }
};

static List build_continuation_ratio_augmented_data(const Eigen::Ref<const MatrixXd>& X,
													const Eigen::Ref<const VectorXd>& y) {
	int n = X.rows();
	int p = X.cols();

	std::vector<double> levels;
	for (int i = 0; i < y.size(); ++i) {
		if (std::find(levels.begin(), levels.end(), y[i]) == levels.end()) {
			levels.push_back(y[i]);
		}
	}
	std::sort(levels.begin(), levels.end());
	int K = levels.size();
	if (K < 2) {
		return List::create(Named("X_aug") = MatrixXd(0, p), Named("z") = VectorXd(0), Named("n_alpha") = 0);
	}
	int n_alpha = K - 1;

	std::vector<int> y_level(n);
	int total_rows = 0;
	for (int i = 0; i < n; ++i) {
		y_level[i] = static_cast<int>(
			std::lower_bound(levels.begin(), levels.end(), y[i]) - levels.begin());
		total_rows += std::min(y_level[i] + 1, n_alpha);
	}

	MatrixXd X_aug = MatrixXd::Zero(total_rows, n_alpha + p);
	VectorXd z(total_rows);
	int row = 0;
	for (int i = 0; i < n; ++i) {
		const int yi_level = y_level[i];
		const int rows_i = std::min(yi_level + 1, n_alpha);
		for (int j = 0; j < rows_i; ++j, ++row) {
			X_aug(row, j) = 1.0;
			if (p > 0) X_aug.row(row).tail(p) = X.row(i);
			z[row] = (yi_level == j) ? 1.0 : 0.0;
		}
	}
	return List::create(Named("X_aug") = X_aug, Named("z") = z, Named("n_alpha") = n_alpha);
}

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_continuation_ratio_regression_score_cpp(SEXP X_sexp,
															SEXP y_sexp,
															SEXP params_sexp) {
	NumericMatrix X_r(X_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	List aug = build_continuation_ratio_augmented_data(X, y);
	MatrixXd X_aug = aug["X_aug"];
	VectorXd z = aug["z"];
	if (X_aug.rows() == 0) return VectorXd::Zero(params.size());
	VectorXd eta = X_aug * params;
	VectorXd mu = plogis_array_clamped(eta.array()).matrix();
	return X_aug.transpose() * (z - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_continuation_ratio_regression_hessian_cpp(SEXP X_sexp,
															  SEXP y_sexp,
															  SEXP params_sexp) {
	NumericMatrix X_r(X_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	List aug = build_continuation_ratio_augmented_data(X, y);
	MatrixXd X_aug = aug["X_aug"];
	if (X_aug.rows() == 0) return MatrixXd::Zero(params.size(), params.size());
	VectorXd eta = X_aug * params;
	VectorXd mu = plogis_array_clamped(eta.array()).matrix();
	VectorXd w = mu.array() * (1.0 - mu.array());
	return -weighted_crossprod(X_aug, w);
}

//' @title Fast Continuation-Ratio Regression (C++)
//' @description High-performance continuation-ratio logit model fitting.
//' @param X A numeric matrix of predictors (no intercept column; threshold intercepts are estimated internally).
//' @param y A numeric vector of ordinal responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when no warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @return A list containing coefficients, alpha, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_continuation_ratio_regression_cpp(SEXP X_sexp, SEXP y_sexp, int maxit = 100, double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                             bool smart_cold_start = true,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "lbfgs",
                                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    int p = X.cols();
    List aug = build_continuation_ratio_augmented_data(X, y);
    MatrixXd X_aug = aug["X_aug"];
    VectorXd z = aug["z"];
    int n_alpha = aug["n_alpha"];
    if (n_alpha == 0) {
        return List::create(Named("b") = VectorXd::Zero(p), Named("alpha") = VectorXd::Zero(0));
    }
    
    int p_aug = n_alpha + p;
    ContinuationRatioObjective fun(X_aug, z);
    VectorXd beta = VectorXd::Zero(p_aug);
    if (warm_start_beta.isNotNull()) {
        beta = as<VectorXd>(warm_start_beta);
        if (beta.size() != p_aug) stop("warm_start_beta size mismatch");
    } else if (smart_cold_start) {
        // Smart warm_start_params: OLS on z (the augmented binary response)
        beta = ols_smart_cold_start_beta(X_aug, z);
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p_aug, fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_start_ptr);
    
    return List::create(
        Named("b") = fit.params.tail(p),
        Named("alpha") = fit.params.head(n_alpha),
        Named("beta_full") = fit.params,
        Named("params") = fit.params,
        Named("neg_loglik") = fit.value,
        Named("X_aug") = X_aug,
        Named("z") = z,
        Named("converged") = fit.converged,
        Named("fisher_information") = fun.hessian(fit.params)
    );
}

// [[Rcpp::export]]
List fast_continuation_ratio_regression_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int maxit = 100, double tol = 1e-8,
                                                      Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                      bool smart_cold_start = true,
                                                      Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                      Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                      std::string optimization_alg = "lbfgs",
                                                      Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    int p = X.cols();
    List aug = build_continuation_ratio_augmented_data(X, y);
    MatrixXd X_aug = aug["X_aug"];
    VectorXd z = aug["z"];
    int n_alpha = aug["n_alpha"];
    if (n_alpha == 0) {
         return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_j") = NA_REAL, Named("converged") = false);
    }
    
    int p_aug = n_alpha + p;
    ContinuationRatioObjective fun(X_aug, z);
    VectorXd beta = VectorXd::Zero(p_aug);
    if (warm_start_beta.isNotNull()) {
        beta = as<VectorXd>(warm_start_beta);
        if (beta.size() != p_aug) stop("warm_start_beta size mismatch");
    } else if (smart_cold_start) {
        beta = ols_smart_cold_start_beta(X_aug, z);
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p_aug, fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_start_ptr);
    
    MatrixXd info = fun.hessian(fit.params);
    MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
    int free_j = -1;
    for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
        if (fixed_spec.free_idx[jj] == n_alpha) { free_j = jj + 1; break; }
    double ssq_b_j = (p >= 1 && free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;

    SEXP vcov_sexp = R_NilValue;
    if (fit.converged) {
        MatrixXd cov_free = covariance_from_information(info_free);
        MatrixXd vcov = expand_free_covariance(n_alpha + p, fixed_spec, cov_free, true);
        vcov_sexp = Rcpp::wrap(vcov);
    }
    return List::create(
        Named("b") = fit.params.tail(p),
        Named("ssq_b_j") = ssq_b_j,
        Named("neg_loglik") = fit.value,
        Named("vcov") = vcov_sexp,
        Named("converged") = fit.converged,
        Named("params") = fit.params,
        Named("fisher_information") = info
    );
}
