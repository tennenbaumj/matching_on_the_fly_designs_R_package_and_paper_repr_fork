#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

struct DigammaFunctor {
	double operator()(double x) const {
		return fast_digamma(x);
	}
};

struct TrigammaFunctor {
	double operator()(double x) const {
		return R::trigamma(x);
	}
};

class BetaRegression {
private:
	const Eigen::Ref<const Eigen::VectorXd> m_y;
	const Eigen::Ref<const Eigen::MatrixXd> m_X;
	const Eigen::VectorXd m_weights;
	const int m_n;
	const int m_p;
	const double m_weight_sum;
	const Eigen::VectorXd m_log_y;
	const Eigen::VectorXd m_log1_y;
	Eigen::VectorXd m_eta;
	Eigen::VectorXd m_mu;
	Eigen::VectorXd m_w_grad;

public:
	BetaRegression(const Eigen::Ref<const Eigen::VectorXd>& y, const Eigen::Ref<const Eigen::MatrixXd>& X) :
		m_y(y), m_X(X), m_weights(Eigen::VectorXd::Ones(X.rows())), m_n(X.rows()), m_p(X.cols()),
		m_weight_sum(static_cast<double>(X.rows())),
		m_log_y(y.array().log().matrix()),
		m_log1_y((1.0 - y.array()).log().matrix()),
		m_eta(X.rows()),
		m_mu(X.rows()),
		m_w_grad(X.rows()) {}

	BetaRegression(const Eigen::Ref<const Eigen::VectorXd>& y, const Eigen::Ref<const Eigen::MatrixXd>& X,
	               const Eigen::Ref<const Eigen::VectorXd>& weights) :
		m_y(y), m_X(X), m_weights(weights), m_n(X.rows()), m_p(X.cols()),
		m_weight_sum(weights.sum()),
		m_log_y(y.array().log().matrix()),
		m_log1_y((1.0 - y.array()).log().matrix()),
		m_eta(X.rows()),
		m_mu(X.rows()),
		m_w_grad(X.rows()) {}

	double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
		const auto beta = params.head(m_p);
		const double phi = std::exp(params[m_p]);
		const double lgamma_phi = fast_lgamma(phi);
		const double digamma_phi = fast_digamma(phi);
		const double epsilon = 1e-8;

		m_eta.noalias() = m_X * beta;

		double neg_ll = 0.0;
		double d_neg_ll_d_phi = -m_weight_sum * digamma_phi;
		for (int i = 0; i < m_n; ++i) {
			double mui = 1.0 / (1.0 + std::exp(-m_eta[i]));
			if (mui < epsilon) {
				mui = epsilon;
			} else if (mui > 1.0 - epsilon) {
				mui = 1.0 - epsilon;
			}
			m_mu[i] = mui;

			const double one_minus_mu = 1.0 - mui;
			const double a = mui * phi;
			const double b = one_minus_mu * phi;
			const double log_y = m_log_y[i];
			const double log1_y = m_log1_y[i];
			const double weight = m_weights[i];

			neg_ll += weight * (
				-lgamma_phi +
				fast_lgamma(a) +
				fast_lgamma(b) -
				(a - 1.0) * log_y -
				(b - 1.0) * log1_y
			);

			const double dig_a = fast_digamma(a);
			const double dig_b = fast_digamma(b);
			const double C = dig_a - dig_b - log_y + log1_y;
			const double d_mu_d_eta = mui * one_minus_mu;
			m_w_grad[i] = weight * phi * C * d_mu_d_eta;

			d_neg_ll_d_phi += weight * (
				mui * dig_a +
				one_minus_mu * dig_b -
				mui * log_y -
				one_minus_mu * log1_y
			);
		}

		grad.resize(m_p + 1);
		grad.head(m_p).noalias() = m_X.transpose() * m_w_grad;
		grad[m_p] = d_neg_ll_d_phi * phi;

		return neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
		int total_p = m_p + 1;
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
		Eigen::VectorXd beta = params.head(m_p);
		double phi = std::exp(params[m_p]);
		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
		double epsilon = 1e-8;
		for (int i = 0; i < m_n; ++i) {
			if (mu[i] < epsilon) mu[i] = epsilon;
			if (mu[i] > 1.0 - epsilon) mu[i] = 1.0 - epsilon;
		}

		Eigen::VectorXd a = mu.array() * phi;
		Eigen::VectorXd b = (1.0 - mu.array()) * phi;
		Eigen::VectorXd dig_a = a.unaryExpr(DigammaFunctor());
		Eigen::VectorXd dig_b = b.unaryExpr(DigammaFunctor());
		Eigen::VectorXd tri_a = a.unaryExpr(TrigammaFunctor());
		Eigen::VectorXd tri_b = b.unaryExpr(TrigammaFunctor());

		double* H_data = H.data();
		for (int i = 0; i < m_n; ++i) {
			double mui = mu[i];
			double dmu = mui * (1.0 - mui);
			double d2mu = dmu * (1.0 - 2.0 * mui);
			double C = dig_a[i] - dig_b[i] - m_log_y[i] + m_log1_y[i];
			double B = phi * C;
			double B_mu = phi * phi * (tri_a[i] + tri_b[i]);
			double obs_weight = m_weights[i];
			double w_beta = obs_weight * (B_mu * dmu * dmu + B * d2mu);
			const double* xi = m_X.data() + i;  // xi[j * m_n] == X(i,j)

			for (int c = 0; c < m_p; ++c) {
				const double wxi_c = w_beta * xi[c * m_n];
				for (int r = 0; r <= c; ++r)
					H_data[r + c * total_p] += wxi_c * xi[r * m_n];
			}

			double B_log_phi = obs_weight * phi * (C + a[i] * tri_a[i] - b[i] * tri_b[i]);
			const double s = B_log_phi * dmu;
			for (int r = 0; r < m_p; ++r)
				H_data[r + m_p * total_p] += s * xi[r * m_n];
		}

		double D = -m_weight_sum * fast_digamma(phi);
		double D_phi = -m_weight_sum * R::trigamma(phi);
		for (int i = 0; i < m_n; ++i) {
			double mui = mu[i];
			double obs_weight = m_weights[i];
			D += obs_weight * (
				mui * dig_a[i] + (1.0 - mui) * dig_b[i] -
				mui * m_log_y[i] - (1.0 - mui) * m_log1_y[i]
			);
			D_phi += obs_weight * (mui * mui * tri_a[i] + (1.0 - mui) * (1.0 - mui) * tri_b[i]);
		}
		H(m_p, m_p) = phi * D + phi * phi * D_phi;
		for (int c = 0; c < total_p; ++c)
			for (int r = 0; r < c; ++r)
				H_data[c + r * total_p] = H_data[r + c * total_p];
		return H;
	}

	Eigen::MatrixXd expected_hessian(const Eigen::VectorXd& params) {
		int total_p = m_p + 1;
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
		Eigen::VectorXd beta = params.head(m_p);
		double phi = std::exp(params[m_p]);
		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
		double epsilon = 1e-8;
		for (int i = 0; i < m_n; ++i) {
			if (mu[i] < epsilon) mu[i] = epsilon;
			if (mu[i] > 1.0 - epsilon) mu[i] = 1.0 - epsilon;
		}

		Eigen::VectorXd a = mu.array() * phi;
		Eigen::VectorXd b = (1.0 - mu.array()) * phi;
		Eigen::VectorXd tri_a = a.unaryExpr(TrigammaFunctor());
		Eigen::VectorXd tri_b = b.unaryExpr(TrigammaFunctor());

		const double trigamma_phi = R::trigamma(phi);
		double* H_data = H.data();
		for (int i = 0; i < m_n; ++i) {
			const double mui = mu[i];
			const double obs_weight = m_weights[i];
			const double dmu = mui * (1.0 - mui);
			const double w_beta = obs_weight * phi * phi * (tri_a[i] + tri_b[i]) * dmu * dmu;
			const double cross = obs_weight * phi * (a[i] * tri_a[i] - b[i] * tri_b[i]) * dmu;
			const double* xi = m_X.data() + i;

			for (int c = 0; c < m_p; ++c) {
				const double wxi_c = w_beta * xi[c * m_n];
				for (int r = 0; r <= c; ++r)
					H_data[r + c * total_p] += wxi_c * xi[r * m_n];
			}
			for (int r = 0; r < m_p; ++r)
				H_data[r + m_p * total_p] += cross * xi[r * m_n];
			H(m_p, m_p) += obs_weight * phi * phi * (
				-trigamma_phi + mui * mui * tri_a[i] + (1.0 - mui) * (1.0 - mui) * tri_b[i]
			);
		}

		for (int c = 0; c < total_p; ++c)
			for (int r = 0; r < c; ++r)
				H_data[c + r * total_p] = H_data[r + c * total_p];
		return H;
	}
};

ModelResult fast_beta_regression_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                        const Eigen::Ref<const Eigen::VectorXd>& y,
                                        const Eigen::VectorXd* weights = nullptr,
                                        const Eigen::VectorXd* warm_start_beta = nullptr,
                                        bool smart_cold_start = true,
                                        double start_phi = 10.0,
                                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                        std::string optimization_alg = "lbfgs",
                                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                        bool estimate_only = false) {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    if (warm_start_beta) {
        params.head(p) = *warm_start_beta;
    } else if (smart_cold_start) {
        // Smart warm_start_params: OLS on logit(y)
        Eigen::VectorXd y_logit = (y.array() / (1.0 - y.array())).log().matrix();
        // Handle potential INF/NA from 0/1 in y if any (though beta regr assumes (0,1))
        for(int i=0; i<y_logit.size(); ++i) {
            if (!std::isfinite(y_logit[i])) {
                double yi = std::max(1e-4, std::min(1.0 - 1e-4, y[i]));
                y_logit[i] = std::log(yi / (1.0 - yi));
            }
        }
        params.head(p) = safe_ols_solve(X, y_logit);
    } else {
        params.head(p).setZero();
    }
    params[p] = std::log(start_phi);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

    Eigen::VectorXd weights_work = weights == nullptr ? Eigen::VectorXd::Ones(X.rows()) : *weights;
    BetaRegression fun(y, X, weights_work);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, 1000, 1e-6, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // phi
    res.XtWX = estimate_only ? Eigen::MatrixXd::Zero(p+1, p+1) : fun.hessian(params);
    res.converged = fit.converged;
    return res;
}

} // namespace

//' @title Compute Beta Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a beta regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param params A numeric vector of parameters [beta, log_phi].
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_beta_regression_score_cpp(SEXP X_sexp,
                                              SEXP y_sexp,
                                              SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    BetaRegression fun(y, X);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad; // Return the actual score (gradient of log-likelihood)
}

//' @title Compute Beta Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a beta regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters [beta, log_phi].
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_beta_regression_hessian_cpp(SEXP X_sexp,
                                                SEXP y_sexp,
                                                SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    BetaRegression fun(y, X);
    return -fun.hessian(params); // Return the actual Hessian of log-likelihood (Fisher Information is -Hessian)
}

//' @title Fast Beta Regression (C++)
//' @description High-performance beta regression fitting using Newton-Raphson or L-BFGS.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param start_phi Optional starting value for precision parameter phi.
//' @param compute_std_errs Deprecated.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, phi, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = runif(10)
//' fast_beta_regression_cpp(X, y)
// [[Rcpp::export]]
List fast_beta_regression_cpp(SEXP X_sexp,
								SEXP y_sexp,
								Nullable<NumericVector> warm_start_beta = R_NilValue,
								bool smart_cold_start = true,
								double start_phi = 10.0,
                                bool compute_std_errs = false,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "lbfgs",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                bool estimate_only = false) {

    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (warm_start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(warm_start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y, nullptr, sb_ptr, smart_cold_start, start_phi, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);

    Eigen::VectorXd params_full(fit.b.size() + 1);
    params_full.head(fit.b.size()) = fit.b;
    params_full[fit.b.size()] = std::log(fit.dispersion);
    BetaRegression fun_neg_ll(y, X);
    Eigen::VectorXd dummy_grad(params_full.size());
    double neg_loglik = fun_neg_ll(params_full, dummy_grad);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("neg_loglik") = neg_loglik,
		Named("converged") = fit.converged,
		Named("fisher_information") = fit.XtWX
	);
}

//' @title Fast Weighted Beta Regression (C++)
//' @description High-performance beta regression fitting with nonnegative row weights.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param weights A nonnegative numeric vector of row weights.
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when no warm start is provided.
//' @param start_phi Optional starting value for precision parameter phi.
//' @param compute_std_errs Deprecated.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @param estimate_only If TRUE, skip Fisher information calculation.
//' @return A list containing coefficients, phi, negative log-likelihood, convergence status, and Fisher information.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_beta_regression_weighted_cpp(SEXP X_sexp,
								SEXP y_sexp,
								SEXP weights_sexp,
								Nullable<NumericVector> warm_start_beta = R_NilValue,
								bool smart_cold_start = true,
								double start_phi = 10.0,
                                bool compute_std_errs = false,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "lbfgs",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                bool estimate_only = false) {

    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector weights_r(weights_sexp);
    if (weights_r.size() != X_r.nrow()) {
        stop("weights length must equal nrow(X)");
    }
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(weights_r.begin(), weights_r.size());
    if ((weights.array() < 0.0).any() || !weights.allFinite() || weights.sum() <= 0.0) {
        stop("weights must be finite, nonnegative, and have positive sum");
    }
    Eigen::VectorXd weights_vec = weights;

    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (warm_start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(warm_start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y, &weights_vec, sb_ptr, smart_cold_start, start_phi, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);

    Eigen::VectorXd params_full(fit.b.size() + 1);
    params_full.head(fit.b.size()) = fit.b;
    params_full[fit.b.size()] = std::log(fit.dispersion);
    BetaRegression fun_neg_ll(y, X, weights);
    Eigen::VectorXd dummy_grad(params_full.size());
    double neg_loglik = fun_neg_ll(params_full, dummy_grad);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("neg_loglik") = neg_loglik,
		Named("converged") = fit.converged,
		Named("fisher_information") = fit.XtWX
	);
}

//' @title Fast Beta Regression with Variance (C++)
//' @description Beta regression with full variance-covariance matrix and standard error estimation.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param start_phi Optional starting value for precision parameter phi.
//' @param compute_std_errs Deprecated.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, phi, vcov, standard errors, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = runif(10)
//' fast_beta_regression_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_beta_regression_with_var_cpp(SEXP X_sexp,
									 SEXP y_sexp,
									 Nullable<NumericVector> warm_start_beta = R_NilValue,
									 bool smart_cold_start = true,
									 double start_phi = 10.0,

                                     bool compute_std_errs = true,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "lbfgs",
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {

    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (warm_start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(warm_start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y, nullptr, sb_ptr, smart_cold_start, start_phi, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols() + 1, fixed_idx, fixed_values);
    Eigen::MatrixXd H_free = subset_matrix(fit.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
	Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd cov_mat = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    Eigen::VectorXd se = cov_mat.diagonal().array().sqrt();

    Eigen::VectorXd params_full(fit.b.size() + 1);
    params_full.head(fit.b.size()) = fit.b;
    params_full[fit.b.size()] = std::log(fit.dispersion);
    BetaRegression fun_neg_ll(y, X);
    Eigen::VectorXd dummy_grad(params_full.size());
    double neg_loglik = fun_neg_ll(params_full, dummy_grad);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("neg_loglik") = neg_loglik,
		Named("vcov") = cov_mat,
		Named("std_errs") = se,
        Named("converged") = fit.converged,
        Named("fisher_information") = fit.XtWX
		);
	}
