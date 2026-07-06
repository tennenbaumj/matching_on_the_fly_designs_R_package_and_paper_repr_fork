#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <unordered_map>

using namespace Rcpp;

namespace {

double smart_negbin_theta_start_from_beta(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                          const Eigen::Ref<const Eigen::VectorXi>& y,
                                          const Eigen::Ref<const Eigen::VectorXd>& beta,
                                          double legacy_theta_start) {
    const int n = X.rows();
    const int p = X.cols();
    const int df = std::max(1, n - p);
    if (beta.size() != p || !beta.allFinite()) return legacy_theta_start;

    Eigen::VectorXd eta = (X * beta).array().min(700.0).matrix();
    Eigen::VectorXd mu = eta.array().exp().max(1e-8).matrix();
    double alpha_sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double yi = static_cast<double>(y[i]);
        const double mui = mu[i];
        alpha_sum += ((yi - mui) * (yi - mui) - mui) / (mui * mui);
    }

    const double alpha_hat = alpha_sum / static_cast<double>(df);
    if (!std::isfinite(alpha_hat) || alpha_hat <= 0.0) return legacy_theta_start;

    const double theta_hat = 1.0 / alpha_hat;
    if (!std::isfinite(theta_hat) || theta_hat <= 0.0) return legacy_theta_start;
    return std::max(0.1, theta_hat);
}

class NBLogLik {
private:
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXi> m_y;
    const int m_n;
    const int m_p;

    // Per-observation slot into distinct-y tables (built once at construction).
    std::vector<int>    m_y_slot;
    std::vector<double> m_distinct_y;
    std::vector<double> m_lgamma_y1;        // lgamma(y+1) = log(y!) per distinct y
    std::vector<double> m_lgamma_yptheta;   // preallocated; filled per operator() call
    std::vector<double> m_digamma_yptheta;  // preallocated; filled per operator() and hessian() call
    std::vector<double> m_trigamma_yptheta; // preallocated; filled per hessian() call
    Eigen::VectorXd     m_coef_vec;         // per-obs gradient scalar; GEMV after obs loop

public:
    NBLogLik(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXi>& y) :
        m_X(X), m_y(y), m_n(X.rows()), m_p(X.cols()), m_y_slot(X.rows(), -1),
        m_coef_vec(X.rows())
    {
        std::unordered_map<int, int> seen;
        for (int i = 0; i < m_n; ++i) {
            const int yi = m_y[i];
            auto it = seen.find(yi);
            if (it == seen.end()) {
                const int slot = static_cast<int>(m_distinct_y.size());
                seen[yi] = slot;
                m_distinct_y.push_back(static_cast<double>(yi));
                m_lgamma_y1.push_back(std::lgamma(static_cast<double>(yi) + 1.0));
                m_y_slot[i] = slot;
            } else {
                m_y_slot[i] = it->second;
            }
        }
        const int nd = static_cast<int>(m_distinct_y.size());
        m_lgamma_yptheta.resize(nd);
        m_digamma_yptheta.resize(nd);
        m_trigamma_yptheta.resize(nd);
    }

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        const Eigen::VectorXd beta = params.head(m_p);
        const double log_theta = params[m_p];
        const double theta = std::exp(log_theta);

        const Eigen::VectorXd eta = m_X * beta;

        // Hoist theta-only constants out of the observation loop.
        const double digamma_theta = fast_digamma(theta);
        const double lgamma_theta  = std::lgamma(theta);
        const double log_theta_val = std::log(theta);

        // Fill preallocated per-distinct-y tables for lgamma(y+theta) and digamma(y+theta).
        const int nd = static_cast<int>(m_distinct_y.size());
        for (int k = 0; k < nd; ++k) {
            const double ypt = m_distinct_y[k] + theta;
            m_lgamma_yptheta[k]  = std::lgamma(ypt);
            m_digamma_yptheta[k] = fast_digamma(ypt);
        }

        double neg_ll = 0.0;
        double score_log_theta = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const double eta_i  = eta[i];
            const double mu_i   = std::exp(eta_i);
            const double yi     = m_distinct_y[m_y_slot[i]];
            const double denom  = theta + mu_i;
            const double log_denom = std::log(denom);
            const int slot = m_y_slot[i];

            // log dnbinom_mu via explicit formula — avoids R::dnbinom_mu overhead.
            neg_ll -= m_lgamma_yptheta[slot] - lgamma_theta - m_lgamma_y1[slot]
                    + theta * (log_theta_val - log_denom)
                    + yi * (eta_i - log_denom);

            m_coef_vec[i] = yi - mu_i * (yi + theta) / denom;

            const double dlogf = m_digamma_yptheta[slot] - digamma_theta
                               + log_theta_val - log_denom + 1.0 - (yi + theta) / denom;
            score_log_theta += theta * dlogf;
        }
        grad.head(m_p).noalias() = -(m_X.transpose() * m_coef_vec);
        grad[m_p] = -score_log_theta;
        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double theta = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;

        const double digamma_theta = fast_digamma(theta);
        const double trigamma_theta = R::trigamma(theta);
        const double log_theta = std::log(theta);
        const int nd = static_cast<int>(m_distinct_y.size());
        for (int k = 0; k < nd; ++k) {
            const double ypt = m_distinct_y[k] + theta;
            m_digamma_yptheta[k] = fast_digamma(ypt);
            m_trigamma_yptheta[k] = R::trigamma(ypt);
        }

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            const int slot = m_y_slot[i];
            double yi = m_distinct_y[slot];
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;  // xi[j * m_n] == X(i,j)

            double w_beta = mu_i * theta * (yi + theta) / (denom * denom);
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            double d_score_beta_d_log_theta = theta * mu_i * (yi - mu_i) / (denom * denom);
            for (int r = 0; r < m_p; ++r)
                H_data[r + m_p * total_p] -= d_score_beta_d_log_theta * xi[r * m_n];

            double A = m_digamma_yptheta[slot] - digamma_theta +
                log_theta - std::log(denom) + 1.0 - (yi + theta) / denom;
            double dA_dtheta = m_trigamma_yptheta[slot] - trigamma_theta +
                1.0 / theta - 1.0 / denom + (yi - mu_i) / (denom * denom);
            H(m_p, m_p) -= theta * A + theta * theta * dA_dtheta;
        }
        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

    Eigen::MatrixXd expected_hessian(const Eigen::VectorXd& params) {
        int total_p = m_p + 1;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta = params.head(m_p);
        double theta = std::exp(params[m_p]);
        Eigen::VectorXd eta = m_X * beta;
        const double trigamma_theta = R::trigamma(theta);

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double mu_i = std::exp(eta[i]);
            double denom = theta + mu_i;
            const double* xi = m_X.data() + i;

            double w_beta = mu_i * theta / denom;
            for (int c = 0; c < m_p; ++c) {
                const double wxi_c = w_beta * xi[c * m_n];
                for (int r = 0; r <= c; ++r)
                    H_data[r + c * total_p] += wxi_c * xi[r * m_n];
            }

            const double e_trigamma = expected_trigamma_y_plus_theta(mu_i, theta, trigamma_theta);
            H(m_p, m_p) += -theta * theta * (
                e_trigamma - trigamma_theta + 1.0 / theta - 1.0 / denom
            );
        }

        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

private:
    static double expected_trigamma_y_plus_theta(double mu,
                                                  double theta,
                                                  double trigamma_theta) {
        const double prob_success = theta / (theta + mu);
        double pk = std::exp(theta * std::log(prob_success));
        double trigamma_yptheta = trigamma_theta;
        double sum = pk * trigamma_yptheta;
        double cdf = pk;
        const double ratio_base = mu / (theta + mu);
        const double mean = mu;
        const double sd = std::sqrt(mu + mu * mu / theta);
        const int min_iter = static_cast<int>(std::ceil(mean + 10.0 * sd));
        const int max_iter = 100000;

        for (int k = 0; k < max_iter; ++k) {
            pk *= (static_cast<double>(k) + theta) / static_cast<double>(k + 1) * ratio_base;
            const int y = k + 1;
            const double x = static_cast<double>(k) + theta;
            trigamma_yptheta -= 1.0 / (x * x);
            sum += pk * trigamma_yptheta;
            cdf += pk;
            if (y > min_iter && pk < 1e-14 && 1.0 - cdf < 1e-12) break;
        }
        return sum;
    }
};

ModelResult fast_neg_bin_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                  const Eigen::Ref<const Eigen::VectorXi>& y,
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                  bool smart_cold_start = true,
                                  int maxit = 1000,
                                  double eps_g = 1e-6,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "lbfgs",
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                  bool estimate_only = false) {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params = Eigen::VectorXd::Zero(p + 1);
    const Eigen::VectorXd y_double = y.cast<double>();
    const double mean_y = y_double.mean();
    const double var_y = (y_double.array() - mean_y).square().sum() /
        static_cast<double>(std::max(1, static_cast<int>(y.size()) - 1));
    const double theta_start = (var_y > mean_y && mean_y > 0.0) ?
        std::max(0.1, mean_y * mean_y / (var_y - mean_y)) : 10.0;
    Eigen::VectorXd legacy_params = Eigen::VectorXd::Zero(p + 1);
    if (p > 0 && X.col(0).array().isApprox(Eigen::ArrayXd::Ones(X.rows()), 1e-12)) {
        legacy_params[0] = std::log(std::max(mean_y, 1e-8));
    }
    legacy_params[p] = std::log(theta_start);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (params.size() != p + 1) stop("warm_start_params must have length equal to the number of model parameters");
    } else if (smart_cold_start) {
        Eigen::VectorXd beta_smart = ols_smart_cold_start_beta_on_log1p(X, y_double);
        Eigen::VectorXd beta_start = vector_is_usable_start(beta_smart, p) ? beta_smart : legacy_params.head(p);
        params.head(p) = beta_start;
        params[p] = std::log(smart_negbin_theta_start_from_beta(X, y, beta_start, theta_start));
    } else {
        params = legacy_params;
    }
    params = apply_fixed_values(params, fixed_spec);

    NBLogLik fun(X, y);

    Eigen::MatrixXd H_start_val;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start_val = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start_val;
    } else if (smart_cold_start) {
        H_start_val = edi_opt::negbin_smart_hessian(X, params, y);
        h_ptr = &H_start_val;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // theta
    res.XtWX = estimate_only ? Eigen::MatrixXd::Zero(p+1, p+1) : fun.hessian(params);
    res.iterations = fit.niter;
    res.converged = fit.converged;
    res.sigma2_hat = -fit.value; // using sigma2_hat to store logLik temporarily
    return res;
}

} // namespace

//' @title Compute Negative Binomial Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a negative binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (non-negative integers).
//' @param params A numeric vector of parameters [beta, log_theta].
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_negbin_regression_score_cpp(SEXP X_sexp,
                                                SEXP y_sexp,
                                                SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    NBLogLik fun(X, y);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

//' @title Compute Negative Binomial Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a negative binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters [beta, log_theta].
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_negbin_regression_hessian_cpp(SEXP X_sexp,
                                                  SEXP y_sexp,
                                                  SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    NBLogLik fun(X, y);
    return -fun.hessian(params);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_negbin_regression_expected_hessian_cpp(SEXP X_sexp,
                                                           SEXP y_sexp,
                                                           SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    NBLogLik fun(X, y);
    return fun.expected_hessian(params);
}

//' @title Fast Negative Binomial Regression with Variance (C++)
//' @description Negative binomial regression fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (non-negative integers).
//' @param warm_start_params Optional starting values for coefficients and dispersion. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, theta, vcov, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_neg_bin_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_neg_bin_with_var_cpp(SEXP X_sexp,
                                SEXP y_sexp,
                                Nullable<NumericVector> warm_start_params = R_NilValue,
                                bool smart_cold_start = false,
                                int maxit = 1000,
                                double eps_f = 1e-8,
                                double eps_g = 1e-6,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "lbfgs",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                bool estimate_only = false) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());

    ModelResult res = fast_neg_bin_internal(X, y, warm_start_params, smart_cold_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols() + 1, fixed_idx, fixed_values);
    Eigen::MatrixXd H_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("hess_fisher_info_matrix") = res.XtWX,
        Named("vcov") = vcov
    );
}

//' @title Fast Negative Binomial Regression (C++)
//' @description High-performance negative binomial regression fitting.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_params Optional starting values for coefficients and dispersion. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param eps_f Convergence tolerance for function value.
//' @param eps_g Convergence tolerance for gradient.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, theta, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_neg_bin_cpp(X, y)
// [[Rcpp::export]]
List fast_neg_bin_cpp(SEXP X_sexp,
                        SEXP y_sexp,
                        Nullable<NumericVector> warm_start_params = R_NilValue,
                        bool smart_cold_start = false,
                        int maxit = 1000,
                        double eps_f = 1e-8,
                        double eps_g = 1e-6,
                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                        std::string optimization_alg = "lbfgs",
                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                        bool estimate_only = false) {
    NumericMatrix X_r(X_sexp);
    IntegerVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());

    ModelResult res = fast_neg_bin_internal(X, y, warm_start_params, smart_cold_start, maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
    return List::create(
        Named("b") = res.b,
        Named("theta_hat") = res.dispersion,
        Named("logLik") = res.sigma2_hat,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX
    );
}
