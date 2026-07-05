#include "_helper_functions.h"
#include "ordinal_fixed_link_helpers.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

inline double plogis_stable_cpp(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

class OrdinalRegression {
private:
    edi_ordinal::FixedOrdinalRegression m_model;

public:
    OrdinalRegression(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXd>& y, const Eigen::Ref<const Eigen::VectorXd>& weights = Eigen::VectorXd()) :
        m_model(X, y, edi_ordinal::Link::Logit, -1.0, weights) {}

    static std::vector<double> init_levels(const Eigen::Ref<const Eigen::VectorXd>& y) {
        return edi_ordinal::init_levels(y);
    }

    double neg_log_likelihood(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.neg_log_likelihood(params);
    }

    double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad) const {
        return m_model(params, grad);
    }

    MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.hessian(params);
    }

    MatrixXd expected_hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.expected_hessian(params);
    }
};

} // namespace

//' @title Compute Ordinal Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for an ordinal regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (ordinal categories 1, 2, ...).
//' @param params A numeric vector of parameters [alpha, beta].
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ordinal_regression_score_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, const Rcpp::NumericVector& params) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());
    
    OrdinalRegression model(map_X, map_y);
    Eigen::VectorXd grad(map_params.size());
    model(map_params, grad);
    return wrap(-grad);
}

//' @title Compute Ordinal Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for an ordinal regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters [alpha, beta].
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ordinal_regression_hessian_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, const Rcpp::NumericVector& params) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());
    
    OrdinalRegression model(map_X, map_y);
    return wrap(-model.hessian(map_params));
}

// Simple solver using Newton-Raphson as we have a small number of parameters (usually)
//' @title Fast Ordinal Regression (C++)
//' @description High-performance ordinal regression fitting using Newton-Raphson.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_params Optional starting values for [alpha, beta]. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients (beta), thresholds (alpha), and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_regression_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, Nullable<NumericVector> warm_start_params = R_NilValue, bool smart_cold_start = true, int maxit = 100, double tol = 1e-6,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "lbfgs",
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                  bool estimate_only = false) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    OrdinalRegression model(map_X, map_y);
    int p = map_X.cols();
    int K = OrdinalRegression::init_levels(map_y).size();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (params.size() != n_params) stop("warm_start_params must have length equal to the number of model parameters");
    } else {
        OrdinalStart legacy_start;
        legacy_start.alpha = VectorXd(n_alpha);
        for (int k = 0; k < n_alpha; ++k) {
            legacy_start.alpha[k] = -1.0 + 2.0 * (k + 1) / K;
        }
        legacy_start.beta = VectorXd::Zero(p);
        
        if (smart_cold_start) {
            // Use OLS on rank-transformed y
            Eigen::VectorXd y_rank = (map_y.array() - 1.0) / static_cast<double>(K - 1);
            Eigen::VectorXd beta;
            if (edi_opt::robust_ols_solve(map_X, y_rank, beta)) {
                legacy_start.beta = beta;
            }
        }
        params = ordinal_start_to_params(legacy_start);
    }
    params = apply_fixed_values(params, fixed_spec);

    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    bool has_warm_fisher = warm_start_fisher_info.isNotNull();
    if (has_warm_fisher) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    } else if (smart_cold_start) {
        H_start = model.hessian(params);
        h_ptr = &H_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(model, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    if (estimate_only) {
        return List::create(
            Named("b") = params.tail(p),
            Named("alpha") = params.head(n_alpha),
            Named("n_params") = n_params,
            Named("params") = params,
            Named("converged") = fit.converged,
            Named("iterations") = fit.niter
        );
    }

    MatrixXd H = model.hessian(params);
    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params,
        Named("neg_loglik") = fit.value,
        Named("converged") = fit.converged,
        Named("iterations") = fit.niter,
        Named("observed_information") = H,
        Named("fisher_information") = H,
        Named("information") = H,
        Named("information_type") = "observed"
    );
}

//' @title Fast Weighted Ordinal Regression (C++)
//' @description High-performance weighted ordinal regression fitting using Newton-Raphson.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param weights A numeric vector of nonnegative observation weights.
//' @param warm_start_params Optional starting values for [alpha, beta].
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @return A list containing coefficients (beta), thresholds (alpha), and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_regression_weighted_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, const Rcpp::NumericVector& weights,
                                          Nullable<NumericVector> warm_start_params = R_NilValue, bool smart_cold_start = true, int maxit = 100, double tol = 1e-6,
                                          Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                          std::string optimization_alg = "lbfgs",
                                          Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_weights(weights.begin(), weights.size());

    if (map_weights.size() != map_X.rows()) {
        stop("weights must have length equal to nrow(X)");
    }
    OrdinalRegression model(map_X, map_y, map_weights);
    int p = map_X.cols();
    int K = OrdinalRegression::init_levels(map_y).size();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (params.size() != n_params) stop("warm_start_params must have length equal to the number of model parameters");
    } else {
        OrdinalStart legacy_start;
        legacy_start.alpha = VectorXd(n_alpha);
        for (int k = 0; k < n_alpha; ++k) {
            legacy_start.alpha[k] = -1.0 + 2.0 * (k + 1) / K;
        }
        legacy_start.beta = VectorXd::Zero(p);
        
        if (smart_cold_start) {
            // Use OLS on rank-transformed y
            Eigen::VectorXd y_rank = (map_y.array() - 1.0) / static_cast<double>(K - 1);
            Eigen::VectorXd beta;
            if (edi_opt::robust_ols_solve(map_X, y_rank, beta)) {
                legacy_start.beta = beta;
            }
        }
        params = ordinal_start_to_params(legacy_start);
    }
    params = apply_fixed_values(params, fixed_spec);

    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    } else if (smart_cold_start) {
        H_start = model.hessian(params);
        h_ptr = &H_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(model, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, h_ptr);
    params = fit.params;

    MatrixXd H = model.hessian(params);
    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params,
        Named("neg_loglik") = fit.value,
        Named("converged") = fit.converged,
        Named("iterations") = fit.niter,
        Named("observed_information") = H,
        Named("fisher_information") = H,
        Named("information") = H,
        Named("information_type") = "observed"
    );
}

//' @title Fast Ordinal Regression with Variance (C++)
//' @description Ordinal regression fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_params Optional starting values for [alpha, beta]. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, thresholds, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_regression_with_var_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y,
                                           Nullable<NumericVector> warm_start_params = R_NilValue,
                                           bool smart_cold_start = true,
                                           int maxit = 100, double tol = 1e-6,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "lbfgs",
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    List res = fast_ordinal_regression_cpp(X, y, warm_start_params, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
    VectorXd params = res["params"];
    bool converged = res["converged"];
    
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    MatrixXd H = as<MatrixXd>(res["observed_information"]);
    
    int n_params = params.size();
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
    FullPivLU<MatrixXd> lu(H_free);

    List output = List::create(
        Named("b") = res["b"],
        Named("alpha") = res["alpha"],
        Named("params") = params,
        Named("neg_loglik") = res["neg_loglik"],
        Named("converged") = converged,
        Named("iterations") = res["iterations"],
        Named("fisher_information") = H
    );
    if (!lu.isInvertible()) {
        output["ssq_b_j"] = NA_REAL;
        output["vcov"] = R_NilValue;
        return output;
    }
    
    MatrixXd vcov_free = lu.inverse();
    MatrixXd vcov_full = expand_free_covariance(n_params, fixed_spec, vcov_free, true);

    int p = map_X.cols();
    int n_alpha = n_params - p;
    
    // We want the variance of the first covariate after alphas (which is often the treatment)
    // In our params, it's at index n_alpha
    double ssq_b_j = (p >= 1) ? vcov_full(n_alpha, n_alpha) : NA_REAL;

    output["vcov"] = vcov_full;
    output["ssq_b_j"] = ssq_b_j;
    return output;
}

//' @title Ordinal G-Computation Post-Fit (C++)
//' @description Performs G-computation for ordinal outcomes using a fitted model's parameters.
//' @param X_fit Matrix of predictors used in fit.
//' @param y Vector of responses.
//' @param coef_hat Estimated coefficients (beta).
//' @param alpha_hat Estimated thresholds (alpha).
//' @param j_treat 1-based index of treatment column.
//' @return A list containing G-computation results (means, difference, SE).
//' @export
//' @keywords internal
// [[Rcpp::export]]
List ordinal_gcomp_post_fit_cpp(const Rcpp::NumericMatrix& X_fit,
                                const Rcpp::NumericVector& y,
                                const Rcpp::NumericVector& coef_hat,
                                const Rcpp::NumericVector& alpha_hat,
                                int j_treat) {
    Eigen::Map<const Eigen::MatrixXd> map_X_fit(X_fit.begin(), X_fit.rows(), X_fit.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_coef_hat(coef_hat.begin(), coef_hat.size());
    Eigen::Map<const Eigen::VectorXd> map_alpha_hat(alpha_hat.begin(), alpha_hat.size());

    const int n = map_X_fit.rows();
    const int p = map_X_fit.cols();
    const int n_alpha = map_alpha_hat.size();
    const int j_treat0 = j_treat - 1;

    if (j_treat0 < 0 || j_treat0 >= p) {
        stop("treatment column index is out of bounds");
    }
    if (map_y.size() != n || map_coef_hat.size() != p) {
        stop("dimension mismatch in ordinal_gcomp_post_fit_cpp");
    }

    VectorXd params(n_alpha + p);
    params.head(n_alpha) = map_alpha_hat;
    params.tail(p) = map_coef_hat;

    OrdinalRegression model(map_X_fit, map_y);
    MatrixXd H = model.hessian(params);
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) {
        stop("failed to invert ordinal Hessian");
    }
    MatrixXd vcov_full = lu.inverse();
    vcov_full = 0.5 * (vcov_full + vcov_full.transpose());

    for (int j = 0; j < vcov_full.rows(); ++j) {
        for (int k = 0; k < vcov_full.cols(); ++k) {
            if (!R_finite(vcov_full(j, k))) {
                stop("non-finite ordinal covariance");
            }
        }
    }

    MatrixXd X1 = map_X_fit;
    MatrixXd X0 = map_X_fit;
    X1.col(j_treat0).setOnes();
    X0.col(j_treat0).setZero();
    VectorXd eta1 = X1 * map_coef_hat;
    VectorXd eta0 = X0 * map_coef_hat;

    auto compute_mean = [&](const VectorXd& eta_vec) {
        double total_mean = 0.0;
        for (int i = 0; i < n; ++i) {
            double mean_i = 1.0;
            for (int k = 0; k < n_alpha; ++k) {
                mean_i += plogis_stable_cpp(eta_vec[i] - map_alpha_hat[k]);
            }
            total_mean += mean_i;
        }
        return total_mean / static_cast<double>(n);
    };

    VectorXd eta1_scratch(n);
    VectorXd eta0_scratch(n);
    auto compute_md_from_alpha_eta = [&](const auto& alpha,
                                         const VectorXd& eta1_values,
                                         const VectorXd& eta0_values) {
        for (int k = 1; k < n_alpha; ++k) {
            if (alpha[k] <= alpha[k - 1]) return NA_REAL;
        }
        double mean1_loc = 0.0;
        double mean0_loc = 0.0;
        for (int i = 0; i < n; ++i) {
            double m1 = 1.0;
            double m0 = 1.0;
            for (int k = 0; k < n_alpha; ++k) {
                m1 += plogis_stable_cpp(eta1_values[i] - alpha[k]);
                m0 += plogis_stable_cpp(eta0_values[i] - alpha[k]);
            }
            mean1_loc += m1;
            mean0_loc += m0;
        }
        return (mean1_loc - mean0_loc) / static_cast<double>(n);
    };

    const double mean1 = compute_mean(eta1);
    const double mean0 = compute_mean(eta0);
    const double md = mean1 - mean0;

    const double h_step = 1e-5;
    VectorXd grad_md(params.size());
    VectorXd p_plus = params;
    VectorXd p_minus = params;
    for (int j = 0; j < params.size(); ++j) {
        p_plus[j] = params[j] + h_step;
        p_minus[j] = params[j] - h_step;
        double f_plus;
        double f_minus;
        if (j < n_alpha) {
            f_plus = compute_md_from_alpha_eta(p_plus.head(n_alpha), eta1, eta0);
            f_minus = compute_md_from_alpha_eta(p_minus.head(n_alpha), eta1, eta0);
        } else {
            const int beta_j = j - n_alpha;
            eta1_scratch.noalias() = eta1 + h_step * X1.col(beta_j);
            eta0_scratch.noalias() = eta0 + h_step * X0.col(beta_j);
            f_plus = compute_md_from_alpha_eta(params.head(n_alpha), eta1_scratch, eta0_scratch);
            eta1_scratch.noalias() = eta1 - h_step * X1.col(beta_j);
            eta0_scratch.noalias() = eta0 - h_step * X0.col(beta_j);
            f_minus = compute_md_from_alpha_eta(params.head(n_alpha), eta1_scratch, eta0_scratch);
        }
        if (!R_finite(f_plus) || !R_finite(f_minus)) {
            grad_md[j] = NA_REAL;
        } else {
            grad_md[j] = (f_plus - f_minus) / (2.0 * h_step);
        }
        p_plus[j] = params[j];
        p_minus[j] = params[j];
    }

    double se_md = NA_REAL;
    bool grad_is_finite = true;
    for (int j = 0; j < grad_md.size(); ++j) {
        if (!R_finite(grad_md[j])) {
            grad_is_finite = false;
            break;
        }
    }
    if (grad_is_finite) {
        const double var_md = (grad_md.transpose() * vcov_full * grad_md)(0, 0);
        if (R_finite(var_md) && var_md >= 0.0) {
            se_md = std::sqrt(var_md);
        }
    }

    MatrixXd vcov_beta = vcov_full.block(n_alpha, n_alpha, p, p);
    VectorXd std_err(p);
    VectorXd z_vals(p);
    for (int j = 0; j < p; ++j) {
        const double var_j = vcov_beta(j, j);
        std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
        z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? map_coef_hat[j] / std_err[j] : NA_REAL;
    }

    return List::create(
        Named("vcov") = vcov_beta,
        Named("std_err") = std_err,
        Named("z_vals") = z_vals,
        Named("mean1") = mean1,
        Named("mean0") = mean0,
        Named("md") = md,
        Named("se_md") = se_md
    );
}

//' @title Expand Continuation Ratio Data (C++)
//' @description Utility to expand ordinal data for continuation ratio regression.
//' @param y Vector of responses.
//' @param w Vector of treatment indicators.
//' @param strata Vector of strata.
//' @param K Number of categories.
//' @return A list with expanded y, w, and strata.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List expand_continuation_ratio_data_cpp(const Rcpp::IntegerVector& y, const Rcpp::IntegerVector& w, const Rcpp::IntegerVector& strata, int K) {
    Eigen::Map<const Eigen::VectorXi> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXi> map_w(w.begin(), w.size());
    Eigen::Map<const Eigen::VectorXi> map_strata(strata.begin(), strata.size());

    int n = map_y.size();
    int n_alpha = K - 1;
    int num_strata = map_strata.maxCoeff();
    
    std::vector<int> y_stack;
    std::vector<int> w_stack;
    std::vector<int> strata_stack;
    
    for (int i = 0; i < n; ++i) {
        int yi = map_y[i];
        int wi = map_w[i];
        int si = map_strata[i];
        
        for (int j = 1; j <= std::min(yi, n_alpha); ++j) {
            y_stack.push_back((yi == j) ? 1 : 0);
            w_stack.push_back(wi);
            strata_stack.push_back(si + (j - 1) * num_strata);
        }
    }
    
    return List::create(
        Named("y") = wrap(y_stack),
        Named("w") = wrap(w_stack),
        Named("strata") = wrap(strata_stack)
    );
}

//' @title Expand Adjacent Category Data (C++)
//' @description Utility to expand ordinal data for adjacent category logit regression.
//' @param y Vector of responses.
//' @param w Vector of treatment indicators.
//' @param strata Vector of strata.
//' @param K Number of categories.
//' @return A list with expanded y, w, and strata.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List expand_adjacent_category_data_cpp(const Rcpp::IntegerVector& y, const Rcpp::IntegerVector& w, const Rcpp::IntegerVector& strata, int K) {
    Eigen::Map<const Eigen::VectorXi> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXi> map_w(w.begin(), w.size());
    Eigen::Map<const Eigen::VectorXi> map_strata(strata.begin(), strata.size());

    int n = map_y.size();
    int n_alpha = K - 1;
    int num_strata = map_strata.maxCoeff();
    
    std::vector<int> y_stack;
    std::vector<int> w_stack;
    std::vector<int> strata_stack;
    
    for (int i = 0; i < n; ++i) {
        int yi = map_y[i];
        int wi = map_w[i];
        int si = map_strata[i];
        
        for (int j = 1; j <= n_alpha; ++j) {
            if (yi == j || yi == j + 1) {
                y_stack.push_back((yi == j + 1) ? 1 : 0);
                w_stack.push_back(wi);
                strata_stack.push_back(si + (j - 1) * num_strata);
            }
        }
    }
    
    return List::create(
        Named("y") = wrap(y_stack),
        Named("w") = wrap(w_stack),
        Named("strata") = wrap(strata_stack)
    );
}
