#include "_helper_functions.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

class AdjacentCategoryLogitNegLogLik {
private:
    const Eigen::Ref<const MatrixXd> m_X;
    const std::vector<int>& m_y;
    int m_n;
    int m_p;
    int m_K;

public:
    AdjacentCategoryLogitNegLogLik(const Eigen::Ref<const MatrixXd>& X, const std::vector<int>& y, int K) :
        m_X(X), m_y(y), m_n(X.rows()), m_p(X.cols()), m_K(K) {}

    double operator()(const VectorXd& params, VectorXd& grad) const {
        const int n_alpha = m_K - 1;
        const VectorXd alpha = params.head(n_alpha);
        const VectorXd beta = params.tail(m_p);

        grad.setZero(params.size());

        double neg_ll = 0.0;
        std::vector<double> alpha_suffix(m_K, 0.0);
        for (int j = m_K - 2; j >= 0; --j) {
            alpha_suffix[j] = alpha_suffix[j + 1] + alpha[j];
        }

        VectorXd score_offsets(m_K);
        for (int j = 0; j < m_K - 1; ++j) {
            score_offsets[j] = -static_cast<double>(m_K - 1 - j);
        }
        score_offsets[m_K - 1] = 0.0;
        VectorXd scores = VectorXd::Zero(m_K);
        VectorXd probs = VectorXd::Zero(m_K);
        VectorXd cdf = VectorXd::Zero(m_K - 1);
        VectorXd y_levels = VectorXd::LinSpaced(m_K, 1.0, static_cast<double>(m_K));

        for (int i = 0; i < m_n; ++i) {
            const double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;

            for (int j = 0; j < m_K - 1; ++j) {
                scores[j] = alpha_suffix[j] - static_cast<double>(m_K - 1 - j) * eta;
            }
            scores[m_K - 1] = 0.0;

            const double max_score = scores.maxCoeff();
            probs = (scores.array() - max_score).exp().matrix();
            const double denom = probs.sum();
            probs /= denom;

            const int y_i = m_y[i];
            neg_ll -= (scores[y_i - 1] - max_score - std::log(denom));

            double running_cdf = 0.0;
            double ey = probs.dot(y_levels);
            for (int j = 0; j < m_K; ++j) {
                if (j < m_K - 1) {
                    running_cdf += probs[j];
                    cdf[j] = running_cdf;
                }
            }

            for (int j = 0; j < m_K - 1; ++j) {
                grad[j] -= (((y_i <= (j + 1)) ? 1.0 : 0.0) - cdf[j]);
            }
            if (m_p > 0) {
                grad.tail(m_p).noalias() -= m_X.row(i).transpose() * (static_cast<double>(y_i) - ey);
            }
        }

        return neg_ll;
    }

    MatrixXd hessian(const VectorXd& params) const {
        const int n_alpha = m_K - 1;
        const VectorXd beta = params.tail(m_p);
        
        MatrixXd hess = MatrixXd::Zero(params.size(), params.size());

        std::vector<double> alpha_suffix(m_K, 0.0);
        for (int j = m_K - 2; j >= 0; --j) {
            alpha_suffix[j] = alpha_suffix[j + 1] + params[j];
        }

        VectorXd y_levels = VectorXd::LinSpaced(m_K, 1.0, static_cast<double>(m_K));
        VectorXd scores = VectorXd::Zero(m_K);
        VectorXd probs = VectorXd::Zero(m_K);
        VectorXd cdf = VectorXd::Zero(m_K - 1);
        VectorXd prefix_first_moment = VectorXd::Zero(m_K - 1);
        const int total_p = params.size();
        double* H_data = hess.data();

        for (int i = 0; i < m_n; ++i) {
            const double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;

            for (int j = 0; j < m_K - 1; ++j) {
                scores[j] = alpha_suffix[j] - static_cast<double>(m_K - 1 - j) * eta;
            }
            scores[m_K - 1] = 0.0;

            const double max_score = scores.maxCoeff();
            probs = (scores.array() - max_score).exp().matrix();
            probs /= probs.sum();

            const double ey = probs.dot(y_levels);
            const double ey2 = probs.dot(y_levels.array().square().matrix());
            double running_cdf = 0.0;
            double running_first_moment = 0.0;
            for (int j = 0; j < m_K; ++j) {
                if (j < m_K - 1) {
                    running_cdf += probs[j];
                    running_first_moment += y_levels[j] * probs[j];
                    cdf[j] = running_cdf;
                    prefix_first_moment[j] = running_first_moment;
                }
            }
            double var_y = std::max(0.0, ey2 - ey * ey);

            for (int j = 0; j < m_K - 1; ++j) {
                for (int k = j; k < m_K - 1; ++k) {
                    double f_min = cdf[std::min(j, k)];
                    double val = (f_min - cdf[j] * cdf[k]);
                    hess(j, k) += val;
                    if (j != k) hess(k, j) += val;
                }
            }

            if (m_p > 0) {
                const double* xi = m_X.data() + i;  // xi[b * m_n] == X(i,b)
                // Cross block: alpha-beta and beta-alpha (both triangles filled directly)
                for (int j = 0; j < m_K - 1; ++j) {
                    double cov_ind_y = prefix_first_moment[j] - ey * cdf[j];
                    for (int b = 0; b < m_p; ++b) {
                        double val = cov_ind_y * xi[b * m_n];
                        H_data[j + (n_alpha + b) * total_p] += val;  // hess(j, n_alpha+b)
                        H_data[(n_alpha + b) + j * total_p] += val;  // hess(n_alpha+b, j)
                    }
                }
                // Beta-beta block (upper triangle only)
                for (int c = 0; c < m_p; ++c) {
                    const double s = var_y * xi[c * m_n];
                    for (int r = 0; r <= c; ++r)
                        H_data[(n_alpha + r) + (n_alpha + c) * total_p] += s * xi[r * m_n];
                }
            }
        }

        // Reflect upper triangle to lower for beta-beta block only
        if (m_p > 0) {
            for (int c = 0; c < m_p; ++c)
                for (int r = 0; r < c; ++r)
                    H_data[(n_alpha + c) + (n_alpha + r) * total_p] = H_data[(n_alpha + r) + (n_alpha + c) * total_p];
        }

        return hess;
    }
};

std::vector<double> get_levels(const Eigen::Ref<const VectorXd>& y) {
    std::vector<double> levels(y.data(), y.data() + y.size());
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
    return levels;
}

std::vector<int> map_y_to_1K(const Eigen::Ref<const VectorXd>& y, const std::vector<double>& levels) {
    int n = y.size();
    int K = levels.size();
    std::vector<int> y_mapped(n);
    for (int i = 0; i < n; ++i) {
        double yi = y[i];
        auto it = std::lower_bound(levels.begin(), levels.end(), yi);
        y_mapped[i] = static_cast<int>(std::distance(levels.begin(), it)) + 1;
    }
    return y_mapped;
}

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_adjacent_category_logit_score_cpp(SEXP X_sexp,
                                                      SEXP y_sexp,
                                                      SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    std::vector<double> levels = get_levels(y);
    std::vector<int> y_mapped = map_y_to_1K(y, levels);
    AdjacentCategoryLogitNegLogLik fun(X, y_mapped, levels.size());
    VectorXd grad(params.size());
    fun(params, grad);
    return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_adjacent_category_logit_hessian_cpp(SEXP X_sexp,
                                                        SEXP y_sexp,
                                                        SEXP params_sexp) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

    std::vector<double> levels = get_levels(y);
    std::vector<int> y_mapped = map_y_to_1K(y, levels);
    AdjacentCategoryLogitNegLogLik fun(X, y_mapped, levels.size());
    return -fun.hessian(params);
}

//' @title Fast Adjacent-Category Logit (C++)
//' @description High-performance adjacent-category logit model fitting.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (categorical).
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @return A list containing coefficients, alpha, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_adjacent_category_logit_cpp(SEXP X_sexp, SEXP y_sexp, int maxit = 100, double tol = 1e-8,
                                        bool smart_cold_start = true,
                                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                        std::string optimization_alg = "lbfgs",
                                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    std::vector<double> levels = get_levels(y);
    int K = levels.size();
    if (K < 2) {
        stop("Adjacent-category logits require at least two observed outcome categories.");
    }
    std::vector<int> y_mapped = map_y_to_1K(y, levels);
    AdjacentCategoryLogitNegLogLik fun(X, y_mapped, K);
    
    int n_alpha = K - 1;
    int p = X.cols();
    int n_par = n_alpha + p;
    VectorXd params = VectorXd::Zero(n_par);
    
    if (warm_start_params.isNotNull()) {
        params = as<VectorXd>(warm_start_params);
        if (params.size() != n_par) stop("warm_start_params size mismatch");
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == p) {
            params.tail(p) = sb;
        }
    } else if (smart_cold_start) {
        // Smart warm_start_params: OLS on y_mapped
        Eigen::VectorXd y_double(y_mapped.size());
        for(size_t i=0; i<y_mapped.size(); ++i) y_double[i] = (double)y_mapped[i];
        params.tail(p) = ols_smart_cold_start_beta(X, y_double);
    }
    
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_start_ptr);

    return List::create(
        Named("b") = fit.params.tail(X.cols()),
        Named("alpha") = fit.params.head(K - 1),
        Named("params") = fit.params,
        Named("converged") = fit.converged
    );
}

//' @title Fast Adjacent-Category Logit with Variance (C++)
//' @description Adjacent-category logit model fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (categorical).
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_adjacent_category_logit_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int maxit = 100, double tol = 1e-8,
                                                bool smart_cold_start = true,
                                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                std::string optimization_alg = "lbfgs",
                                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    std::vector<double> levels = get_levels(y);
    int K = levels.size();
    if (K < 2) {
        stop("Adjacent-category logits require at least two observed outcome categories.");
    }
    std::vector<int> y_mapped = map_y_to_1K(y, levels);
    AdjacentCategoryLogitNegLogLik fun(X, y_mapped, K);
    
    int n_alpha = K - 1;
    int p = X.cols();
    int n_par = n_alpha + p;
    VectorXd params = VectorXd::Zero(n_par);

    if (warm_start_params.isNotNull()) {
        params = as<VectorXd>(warm_start_params);
        if (params.size() != n_par) stop("warm_start_params size mismatch");
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == p) {
            params.tail(p) = sb;
        }
    } else if (smart_cold_start) {
        // Smart warm_start_params: OLS on y_mapped
        Eigen::VectorXd y_double(y_mapped.size());
        for(size_t i=0; i<y_mapped.size(); ++i) y_double[i] = (double)y_mapped[i];
        params.tail(p) = ols_smart_cold_start_beta(X, y_double);
    }

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_start_ptr);

    MatrixXd info = fun.hessian(fit.params);
    MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
    int free_j = -1;
    for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
        if (fixed_spec.free_idx[jj] == n_alpha) { free_j = jj; break; }

    // Adjacent-category fits can have estimable treatment effects even when
    // nuisance columns make the information matrix rank deficient.  LDLT may
    // report success for such matrices while returning a finite but invalid
    // (often negative) diagonal entry, so use the rank-aware inverse here.
    MatrixXd cov_free = symmetric_pseudo_inverse(info_free);
    double ssq_b_1 = NA_REAL;
    if (X.cols() >= 1 && free_j >= 0 && cov_free.allFinite()) {
        const double treatment_variance = cov_free(free_j, free_j);
        if (R_finite(treatment_variance) && treatment_variance > 0.0) {
            ssq_b_1 = treatment_variance;
        }
    }

    SEXP vcov_sexp = R_NilValue;
    if (fit.converged) {
        MatrixXd vcov = expand_free_covariance(n_par, fixed_spec, cov_free, true);
        vcov_sexp = Rcpp::wrap(vcov);
    }
    return List::create(
        Named("b") = fit.params.tail(X.cols()),
        Named("alpha") = fit.params.head(K - 1),
        Named("params") = fit.params,
        Named("ssq_b_1") = ssq_b_1,
        Named("ssq_b_j") = ssq_b_1,
        Named("vcov") = vcov_sexp,
        Named("fisher_information") = info,
        Named("converged") = fit.converged
    );
}

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
NumericVector compute_adj_cat_logit_distr_parallel_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    const Rcpp::IntegerMatrix& w_mat,
    double delta,
    int num_cores
) {
    NumericMatrix X_r(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    int nsim = w_mat.cols();
    int n = y.size();
    int p_covars = X.cols();
    int p_full = p_covars + 1;

    std::vector<double> results(nsim, NA_REAL);
    const int* w_ptr = w_mat.begin();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;

        Eigen::MatrixXd X_full(n, p_full);
        Eigen::VectorXd y_shifted(n);

        for (int i = 0; i < n; ++i) {
            X_full(i, 0) = (double)w_col[i];
            for (int k = 0; k < p_covars; ++k) {
                X_full(i, 1 + k) = X(i, k);
            }
            y_shifted[i] = (w_col[i] == 1) ? y[i] + delta : y[i];
        }

        std::vector<double> levels = get_levels(y_shifted);
        int K = levels.size();
        if (K < 2) continue;

        std::vector<int> y_mapped = map_y_to_1K(y_shifted, levels);
        AdjacentCategoryLogitNegLogLik fun(X_full, y_mapped, K);
        VectorXd params = VectorXd::Zero((K - 1) + p_full);

        LikelihoodFitResult fit = optimize_likelihood(fun, params, 100, 1e-8, "newton_raphson", "newton_raphson");

        int n_alpha = K - 1;
        if ((int)fit.params.size() >= n_alpha + 1 && std::isfinite(fit.params[n_alpha])) {
            results[b] = fit.params[n_alpha];
        }
    }

    return wrap(results);
}
