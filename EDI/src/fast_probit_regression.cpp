#include "_helper_functions.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>

using namespace Rcpp;

namespace {

// Log-scale pnorm: falls back to R for |x| > 6 where direct log(erfc) loses precision.
inline double log_pnorm_lower(double x) {
    if (x > -6.0 && x < 6.0) return std::log(0.5 * fast_erfc(-x * kSqrt1_2));
    return R::pnorm5(x, 0.0, 1.0, 1, 1);
}
inline double log_pnorm_upper(double x) {
    if (x > -6.0 && x < 6.0) return std::log(0.5 * fast_erfc(x * kSqrt1_2));
    return R::pnorm5(x, 0.0, 1.0, 0, 1);
}

// Generalized residual for Probit: y*phi/Phi - (1-y)*phi/(1-Phi)
inline double probit_gen_residual_optimized(double y, double phi, double Phi, double Phi_inv) {
    if (y > 0.5) {
        return phi / Phi;
    } else {
        return -phi / Phi_inv;
    }
}

template<typename RDerived, typename WDerived>
inline void score_weighted_crossprod_colwise_assign(const Eigen::MatrixXd& X,
                                                    const Eigen::MatrixBase<RDerived>& residual,
                                                    const Eigen::MatrixBase<WDerived>& w,
                                                    Eigen::VectorXd& score,
                                                    Eigen::MatrixXd& out) {
    const int n = X.rows();
    const int p = X.cols();
    score.setZero();
    out.setZero();
    for (int j = 0; j < p; ++j) {
        const double* xj = X.col(j).data();
        for (int k = j; k < p; ++k) {
            const double* xk = X.col(k).data();
            double acc = 0.0;
            if (k == j) {
                double score_acc = 0.0;
                for (int i = 0; i < n; ++i) {
                    acc += xj[i] * w[i] * xj[i];
                    score_acc += xj[i] * residual[i];
                }
                score[j] = score_acc;
            } else {
                for (int i = 0; i < n; ++i) acc += xj[i] * w[i] * xk[i];
            }
            out(j, k) = acc;
            if (k != j) out(k, j) = acc;
        }
    }
}

class ProbitLbfgsObjective : public Numer::MFuncGrad {
private:
    const Eigen::Ref<const RowMajorMatrixXd> m_X;
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::VectorXd> m_weights;
    const Eigen::Ref<const Eigen::VectorXd> m_eta_fixed;
    bool m_use_weights;
    int m_n;

public:
    ProbitLbfgsObjective(const Eigen::Ref<const RowMajorMatrixXd>& X,
                         const Eigen::Ref<const Eigen::VectorXd>& y,
                         const Eigen::Ref<const Eigen::VectorXd>& weights,
                         const Eigen::Ref<const Eigen::VectorXd>& eta_fixed,
                         bool use_weights) :
        m_X(X), m_y(y), m_weights(weights), m_eta_fixed(eta_fixed), 
        m_use_weights(use_weights), m_n(X.rows()) {}

    virtual double f_grad(Numer::Constvec& beta, Numer::Refvec grad) override {
        const Eigen::VectorXd eta = m_X * beta + m_eta_fixed;
        Eigen::VectorXd gen_res(m_n);
        double f = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const double ei = eta[i];
            const double wi = m_use_weights ? m_weights[i] : 1.0;
            
            const double lp = log_pnorm_lower(ei);
            const double lq = log_pnorm_upper(ei);
            f -= wi * (m_y[i] * lp + (1.0 - m_y[i]) * lq);
            
            const double phi = dnorm_fast(ei);
            const double Phi = pnorm_fast(ei);
            const double Phi_inv = 1.0 - Phi;
            gen_res[i] = wi * probit_gen_residual_optimized(m_y[i], phi, Phi, Phi_inv);
        }

        grad.noalias() = -m_X.transpose() * gen_res;
        return f;
    }
};

} // namespace

// Internal probit fitting core.
ModelResult fast_probit_regression_internal(
        const Eigen::Ref<const Eigen::MatrixXd>& X_eigen,
        const Eigen::Ref<const Eigen::VectorXd>& y_eigen,
        const Eigen::Ref<const Eigen::VectorXd>& weights_eigen,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100,
        double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
        bool estimate_only = false) {

    const int n = X_eigen.rows();
    const int p = X_eigen.cols();
    const bool use_weights = (weights_eigen.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = static_cast<int>(fixed_spec.free_idx.size());

    Eigen::VectorXd beta_start = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (static_cast<int>(beta_start.size()) != p)
            Rcpp::stop("warm_start_beta must have length equal to ncol(X)");
    } else if (smart_cold_start) {
        beta_start = ols_smart_cold_start_beta(X_eigen, y_eigen.array().unaryExpr([](double v){
            return R::qnorm5((v + 0.5) / 2.0, 0.0, 1.0, 1, 0);
        }));
    }
    beta_start = apply_fixed_values(beta_start, fixed_spec);

    Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
    for (int k = 0; k < static_cast<int>(fixed_spec.fixed_idx.size()); ++k) {
        eta_fixed += X_eigen.col(fixed_spec.fixed_idx[k]) * fixed_spec.fixed_values[k];
    }

    Eigen::MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X_eigen.col(fixed_spec.free_idx[j]);

    Eigen::VectorXd beta_free = subset_vector(beta_start, fixed_spec.free_idx);

    ModelResult res;
    res.iterations = 0;
    res.converged = false;

    if (optimization_alg == "lbfgs") {
        double fopt = 0.0;
        ProbitLbfgsObjective obj(X_free, y_eigen, weights_eigen, eta_fixed, use_weights);
        int status = Numer::optim_lbfgs(obj, beta_free, fopt, maxit, tol, tol);
        res.converged = (status >= 0);
        res.iterations = NA_INTEGER;
        res.neg_ll = fopt;
    } else {
        Eigen::VectorXd mu(n);
        Eigen::VectorXd w(n);
        Eigen::VectorXd gen_res(n);
        Eigen::MatrixXd XtWX(p_free, p_free);
        Eigen::VectorXd score_free(p_free);

        for (int iter = 0; iter < maxit; ++iter) {
            res.iterations++;
            const Eigen::VectorXd eta = X_free * beta_free + eta_fixed;

            for (int i = 0; i < n; ++i) {
                const double ei = eta[i];
                const double wi = use_weights ? weights_eigen[i] : 1.0;
                
                const double phi = dnorm_fast(ei);
                const double Phi = pnorm_fast(ei);
                const double Phi_inv = 1.0 - Phi;
                const double vm = std::max(1e-15, Phi * Phi_inv);
                
                mu[i] = Phi;
                w[i] = wi * (phi * phi / vm);
                gen_res[i] = wi * probit_gen_residual_optimized(y_eigen[i], phi, Phi, Phi_inv);
            }

            const bool use_warm_xtwx = (iter == 0) && warm_start_fisher_info.isNotNull();
            if (!use_warm_xtwx) {
                score_weighted_crossprod_colwise_assign(X_free, gen_res, w, score_free, XtWX);
            } else {
                score_free.noalias() = X_free.transpose() * gen_res;
                Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(Rcpp::NumericMatrix(warm_start_fisher_info));
                XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
            }
            if (score_free.norm() < tol) { res.converged = true; break; }

            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
            if (ldlt.info() != Eigen::Success) break;
            const Eigen::VectorXd delta = ldlt.solve(score_free);
            if (!delta.allFinite()) break;

            beta_free += delta;
            if (delta.norm() < tol) { res.converged = true; break; }
        }
        res.mu = mu;
        res.XtWX = expand_free_covariance(p, fixed_spec, XtWX, false);
        res.score = expand_free_params(score_free, Eigen::VectorXd::Zero(p), fixed_spec);
        
        // Reuse mu[] from last IRLS iteration — avoids 1 GEMV + n erfc calls
        double nl = 0.0;
        for (int i = 0; i < n; ++i) {
            const double wi = use_weights ? weights_eigen[i] : 1.0;
            nl -= wi * (y_eigen[i] > 0.5 ? std::log(mu[i]) : std::log1p(-mu[i]));
        }
        res.neg_ll = nl;
    }

    res.b = expand_free_params(beta_free, beta_start, fixed_spec);
    
    if (!estimate_only) {
        const Eigen::VectorXd eta = X_eigen * res.b;
        if (res.mu.size() == 0) {
            res.mu = eta.array().unaryExpr([](double e){ return pnorm_fast(e); }).matrix();
        }
        if (res.XtWX.size() == 0) {
            Eigen::VectorXd weights_vec(n);
            for (int i = 0; i < n; ++i) {
                const double wi = use_weights ? weights_eigen[i] : 1.0;
                const double ei = eta[i];
                const double phi = dnorm_fast(ei);
                const double Phi = pnorm_fast(ei);
                const double vm = std::max(1e-15, Phi * (1.0 - Phi));
                weights_vec[i] = wi * (phi * phi / vm);
            }
            res.XtWX = X_eigen.transpose() * weights_vec.asDiagonal() * X_eigen;
        }
        if (res.score.size() == 0) {
            Eigen::VectorXd full_gen_res(n);
            for (int i = 0; i < n; ++i) {
                const double wi = use_weights ? weights_eigen[i] : 1.0;
                const double ei = eta[i];
                const double phi = dnorm_fast(ei);
                const double Phi = pnorm_fast(ei);
                full_gen_res[i] = wi * probit_gen_residual_optimized(y_eigen[i], phi, Phi, 1.0 - Phi);
            }
            res.score = X_eigen.transpose() * full_gen_res;
        }
    }
    
    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd get_probit_regression_score_cpp(SEXP X_sexp, SEXP y_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const Eigen::VectorXd eta = X * beta;
    const int n = X.rows();
    Eigen::VectorXd gen_res(n);
    for (int i = 0; i < n; ++i) {
        const double ei = eta[i];
        gen_res[i] = probit_gen_residual_optimized(y[i], dnorm_fast(ei), pnorm_fast(ei), pnorm_fast(-ei));
    }
    return X.transpose() * gen_res;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_probit_regression_hessian_cpp(SEXP X_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp);
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());

    const Eigen::VectorXd eta = X * beta;
    const int n = X.rows();
    Eigen::VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        const double ei = eta[i];
        const double phi = dnorm_fast(ei);
        const double Phi = pnorm_fast(ei);
        const double vm = std::max(1e-15, Phi * (1.0 - Phi));
        w[i] = phi * phi / vm;
    }
    return -weighted_crossprod(X, w);
}

//' @title Fast Probit Regression (C++)
//' @description High-performance probit GLM fitting via IRLS or L-BFGS.
//' @param X_sexp A numeric matrix of predictors (including intercept column).
//' @param y_sexp A numeric vector of binary responses (0/1).
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use an OLS-based initial guess when no warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm ("irls" or "lbfgs"). Default "irls".
//' @param warm_start_weights Optional initial IRLS working weights.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @param estimate_only Logical. If TRUE, skip variance computation and return only coefficients.
//' @return A list containing coefficients and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_probit_regression_cpp(SEXP X_sexp, SEXP y_sexp,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100, double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
        bool estimate_only = false) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    ModelResult res = fast_probit_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta,
        smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info, estimate_only);
    if (estimate_only) {
        return List::create(
            Named("b") = res.b,
            Named("converged") = res.converged,
            Named("iterations") = res.iterations
        );
    }
    const int n = X.rows();
    const Eigen::VectorXd eta = X * res.b;
    Eigen::VectorXd weights_vec(n);
    for (int i = 0; i < n; ++i) {
        const double ei = eta[i];
        const double phi = dnorm_fast(ei);
        const double Phi = pnorm_fast(ei);
        const double vm = std::max(1e-15, Phi * (1.0 - Phi));
        weights_vec[i] = phi * phi / vm;
    }
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
        Named("converged") = res.converged
    );
}

// [[Rcpp::export]]
List fast_probit_regression_weighted_cpp(SEXP X_sexp, SEXP y_sexp,
        SEXP weights_sexp,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100, double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericVector w_r(weights_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(w_r.begin(), w_r.size());

    ModelResult res = fast_probit_regression_internal(X, y, weights, warm_start_beta,
        smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}

// [[Rcpp::export]]
List fast_probit_regression_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int j = 2,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    ModelResult res = fast_probit_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta,
        smart_cold_start, 100, 1e-8, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info, false);
    
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int overall_j) -> int {
        for (int jj = 0; jj < static_cast<int>(fixed_spec.free_idx.size()); ++jj)
            if (fixed_spec.free_idx[jj] == overall_j) return jj + 1; // 1-based
        return -1;
    };

    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;

    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("params") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("score") = res.score,
        Named("observed_information") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("information") = res.XtWX,
        Named("information_type") = "fisher",
        Named("hessian") = -res.XtWX,
        Named("neg_loglik") = res.neg_ll,
        Named("neg_ll") = res.neg_ll,
        Named("loglik") = R_finite(res.neg_ll) ? -res.neg_ll : NA_REAL,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}
