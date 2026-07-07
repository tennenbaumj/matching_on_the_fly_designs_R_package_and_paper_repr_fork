#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// Internal pure C++ logic
ModelResult fast_ols_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                              const Eigen::Ref<const Eigen::VectorXd>& y,
                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                              bool estimate_only = false) {
    const int n = X.rows();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();
    
    Eigen::VectorXd y_adj = y;
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }

    ModelResult res;
    res.b = Eigen::VectorXd::Zero(p);
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    if (p_free > 0) {
        Eigen::MatrixXd XtX_free(p_free, p_free);
        Eigen::VectorXd Xty_free(p_free);
        Eigen::VectorXd beta_free;

        if (p_free == p) {
            XtX_free = symmetric_crossprod(X);
            Xty_free.noalias() = X.transpose() * y_adj;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
            if (ldlt.info() == Eigen::Success) {
                beta_free = ldlt.solve(Xty_free);
            } else {
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
                beta_free = qr.solve(y_adj);
            }
        } else {
            Eigen::MatrixXd X_free(n, p_free);
            for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
            XtX_free = symmetric_crossprod(X_free);
            Xty_free.noalias() = X_free.transpose() * y_adj;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
            if (ldlt.info() == Eigen::Success) {
                beta_free = ldlt.solve(Xty_free);
            } else {
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_free);
                beta_free = qr.solve(y_adj);
            }
        }

        for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
        if (!estimate_only) res.XtWX = expand_free_covariance(p, fixed_spec, XtX_free, false);
    } else {
        if (!estimate_only) res.XtWX = Eigen::MatrixXd::Zero(p, p);
    }

    if (!res.b.allFinite()) {
        res.b = Eigen::VectorXd::Constant(p, NA_REAL);
    }
    return res;
}

// [[Rcpp::export]]
List fast_ols_cpp(SEXP X_sexp, SEXP y_sexp,
                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values, true);
    return List::create(Named("b") = res.b);
}

// [[Rcpp::export]]
List fast_ols_with_var_cpp(SEXP X_sexp, SEXP y_sexp,
                           int j = 2,
                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_in(y_r.begin(), y_r.size());
    const int n = X.rows();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();

    Eigen::VectorXd y_to_use = y_in;
    for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) {
        y_to_use.noalias() -= X.col(fixed_spec.fixed_idx[k]) * fixed_spec.fixed_values[k];
    }
    const double yTy = y_to_use.squaredNorm();

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (int k = 0; k < fixed_spec.fixed_idx.size(); ++k) beta[fixed_spec.fixed_idx[k]] = fixed_spec.fixed_values[k];

    Eigen::VectorXd beta_free;
    Eigen::MatrixXd XtX_free(p_free, p_free);
    Eigen::VectorXd Xty_free(p_free);
    bool converged = false;

    if (p_free == p) {
        XtX_free = symmetric_crossprod(X);
        Xty_free.noalias() = X.transpose() * y_to_use;
        
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
        if (ldlt.info() == Eigen::Success) {
            beta_free = ldlt.solve(Xty_free);
            beta = beta_free;
            converged = true;
            
            double sse = std::max(0.0, yTy - beta_free.dot(Xty_free));
            double sigma2_hat = sse / (n - p);
            
            auto compute_ssq = [&](int col_idx) {
                if (col_idx < 0 || col_idx >= p) return NA_REAL;
                return sigma2_hat * ldlt.solve(Eigen::VectorXd::Unit(p, col_idx))(col_idx);
            };
            
            double ssq_j = compute_ssq(j - 1);
            double ssq_2 = (j == 2) ? ssq_j : compute_ssq(1);

            return List::create(
                Named("b") = beta,
                Named("XtX") = XtX_free,
                Named("ssq_b_j") = ssq_j,
                Named("ssq_b_2") = ssq_2,
                Named("sigma2_hat") = sigma2_hat,
                Named("converged") = true
            );
        } else {
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
            beta_free = qr.solve(y_to_use);
            beta = beta_free;
        }
    } else {
        Eigen::MatrixXd X_free(n, p_free);
        for (int k = 0; k < p_free; ++k) X_free.col(k) = X.col(fixed_spec.free_idx[k]);
        XtX_free = symmetric_crossprod(X_free);
        Xty_free.noalias() = X_free.transpose() * y_to_use;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX_free);
        if (ldlt.info() == Eigen::Success) {
            beta_free = ldlt.solve(Xty_free);
            for (int k = 0; k < p_free; ++k) beta[fixed_spec.free_idx[k]] = beta_free[k];
            converged = true;
            
            double sse = std::max(0.0, yTy - beta_free.dot(Xty_free));
            double sigma2_hat = sse / (n - p_free);
            
            auto compute_ssq = [&](int f_idx) {
                return sigma2_hat * ldlt.solve(Eigen::VectorXd::Unit(p_free, f_idx))(f_idx);
            };

            auto free_idx_of = [&](int k) -> int {
                for (int jj = 0; jj < p_free; ++jj) if (fixed_spec.free_idx[jj] == k) return jj;
                return -1;
            };

            int f_j = (j > 0 && j <= p) ? free_idx_of(j - 1) : -1;
            double ssq_j = (f_j >= 0) ? compute_ssq(f_j) : NA_REAL;
            int f_2 = (p >= 2) ? free_idx_of(1) : -1;
            double ssq_2 = (f_2 >= 0) ? ((f_2 == f_j) ? ssq_j : compute_ssq(f_2)) : NA_REAL;

            return List::create(
                Named("b") = beta,
                Named("XtX") = expand_free_covariance(p, fixed_spec, XtX_free, false),
                Named("ssq_b_j") = ssq_j,
                Named("ssq_b_2") = ssq_2,
                Named("sigma2_hat") = sigma2_hat,
                Named("converged") = true
            );
        } else {
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_free);
            beta_free = qr.solve(y_to_use);
            for (int k = 0; k < p_free; ++k) beta[fixed_spec.free_idx[k]] = beta_free[k];
        }
    }

    double sse = std::max(0.0, yTy - beta_free.dot(Xty_free));
    return List::create(
        Named("b") = beta,
        Named("converged") = converged,
        Named("sigma2_hat") = sse / (n - p_free)
    );
}
