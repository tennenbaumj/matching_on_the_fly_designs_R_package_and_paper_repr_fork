#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// --- Weight Functions ---

// Huber weight function
double huber_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) return 1.0;
    return c / abs_r;
}

// Tukey's Bisquare (Biweight) weight function
double bisquare_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) {
        double tmp = 1.0 - (r / c) * (r / c);
        return tmp * tmp;
    }
    return 0.0;
}

// --- Internal IRLS Logic ---

struct RobustModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd w;
    Eigen::MatrixXd XtWX;
    Eigen::MatrixXd X_free;
    double XtX_inv_diag_j;
    double scale;
    int iterations;
    bool converged;
    double ssq_b_j;

    RobustModelResult() : XtX_inv_diag_j(NA_REAL), scale(NA_REAL), iterations(0), converged(false), ssq_b_j(NA_REAL) {}
};

RobustModelResult fast_robust_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X, 
    const Eigen::Ref<const Eigen::VectorXd>& y, 
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    std::string method = "MM",
    double c = 1.345, // Huber constant
    double c_bisquare = 4.685, // Bisquare constant
    int maxit = 50,
    double tol = 1e-7,
    double scale_est = -1.0, // If negative, compute MAD
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
    bool estimate_only = false,
    int variance_j = 0
) {
    int n = X.rows();
    int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    Eigen::VectorXd y_adj = y;
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }
    RobustModelResult res;
    res.X_free = X_free;
    int free_variance_j = -1;
    if (variance_j > 0) {
        const int variance_j0 = variance_j - 1;
        for (int jj = 0; jj < fixed_spec.free_idx.size(); ++jj) {
            if (fixed_spec.free_idx[jj] == variance_j0) {
                free_variance_j = jj;
                break;
            }
        }
    }

    // 1. Initial estimate
    Eigen::VectorXd b_free;
    if (warm_start_beta.isNotNull()) {
        Eigen::VectorXd b_start = as<Eigen::VectorXd>(NumericVector(warm_start_beta));
        b_start = apply_fixed_values(b_start, fixed_spec);
        b_free = subset_vector(b_start, fixed_spec.free_idx);
    } else if (smart_cold_start) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_free);
        b_free = qr.solve(y_adj);
        if (!estimate_only && free_variance_j >= 0) {
            Eigen::MatrixXd R = Eigen::MatrixXd::Zero(p_free, p_free);
            R.triangularView<Eigen::Upper>() =
                qr.matrixR().topLeftCorner(p_free, p_free).template triangularView<Eigen::Upper>();
            Eigen::MatrixXd RtR = R.transpose() * R;
            Eigen::VectorXd e_orig = Eigen::VectorXd::Unit(p_free, free_variance_j);
            Eigen::VectorXd e_piv = qr.colsPermutation().transpose() * e_orig;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(RtR);
            if (ldlt.info() == Eigen::Success) {
                Eigen::VectorXd z = ldlt.solve(e_piv);
                if (z.allFinite()) {
                    res.XtX_inv_diag_j = e_piv.dot(z);
                }
            }
        }
    } else {
        b_free = Eigen::VectorXd::Zero(p_free);
    }
    res.b = Eigen::VectorXd::Zero(p);
    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = b_free[j];
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    Eigen::VectorXd r = y - X * res.b;
    
    // 2. Scale estimation (MAD of residuals)
    if (scale_est < 0) {
        std::vector<double> abs_r(n);
        for (int i = 0; i < n; ++i) abs_r[i] = std::abs(r[i]);
        double median_abs_r;
        if (n < 512) {
            // Full sorting is faster for small vectors because selection has higher fixed overhead.
            std::sort(abs_r.begin(), abs_r.end());
            median_abs_r = (n % 2 == 0) ? (abs_r[n / 2 - 1] + abs_r[n / 2]) / 2.0 : abs_r[n / 2];
        } else {
            const auto upper_mid = abs_r.begin() + n / 2;
            std::nth_element(abs_r.begin(), upper_mid, abs_r.end());
            median_abs_r = *upper_mid;
            if (n % 2 == 0) {
                const double lower_mid = *std::max_element(abs_r.begin(), upper_mid);
                median_abs_r = (lower_mid + median_abs_r) / 2.0;
            }
        }
        res.scale = median_abs_r / 0.6745;
    } else {
        res.scale = scale_est;
    }

    if (res.scale < 1e-10) res.scale = 1e-10;

    // 3. IRLS loop
    res.w = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b_old = b_free;
    
    for (int iter = 1; iter <= maxit; ++iter) {
        res.iterations = iter;
        
        // Update weights
        if (iter == 1 && warm_start_weights.isNotNull()) {
            Eigen::VectorXd ww = as<Eigen::VectorXd>(warm_start_weights);
            if (ww.size() != n) stop("warm_start_weights must have length equal to nrow(X)");
            res.w = ww;
        } else {
            const Eigen::ArrayXd u = r.array() / res.scale;
            if (method == "M") {
                const Eigen::ArrayXd abs_u = u.abs();
                res.w = (abs_u <= c).select(1.0, c / abs_u).matrix();
            } else {
                const Eigen::ArrayXd abs_u = u.abs();
                const Eigen::ArrayXd tmp = 1.0 - (u / c_bisquare).square();
                res.w = (abs_u <= c_bisquare).select(tmp.square(), 0.0).matrix();
            }
        }

        // Solve Weighted Least Squares
        Eigen::MatrixXd XtWX;
        if (iter == 1 && warm_start_fisher_info.isNotNull()) {
            Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
            if (info_full.rows() != p || info_full.cols() != p) stop("warm_start_fisher_info must be a p x p matrix");
            XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else {
            XtWX = weighted_crossprod(X_free, res.w);
        }
        Eigen::VectorXd XtWy = weighted_crossprod_rhs(X_free, res.w, y_adj);
        
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) break; // Numerical failure
        
        b_free = ldlt.solve(XtWy);
        for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = b_free[j];
        r = y - X * res.b;

        // Check convergence
        if ((b_free - b_old).norm() / (b_free.norm() + 1e-10) < tol) {
            res.converged = true;
            if (!estimate_only) res.XtWX = expand_free_covariance(p, fixed_spec, XtWX, false);
            break;
        }
        b_old = b_free;
    }

    return res;
}

//' @title Fast Robust Regression (C++)
//' @description High-performance robust regression fitting using IRLS.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param method Robust estimation method ("M" or "MM").
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param c Huber constant.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, weights, and scale estimate.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_robust_regression_cpp(
    SEXP X_sexp, 
    SEXP y_sexp, 
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    std::string method = "MM",
    int j = 2,
    double c = 1.345,
    int maxit = 50,
    double tol = 1e-7,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
    bool estimate_only = false
) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());

    RobustModelResult res = fast_robust_regression_internal(X, y, warm_start_beta, smart_cold_start, method, c, 4.685, maxit, tol, -1.0, fixed_idx, fixed_values, warm_start_weights, warm_start_fisher_info, estimate_only, j);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    
    if (estimate_only) {
        return List::create(
            Named("coefficients") = res.b,
            Named("scale") = res.scale,
            Named("converged") = res.converged,
            Named("iterations") = res.iterations
        );
    }

    if (!res.converged && res.XtWX.rows() == 0) {
        Eigen::MatrixXd XtWX_free = weighted_crossprod(res.X_free, res.w);
        res.XtWX = expand_free_covariance(X.cols(), fixed_spec, XtWX_free, false);
    }

    int n = X.rows();
    int p = X.cols();
    Eigen::VectorXd r = y - X * res.b;
    
    double ssq_j = NA_REAL;
    if (res.converged || res.iterations == maxit) {
        Eigen::VectorXd psi_r(n);
        double sum_psi_prime = 0;
        const Eigen::ArrayXd u = r.array() / res.scale;
        if (method == "M") {
            const Eigen::ArrayXd abs_u = u.abs();
            const Eigen::ArrayXd sign_u = (u > 0.0).select(
                Eigen::ArrayXd::Ones(n),
                -Eigen::ArrayXd::Ones(n)
            );
            psi_r = (abs_u <= c).select(r.array(), c * res.scale * sign_u).matrix();
            sum_psi_prime = (abs_u <= c).template cast<double>().sum();
        } else {
            const double c_b = 4.685;
            const Eigen::ArrayXd abs_u = u.abs();
            const Eigen::ArrayXd u_scaled_sq = (u / c_b).square();
            const Eigen::ArrayXd tmp = 1.0 - u_scaled_sq;
            psi_r = (abs_u <= c_b).select(r.array() * tmp.square(), 0.0).matrix();
            sum_psi_prime = (abs_u <= c_b).select(tmp * (1.0 - 5.0 * u_scaled_sq), 0.0).sum();
        }
        
        double m = sum_psi_prime / n;
        double sum_psi_sq = psi_r.squaredNorm();
        double factor = (n / (double(n - fixed_spec.free_idx.size()))) * sum_psi_sq / (n * m * m);

        if (j > 0 && j <= p) {
            const int j0 = j - 1;
            int free_j = -1;
            for (int jj = 0; jj < fixed_spec.free_idx.size(); ++jj) {
                if (fixed_spec.free_idx[jj] == j0) {
                    free_j = jj;
                    break;
                }
            }
            if (free_j >= 0) {
                double inv_diag = res.XtX_inv_diag_j;
                if (!R_finite(inv_diag)) {
                    Eigen::MatrixXd XtX = symmetric_crossprod(res.X_free);
                    inv_diag = compute_diagonal_inverse_entry(XtX, free_j + 1);
                }
                if (R_finite(inv_diag)) {
                    ssq_j = factor * inv_diag;
                }
            }
        }
    }

    return List::create(
        Named("coefficients") = res.b,
        Named("scale") = res.scale,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("ssq_b_j") = ssq_j,
        Named("fisher_information") = res.XtWX
    );
}
