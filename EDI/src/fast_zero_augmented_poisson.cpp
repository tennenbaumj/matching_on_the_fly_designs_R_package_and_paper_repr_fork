#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

// Helper for log1pexp(x) = log(1 + exp(x))
double log1pexp(double x) {
    if (x > 0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

class ZeroAugmentedPoisson {
private:
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    const Eigen::Ref<const Eigen::MatrixXd> m_Xzi;
    const int m_n;
    const int m_p_cond;
    const int m_p_zi;
    const bool m_is_hurdle;
    // Preallocated scratch vectors — avoid heap allocation on every operator() call.
    // operator() is non-const so no mutable needed; optimizer holds functor by non-const ref.
    Eigen::VectorXd m_eta_cond, m_eta_zi, m_w_cond, m_w_zi;

public:
    ZeroAugmentedPoisson(const Eigen::Ref<const Eigen::VectorXd>& y,
                         const Eigen::Ref<const Eigen::MatrixXd>& X,
                         const Eigen::Ref<const Eigen::MatrixXd>& Xzi,
                         bool is_hurdle) :
        m_y(y), m_X(X), m_Xzi(Xzi), m_n(X.rows()),
        m_p_cond(X.cols()), m_p_zi(Xzi.cols()), m_is_hurdle(is_hurdle),
        m_eta_cond(m_n), m_eta_zi(m_n), m_w_cond(m_n), m_w_zi(m_n) {}

    double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
        m_eta_cond.noalias() = m_X   * params.head(m_p_cond);
        m_eta_zi.noalias()   = m_Xzi * params.tail(m_p_zi);

        double neg_ll = 0.0;
        m_w_cond.setZero();
        m_w_zi.setZero();

        for (int i = 0; i < m_n; ++i) {
            const double eta_c = std::min(m_eta_cond[i], 700.0);
            const double eta_z = std::max(std::min(m_eta_zi[i], 700.0), -700.0);
            const double lam   = std::exp(eta_c);
            // Compute exp(-eta_z) once; reuse for sigmoid AND log1pexp to save one exp() per obs.
            const double en    = std::exp(-eta_z);
            const double p     = 1.0 / (1.0 + en);
            // log1pexp(eta_z) via the stable identity that reuses en:
            //   eta_z > 0: eta_z + log1p(exp(-eta_z))  [exp(-eta_z) < 1]
            //   eta_z <= 0: log1p(1/exp(-eta_z))        [1/exp(-eta_z) < 1]
            const double lse   = (eta_z > 0.0)
                ? eta_z + std::log1p(en)
                : std::log1p(1.0 / en);

            if (m_is_hurdle) {
                if (m_y[i] == 0) {
                    neg_ll -= std::log(std::max(p, 1e-15));
                    m_w_zi[i] = -(1.0 - p);
                } else {
                    const double eml = std::exp(-lam);
                    const double log1m_eml = (lam > 0.6931471805599453)
                        ? std::log1p(-eml)
                        : std::log(-std::expm1(-lam));
                    neg_ll -= (-lse + m_y[i] * eta_c - lam - log1m_eml);
                    m_w_zi[i]   = p;
                    m_w_cond[i] = -(m_y[i] - lam / (1.0 - eml));
                }
            } else {
                if (m_y[i] == 0) {
                    const double eml       = std::exp(-lam);
                    const double prob_zero = p + (1.0 - p) * eml;
                    neg_ll -= std::log(std::max(prob_zero, 1e-15));
                    const double inv_pz = 1.0 / prob_zero;
                    m_w_zi[i]   = -(p * (1.0 - p) * (1.0 - eml)) * inv_pz;
                    m_w_cond[i] =  (1.0 - p) * lam * eml * inv_pz;
                } else {
                    neg_ll -= (-lse + m_y[i] * eta_c - lam);
                    m_w_zi[i]   = p;
                    m_w_cond[i] = -(m_y[i] - lam);
                }
            }
        }

        grad.resize(m_p_cond + m_p_zi);
        grad.head(m_p_cond).noalias() = m_X.transpose()   * m_w_cond;
        grad.tail(m_p_zi).noalias()   = m_Xzi.transpose() * m_w_zi;

        return neg_ll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
        int total_p = m_p_cond + m_p_zi;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta_cond = params.head(m_p_cond);
        Eigen::VectorXd beta_zi = params.tail(m_p_zi);
        Eigen::VectorXd eta_cond = m_X * beta_cond;
        Eigen::VectorXd eta_zi = m_Xzi * beta_zi;

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double eta_c = std::min(eta_cond[i], 700.0);
            double lambda = std::exp(eta_c);
            double eta_z = std::max(std::min(eta_zi[i], 700.0), -700.0);
            double pi = 1.0 / (1.0 + std::exp(-eta_z));
            double pi_prime = pi * (1.0 - pi);
            const double* xci = m_X.data() + i;    // xci[j * m_n] == X(i,j)
            const double* xzi = m_Xzi.data() + i;  // xzi[j * m_n] == Xzi(i,j)

            if (m_is_hurdle) {
                // zi-zi block (upper triangle): rows/cols m_p_cond..total_p-1
                for (int c = 0; c < m_p_zi; ++c) {
                    const double s = pi_prime * xzi[c * m_n];
                    for (int r = 0; r <= c; ++r)
                        H_data[(m_p_cond + r) + (m_p_cond + c) * total_p] += s * xzi[r * m_n];
                }
                if (m_y[i] > 0) {
                    double exp_ml = std::exp(-lambda);
                    double denom = std::max(1.0 - exp_ml, 1e-15);
                    double h_cond = lambda * (denom - lambda * exp_ml) / (denom * denom);
                    // cond-cond block (upper triangle): rows/cols 0..m_p_cond-1
                    for (int c = 0; c < m_p_cond; ++c) {
                        const double s = h_cond * xci[c * m_n];
                        for (int r = 0; r <= c; ++r)
                            H_data[r + c * total_p] += s * xci[r * m_n];
                    }
                }
            } else {
                if (m_y[i] == 0) {
                    double r = std::exp(-lambda);
                    double q = std::max(pi + (1.0 - pi) * r, 1e-15);
                    double pi_second = pi_prime * (1.0 - 2.0 * pi);
                    double q_z = pi_prime * (1.0 - r);
                    double q_zz = pi_second * (1.0 - r);
                    double q_c = -(1.0 - pi) * lambda * r;
                    double q_cc = (1.0 - pi) * r * (lambda * lambda - lambda);
                    double q_zc = pi_prime * lambda * r;

                    double h_cc = q_c * q_c / (q * q) - q_cc / q;
                    double h_zz = q_z * q_z / (q * q) - q_zz / q;
                    double h_cz = q_c * q_z / (q * q) - q_zc / q;

                    // cond-cond block (upper triangle)
                    for (int c = 0; c < m_p_cond; ++c) {
                        const double s = h_cc * xci[c * m_n];
                        for (int r = 0; r <= c; ++r)
                            H_data[r + c * total_p] += s * xci[r * m_n];
                    }
                    // zi-zi block (upper triangle)
                    for (int c = 0; c < m_p_zi; ++c) {
                        const double s = h_zz * xzi[c * m_n];
                        for (int r = 0; r <= c; ++r)
                            H_data[(m_p_cond + r) + (m_p_cond + c) * total_p] += s * xzi[r * m_n];
                    }
                    // cond-zi cross block (full, not symmetric)
                    for (int c = 0; c < m_p_zi; ++c)
                        for (int r = 0; r < m_p_cond; ++r)
                            H_data[r + (m_p_cond + c) * total_p] += h_cz * xci[r * m_n] * xzi[c * m_n];
                } else {
                    // cond-cond block (upper triangle)
                    for (int c = 0; c < m_p_cond; ++c) {
                        const double s = lambda * xci[c * m_n];
                        for (int r = 0; r <= c; ++r)
                            H_data[r + c * total_p] += s * xci[r * m_n];
                    }
                    // zi-zi block (upper triangle)
                    for (int c = 0; c < m_p_zi; ++c) {
                        const double s = pi_prime * xzi[c * m_n];
                        for (int r = 0; r <= c; ++r)
                            H_data[(m_p_cond + r) + (m_p_cond + c) * total_p] += s * xzi[r * m_n];
                    }
                }
            }
        }
        // Reflect upper triangle to lower (handles cond-cond, zi-zi, and cond-zi/zi-cond blocks)
        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }

    Eigen::MatrixXd expected_hessian(const Eigen::VectorXd& params) {
        int total_p = m_p_cond + m_p_zi;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
        Eigen::VectorXd beta_cond = params.head(m_p_cond);
        Eigen::VectorXd beta_zi = params.tail(m_p_zi);
        Eigen::VectorXd eta_cond = m_X * beta_cond;
        Eigen::VectorXd eta_zi = m_Xzi * beta_zi;

        double* H_data = H.data();
        for (int i = 0; i < m_n; ++i) {
            double eta_c = std::min(eta_cond[i], 700.0);
            double lambda = std::exp(eta_c);
            double eta_z = std::max(std::min(eta_zi[i], 700.0), -700.0);
            double pi = 1.0 / (1.0 + std::exp(-eta_z));
            double pi_prime = pi * (1.0 - pi);
            const double* xci = m_X.data() + i;
            const double* xzi = m_Xzi.data() + i;

            if (m_is_hurdle) {
                for (int c = 0; c < m_p_zi; ++c) {
                    const double s = pi_prime * xzi[c * m_n];
                    for (int r = 0; r <= c; ++r)
                        H_data[(m_p_cond + r) + (m_p_cond + c) * total_p] += s * xzi[r * m_n];
                }

                double exp_ml = std::exp(-lambda);
                double denom = std::max(1.0 - exp_ml, 1e-15);
                double h_cond = lambda * (denom - lambda * exp_ml) / (denom * denom);
                double p_positive = 1.0 - pi;
                for (int c = 0; c < m_p_cond; ++c) {
                    const double s = p_positive * h_cond * xci[c * m_n];
                    for (int r = 0; r <= c; ++r)
                        H_data[r + c * total_p] += s * xci[r * m_n];
                }
            } else {
                double r0 = std::exp(-lambda);
                double q0 = std::max(pi + (1.0 - pi) * r0, 1e-15);
                double p_positive = (1.0 - pi) * (1.0 - r0);
                double pi_second = pi_prime * (1.0 - 2.0 * pi);
                double q_z = pi_prime * (1.0 - r0);
                double q_zz = pi_second * (1.0 - r0);
                double q_c = -(1.0 - pi) * lambda * r0;
                double q_cc = (1.0 - pi) * r0 * (lambda * lambda - lambda);
                double q_zc = pi_prime * lambda * r0;

                double h0_cc = q_c * q_c / (q0 * q0) - q_cc / q0;
                double h0_zz = q_z * q_z / (q0 * q0) - q_zz / q0;
                double h0_cz = q_c * q_z / (q0 * q0) - q_zc / q0;

                double h_cc = q0 * h0_cc + p_positive * lambda;
                double h_zz = q0 * h0_zz + p_positive * pi_prime;
                double h_cz = q0 * h0_cz;

                for (int c = 0; c < m_p_cond; ++c) {
                    const double s = h_cc * xci[c * m_n];
                    for (int rr = 0; rr <= c; ++rr)
                        H_data[rr + c * total_p] += s * xci[rr * m_n];
                }
                for (int c = 0; c < m_p_zi; ++c) {
                    const double s = h_zz * xzi[c * m_n];
                    for (int rr = 0; rr <= c; ++rr)
                        H_data[(m_p_cond + rr) + (m_p_cond + c) * total_p] += s * xzi[rr * m_n];
                }
                for (int c = 0; c < m_p_zi; ++c)
                    for (int rr = 0; rr < m_p_cond; ++rr)
                        H_data[rr + (m_p_cond + c) * total_p] += h_cz * xci[rr * m_n] * xzi[c * m_n];
            }
        }
        for (int c = 0; c < total_p; ++c)
            for (int r = 0; r < c; ++r)
                H_data[c + r * total_p] = H_data[r + c * total_p];
        return H;
    }
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_zero_augmented_poisson_score_cpp(SEXP X_sexp,
													 SEXP y_sexp,
													 SEXP Xzi_sexp,
													 SEXP params_sexp,
													 bool is_hurdle) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericMatrix Xzi_r(Xzi_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::MatrixXd> Xzi(Xzi_r.begin(), Xzi_r.nrow(), Xzi_r.ncol());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
	Eigen::VectorXd grad(params.size());
	fun(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_zero_augmented_poisson_hessian_cpp(SEXP X_sexp,
													   SEXP y_sexp,
													   SEXP Xzi_sexp,
													   SEXP params_sexp,
													   bool is_hurdle) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericMatrix Xzi_r(Xzi_sexp);
    NumericVector params_r(params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::MatrixXd> Xzi(Xzi_r.begin(), Xzi_r.nrow(), Xzi_r.ncol());
    Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
	return -fun.hessian(params);
}

//' @title Fast Zero-Augmented Poisson Regression (C++)
//' @description High-performance ZIP or hurdle Poisson regression fitting using Newton-Raphson or L-BFGS.
//' @param X Matrix of predictors for the conditional component.
//' @param y Vector of responses.
//' @param Xzi Matrix of predictors for the zero-inflation/hurdle component.
//' @param is_hurdle If TRUE, fit a hurdle model; if FALSE, fit a zero-inflated model.
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param estimate_only If TRUE, skip variance component calculations.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first iteration.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_zero_augmented_poisson_cpp(SEXP X_sexp,
                                     SEXP y_sexp,
                                     SEXP Xzi_sexp,
                                     bool is_hurdle,
                                     Nullable<NumericVector> warm_start_params = R_NilValue,
                                     bool smart_cold_start = true,
                                     bool estimate_only = false,
                                     int maxit = 1000,
                                     double tol = 1e-8,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "lbfgs",
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_r(X_sexp);
    NumericVector y_r(y_sexp);
    NumericMatrix Xzi_r(Xzi_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::MatrixXd> Xzi(Xzi_r.begin(), Xzi_r.nrow(), Xzi_r.ncol());

    int p_cond = X.cols();
    int p_zi = Xzi.cols();
    int total_p = p_cond + p_zi;

    Eigen::VectorXd params(total_p);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(NumericVector(warm_start_params));
    } else if (smart_cold_start) {
        params = edi_opt::zap_smart_cold_start(X, Xzi, y);
    } else {
        params.setZero();
        // Naive warm_start_params for intercept
        double mean_y = y.mean();
        if (mean_y > 0) params[0] = std::log(mean_y);
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(total_p, fixed_idx, fixed_values);

    ZeroAugmentedPoisson fun(y, X, Xzi, is_hurdle);
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }

    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, h_ptr);
    } catch (...) {
        return List::create(Named("converged") = false);
    }
    params = fit.params;

    if (estimate_only) {
        return List::create(
            Named("params") = params,
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.value,
            Named("neg_loglik") = fit.value
        );
    }

    Eigen::MatrixXd observed_information = fun.hessian(params);

    Eigen::MatrixXd H_free = subset_matrix(observed_information, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd vcov = expand_free_covariance(total_p, fixed_spec, cov_free, true);

    return List::create(
        Named("coefficients") = List::create(
            Named("cond") = params.head(p_cond),
            Named("zi") = params.tail(p_zi)
        ),
        Named("params") = params,
        Named("vcov") = vcov,
        Named("converged") = fit.converged,
        Named("neg_ll") = fit.value,
        Named("neg_loglik") = fit.value,
        Named("observed_information") = observed_information,
        Named("fisher_information") = observed_information,
        Named("information") = observed_information,
        Named("information_type") = "observed",
        Named("hessian") = -observed_information
    );
}
