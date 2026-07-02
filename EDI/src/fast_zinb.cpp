// Zero-Inflated Negative Binomial (ZINB) regression.
//
// Model:
//   P(Y=0)   = pi + (1-pi) * (theta/(theta+mu))^theta
//   P(Y=y>0) = (1-pi) * NegBin(y; mu, theta)
//
// Parameter vector: [beta_cond(p_cond), beta_zi(p_zi), log_theta]
//   mu  = exp(eta_cond),  pi = sigmoid(eta_zi),  theta = exp(log_theta)
//
// Analytic gradient; analytic Hessian.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;

namespace {

// log(1 + exp(x)) — numerically stable
inline double lse_zinb(double x) {
    if (x > 0.0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

// sigmoid(x)
inline double sigmoid_zinb(double x) {
    if (x >  35.0) return 1.0;
    if (x < -35.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

class ZeroInflatedNegBin {
    const Eigen::Ref<const Eigen::VectorXd> m_y;
    const Eigen::Ref<const Eigen::MatrixXd> m_Xc;
    const Eigen::Ref<const Eigen::MatrixXd> m_Xz;
    const int m_n, m_pc, m_pz;

    // Precomputed at construction from fixed data y.
    // m_y_slot[i]: index into m_distinct_y for obs i (or -1 if y_i == 0).
    // m_distinct_y: sorted unique positive y values (as double).
    // m_lgamma_y1: lgamma(y+1) for each entry of m_distinct_y — constant across calls.
    std::vector<int>    m_y_slot;
    std::vector<double> m_distinct_y;
    std::vector<double> m_lgamma_y1;
    std::vector<double> m_lgamma_yptheta;   // preallocated; filled per operator() call
    std::vector<double> m_digamma_yptheta;  // preallocated; filled per operator() call

    // Preallocated scratch vectors — avoid heap allocation on every operator() call.
    // operator() is non-const, so mutating members across calls is safe (matches ZeroAugmentedPoisson).
    Eigen::VectorXd m_eta_c, m_eta_z, m_w_c, m_w_z;

public:
    ZeroInflatedNegBin(const Eigen::Ref<const Eigen::VectorXd>& y,
                       const Eigen::Ref<const Eigen::MatrixXd>& Xc,
                       const Eigen::Ref<const Eigen::MatrixXd>& Xz)
        : m_y(y), m_Xc(Xc), m_Xz(Xz),
          m_n(y.size()), m_pc(Xc.cols()), m_pz(Xz.cols()),
          m_y_slot(y.size(), -1),
          m_eta_c(m_n), m_eta_z(m_n), m_w_c(m_n), m_w_z(m_n)
    {
        std::unordered_map<int, int> seen;
        for (int i = 0; i < m_n; ++i) {
            if (m_y[i] <= 0.0) continue;
            const int yi_int = static_cast<int>(m_y[i]);
            auto it = seen.find(yi_int);
            if (it == seen.end()) {
                const int slot = static_cast<int>(m_distinct_y.size());
                seen[yi_int] = slot;
                m_distinct_y.push_back(static_cast<double>(yi_int));
                m_lgamma_y1.push_back(R::lgammafn(static_cast<double>(yi_int) + 1.0));
                m_y_slot[i] = slot;
            } else {
                m_y_slot[i] = it->second;
            }
        }
        const int nd = static_cast<int>(m_distinct_y.size());
        m_lgamma_yptheta.resize(nd);
        m_digamma_yptheta.resize(nd);
    }

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        const double log_theta = par[m_pc + m_pz];
        const double theta     = std::exp(std::min(log_theta, 700.0));

        m_eta_c.noalias() = m_Xc * par.head(m_pc);
        m_eta_z.noalias() = m_Xz * par.segment(m_pc, m_pz);

        // Hoist theta-only special functions out of the observation loop.
        const double digamma_theta = fast_digamma(theta);
        const double lgamma_theta  = R::lgammafn(theta);

        // Fill preallocated per-distinct-y tables.
        const int nd = static_cast<int>(m_distinct_y.size());
        for (int k = 0; k < nd; ++k) {
            const double ypt = m_distinct_y[k] + theta;
            m_lgamma_yptheta[k]  = R::lgammafn(ypt);
            m_digamma_yptheta[k] = fast_digamma(ypt);
        }

        double nll = 0.0;
        double d_log_theta = 0.0;
        m_w_c.setZero();
        m_w_z.setZero();

        for (int i = 0; i < m_n; ++i) {
            const double yi  = m_y[i];
            const double ec  = m_eta_c[i];
            const double ez  = m_eta_z[i];
            const double p   = sigmoid_zinb(ez);
            const double mu  = std::exp(ec);
            const double denom     = theta + mu;
            const double log_denom = std::log(denom);

            if (yi <= 0.0) {
                // phi = (theta/(theta+mu))^theta = exp(theta*(log_theta - log_denom))
                const double phi = std::exp(theta * (log_theta - log_denom));
                const double den = p + (1.0 - p) * phi;
                nll -= std::log(den);

                const double inv_den = 1.0 / den;
                m_w_z[i] = -(p * (1.0 - p) * (1.0 - phi)) * inv_den;

                const double d_phi_d_ec = -mu * theta * phi / denom;
                m_w_c[i] = -((1.0 - p) * d_phi_d_ec) * inv_den;

                const double d_phi_d_theta = phi * (log_theta - log_denom + 1.0 - theta / denom);
                d_log_theta -= (1.0 - p) * d_phi_d_theta * theta * inv_den;
            } else {
                const int slot = m_y_slot[i];
                nll -= -lse_zinb(ez)  // log(1-p)
                     + m_lgamma_yptheta[slot] - lgamma_theta - m_lgamma_y1[slot]
                     + theta * (log_theta - log_denom) + yi * (ec - log_denom);

                m_w_z[i] = p;
                m_w_c[i] = -(yi - mu * (yi + theta) / denom);

                d_log_theta -= (m_digamma_yptheta[slot] - digamma_theta + log_theta - log_denom + 1.0 - (yi + theta) / denom) * theta;
            }
        }

        grad.resize(m_pc + m_pz + 1);
        grad.head(m_pc).noalias()          = m_Xc.transpose() * m_w_c;
        grad.segment(m_pc, m_pz).noalias() = m_Xz.transpose() * m_w_z;
        grad[m_pc + m_pz]                  = d_log_theta;
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        return numerical_hessian(*this, par);
    }
};

} // namespace

//' @title Fast Zero-Inflated Negative Binomial Regression (C++)
//' @description High-performance zero-inflated negative binomial model fitting via L-BFGS.
//' @param X Numeric matrix of predictors for the count component (including intercept).
//' @param Xzi Numeric matrix of predictors for the zero-inflation component (including intercept).
//' @param y Numeric vector of non-negative integer count responses.
//' @param warm_start_params Optional starting values for all parameters.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm (default "lbfgs").
//' @param smart_cold_start Logical. If TRUE, use a heuristic initial guess.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix.
//' @param estimate_only Logical. If TRUE, skip variance computation and return only coefficients.
//' @return A list containing coefficients and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_zinb_cpp(SEXP X, SEXP Xzi, SEXP y,
                   Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                   int maxit = 1000, double tol = 1e-8,
                   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                   std::string optimization_alg = "lbfgs",
                   bool smart_cold_start = true,
                   Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                   bool estimate_only = false) {
    NumericMatrix Xc_r(X);
    NumericMatrix Xz_r(Xzi);
    NumericVector y_r(y);
    Eigen::Map<const Eigen::MatrixXd> Xc(Xc_r.begin(), Xc_r.nrow(), Xc_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> Xz(Xz_r.begin(), Xz_r.nrow(), Xz_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_vec(y_r.begin(), y_r.size());

    ZeroInflatedNegBin obj(y_vec, Xc, Xz);
    int n_par = Xc.cols() + Xz.cols() + 1;
    Eigen::VectorXd par = Eigen::VectorXd::Zero(n_par);
    
    if (warm_start_params.isNotNull()) {
        par = as<Eigen::VectorXd>(NumericVector(warm_start_params));
        if (par.size() != n_par) stop("warm_start_params must have length equal to the number of model parameters");
    } else if (smart_cold_start) {
        // Simple initialization
        double mean_y = y_vec.mean();
        if (mean_y < 1e-8) mean_y = 1e-8;
        par.head(Xc.cols()).setZero();
        if (Xc.cols() > 0) par[0] = std::log(mean_y);
        par.segment(Xc.cols(), Xz.cols()).setZero();
        par[n_par - 1] = std::log(1.0); // theta = 1
    }

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);
    
    Eigen::MatrixXd info;
    const Eigen::MatrixXd* info_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info = as<Eigen::MatrixXd>(warm_start_fisher_info);
        if (info.rows() == n_par && info.cols() == n_par) {
            info_ptr = &info;
        }
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_ptr);
    
    const int p_cond = Xc.cols();
    const int p_zi = Xz.cols();
    if (estimate_only) {
        return List::create(
            Named("params") = fit.params,
            Named("coefficients") = List::create(
                Named("cond") = fit.params.head(p_cond),
                Named("zi") = fit.params.segment(p_cond, p_zi)
            ),
            Named("converged") = fit.converged,
            Named("iterations") = fit.niter
        );
    }

    Eigen::MatrixXd hess = obj.hessian(fit.params);
    Rcpp::List out = make_uniform_likelihood_fit_result(fit.params, fit.value, fit.converged, -likelihood_score(obj, fit.params), hess, false);
    out["coefficients"] = List::create(
        Named("cond") = fit.params.head(p_cond),
        Named("zi") = fit.params.segment(p_cond, p_zi)
    );
    return out;
}
