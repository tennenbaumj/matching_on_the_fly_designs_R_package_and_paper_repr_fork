#include <RcppEigen.h>
#include <Rmath.h>
#include "_helper_functions.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

namespace {

// Helper for log-sum-exp trick to compute log(exp(a) + exp(b) - exp(c))
// Used for Clayton Copula logA: log(exp(theta*H1) + exp(theta*H2) - 1)
// where 1 = exp(0).
inline double log_sum_exp_clayton(double a, double b) {
    double m = a;
    if (b > m) m = b;
    if (0.0 > m) m = 0.0;
    double inner = std::exp(a - m) + std::exp(b - m) - std::exp(-m);
    if (inner <= 0) return m + std::log(std::numeric_limits<double>::min());
    return m + std::log(inner);
}

// -----------------------------------------------------------------------------
// Clayton Copula Weibull AFT
// -----------------------------------------------------------------------------

class ClaytonWeibullLikelihood {
private:
    Eigen::Ref<const Eigen::VectorXd> m_y;
    Eigen::Ref<const Eigen::VectorXd> m_dead;
    Eigen::Ref<const Eigen::MatrixXd> m_X;
    Eigen::Ref<const Eigen::MatrixXi> m_pair_idx;
    Eigen::Ref<const Eigen::VectorXi> m_singleton_rows;
    const bool m_has_pairs;
    const bool m_has_singletons;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;

public:
    ClaytonWeibullLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y, 
                             const Eigen::Ref<const Eigen::VectorXd>& dead, 
                             const Eigen::Ref<const Eigen::MatrixXd>& X,
                             const Eigen::Ref<const Eigen::MatrixXi>& pair_idx,
                             const Eigen::Ref<const Eigen::VectorXi>& singleton_rows) :
        m_y(y), m_dead(dead), m_X(X), m_pair_idx(pair_idx), 
        m_singleton_rows(singleton_rows),
        m_has_pairs(pair_idx.rows() > 0),
        m_has_singletons(singleton_rows.size() > 0),
        m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()) {}

    double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad) {
        double log_sigma = params[m_p];
        double log_theta = params[m_p + 1];
        double sigma = std::exp(log_sigma);
        double theta = std::exp(std::min(log_theta, 10.0));
        Eigen::VectorXd beta = params.head(m_p);

        Eigen::VectorXd eta = m_X * beta;
        Eigen::VectorXd log_H = (m_log_y - eta) / sigma;
        Eigen::VectorXd H(m_n);
        Eigen::VectorXd log_f(m_n);
        
        for (int i = 0; i < m_n; ++i) {
            double lH = log_H[i];
            if (lH > 700.0) lH = 700.0;
            H[i] = std::exp(lH);
            log_f[i] = lH - log_sigma - m_log_y[i] - H[i];
        }

        double loglik = 0.0;
        grad.setZero();
        
        Eigen::VectorXd d_loglik_d_eta = Eigen::VectorXd::Zero(m_n);
        double d_loglik_d_log_sigma = 0.0;
        double d_loglik_d_log_theta = 0.0;

        if (m_has_pairs) {
            for (int k = 0; k < m_pair_idx.rows(); ++k) {
                int i1 = m_pair_idx(k, 0);
                int i2 = m_pair_idx(k, 1);
                double h1 = H[i1];
                double h2 = H[i2];
                double d1 = m_dead[i1];
                double d2 = m_dead[i2];
                
                double logA = log_sum_exp_clayton(theta * h1, theta * h2);
                double A = std::exp(logA);
                
                // Common terms for derivatives
                double dA_d_h1 = theta * std::exp(theta * h1);
                double dA_d_h2 = theta * std::exp(theta * h2);
                double dA_d_theta = h1 * std::exp(theta * h1) + h2 * std::exp(theta * h2);

                if (d1 < 0.5 && d2 < 0.5) { // mask00
                    loglik -= (1.0 / theta) * logA;
                    // d/d_theta (-1/theta * logA) = 1/theta^2 * logA - 1/theta * 1/A * dA_d_theta
                    d_loglik_d_log_theta += (logA / (theta * theta) - dA_d_theta / (theta * A)) * theta;
                    
                    double d_ll_d_A = -1.0 / (theta * A);
                    d_loglik_d_eta[i1] += d_ll_d_A * dA_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += d_ll_d_A * dA_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += (d_ll_d_A * dA_d_h1 * (-H[i1] * log_H[i1]) + 
                                             d_ll_d_A * dA_d_h2 * (-H[i2] * log_H[i2]));
                } else if (d1 > 0.5 && d2 < 0.5) { // mask10
                    loglik += log_f[i1] + (-1.0/theta - 1.0) * logA + (theta + 1.0) * h1;
                    
                    d_loglik_d_log_theta += (logA / (theta * theta) + (-1.0/theta - 1.0) * dA_d_theta / A + h1) * theta;
                    
                    double d_ll_d_h1 = (theta + 1.0) + (-1.0/theta - 1.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (-1.0/theta - 1.0) * dA_d_h2 / A;
                    
                    // From log_f[i1]: d/d_eta = (H[i1] - 1)/sigma, d/d_log_sigma = -1 + (H[i1] - 1)*log_H[i1]
                    d_loglik_d_eta[i1] += (H[i1] - 1.0) / sigma + d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -1.0 + (H[i1] - 1.0) * log_H[i1] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                } else if (d1 < 0.5 && d2 > 0.5) { // mask01
                    loglik += log_f[i2] + (-1.0/theta - 1.0) * logA + (theta + 1.0) * h2;
                    
                    d_loglik_d_log_theta += (logA / (theta * theta) + (-1.0/theta - 1.0) * dA_d_theta / A + h2) * theta;
                    
                    double d_ll_d_h1 = (-1.0/theta - 1.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (theta + 1.0) + (-1.0/theta - 1.0) * dA_d_h2 / A;
                    
                    d_loglik_d_eta[i1] += d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += (H[i2] - 1.0) / sigma + d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -1.0 + (H[i2] - 1.0) * log_H[i2] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                } else { // mask11
                    loglik += std::log(theta + 1.0) + log_f[i1] + log_f[i2] + (-1.0/theta - 2.0) * logA + (theta + 1.0) * (h1 + h2);
                    
                    d_loglik_d_log_theta += (1.0/(theta + 1.0) + (1.0/(theta*theta)) * logA + (-1.0/theta - 2.0) * dA_d_theta / A + h1 + h2) * theta;
                    
                    double d_ll_d_h1 = (theta + 1.0) + (-1.0/theta - 2.0) * dA_d_h1 / A;
                    double d_ll_d_h2 = (theta + 1.0) + (-1.0/theta - 2.0) * dA_d_h2 / A;
                    
                    d_loglik_d_eta[i1] += (H[i1] - 1.0) / sigma + d_ll_d_h1 * (-H[i1] / sigma);
                    d_loglik_d_eta[i2] += (H[i2] - 1.0) / sigma + d_ll_d_h2 * (-H[i2] / sigma);
                    d_loglik_d_log_sigma += -2.0 + (H[i1] - 1.0) * log_H[i1] + (H[i2] - 1.0) * log_H[i2] + 
                                            d_ll_d_h1 * (-H[i1] * log_H[i1]) + 
                                            d_ll_d_h2 * (-H[i2] * log_H[i2]);
                }
            }
        }

        if (m_has_singletons) {
            for (int k = 0; k < m_singleton_rows.size(); ++k) {
                int i = m_singleton_rows[k];
                if (m_dead[i] > 0.5) {
                    loglik += log_f[i];
                    d_loglik_d_eta[i] += (H[i] - 1.0) / sigma;
                    d_loglik_d_log_sigma += -1.0 + (H[i] - 1.0) * log_H[i];
                } else {
                    loglik -= H[i];
                    d_loglik_d_eta[i] += H[i] / sigma;
                    d_loglik_d_log_sigma += H[i] * log_H[i];
                }
            }
        }

        grad.head(m_p) = - m_X.transpose() * d_loglik_d_eta;
        grad[m_p] = - d_loglik_d_log_sigma;
        grad[m_p + 1] = - d_loglik_d_log_theta;

        return -loglik;
    }
    
    Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params) {
        int total_p = params.size();
        Eigen::MatrixXd H(total_p, total_p);
        H.setZero();
        double h = 1e-6;
        Eigen::VectorXd grad_at_params(total_p);
        operator()(params, grad_at_params);

        for (int i = 0; i < total_p; ++i) {
            Eigen::VectorXd p_plus = params;
            p_plus[i] += h;
            Eigen::VectorXd g_plus(total_p);
            operator()(p_plus, g_plus);
            H.col(i) = (g_plus - grad_at_params) / h;
        }
        H = (H + H.transpose()) / 2.0;
        return H;
    }
};

// -----------------------------------------------------------------------------
// Dependent Censoring Transformation Regression
// -----------------------------------------------------------------------------

class DepCensTransformLikelihood {
private:
    Eigen::Ref<const Eigen::VectorXd> m_y;
    Eigen::Ref<const Eigen::VectorXd> m_dead;
    Eigen::Ref<const Eigen::MatrixXd> m_X;
    const int m_n;
    const int m_p;
    const Eigen::VectorXd m_log_y;
    mutable Eigen::VectorXd m_d_ll_d_mu_event;
    mutable Eigen::VectorXd m_d_ll_d_mu_cens;

public:
    DepCensTransformLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y,
                               const Eigen::Ref<const Eigen::VectorXd>& dead,
                               const Eigen::Ref<const Eigen::MatrixXd>& X) :
        m_y(y), m_dead(dead), m_X(X), m_n(y.size()), m_p(X.cols()),
        m_log_y(y.array().log().matrix()),
        m_d_ll_d_mu_event(y.size()), m_d_ll_d_mu_cens(y.size()) {}

    double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad) {
        // params: [beta_event (p), beta_cens (p), log_sigma_event (1), log_sigma_cens (1), atanh_rho (1)]
        Eigen::VectorXd beta_event = params.head(m_p);
        Eigen::VectorXd beta_cens = params.segment(m_p, m_p);
        double log_sigma_event = params[2 * m_p];
        double log_sigma_cens = params[2 * m_p + 1];
        double atanh_rho = params[2 * m_p + 2];
        atanh_rho = std::min(std::max(atanh_rho, -3.0), 3.0);

        double sigma_event = std::exp(log_sigma_event);
        double sigma_cens = std::exp(log_sigma_cens);
        double rho = std::tanh(atanh_rho);
        double one_minus_rho_sq = std::max(1.0 - rho * rho, 1e-12);
        double sd_cond = std::sqrt(one_minus_rho_sq);

        Eigen::VectorXd mu_event = m_X * beta_event;
        Eigen::VectorXd mu_cens = m_X * beta_cens;
        Eigen::VectorXd z_event = (m_log_y - mu_event) / sigma_event;
        Eigen::VectorXd z_cens = (m_log_y - mu_cens) / sigma_cens;

        double loglik = 0.0;
        grad.setZero();

        m_d_ll_d_mu_event.setZero();
        m_d_ll_d_mu_cens.setZero();
        double d_ll_d_log_sigma_event = 0.0;
        double d_ll_d_log_sigma_cens = 0.0;
        double d_ll_d_atanh_rho = 0.0;

        const double omrs_sqrt = std::sqrt(one_minus_rho_sq);
        const double omrs_1p5  = one_minus_rho_sq * omrs_sqrt;

        for (int i = 0; i < m_n; ++i) {
            double ze = z_event[i];
            double zc = z_cens[i];
            double d = m_dead[i];

            double log_f_event = fast_log_dnorm(ze) - log_sigma_event - m_log_y[i];
            double log_f_cens  = fast_log_dnorm(zc) - log_sigma_cens  - m_log_y[i];

            double w_event = (rho * ze - zc) / sd_cond;
            double w_cens  = (rho * zc - ze) / sd_cond;

            double log_S_cond_cens  = fast_log_pnorm(w_event);
            double log_S_cond_event = fast_log_pnorm(w_cens);

            double li = d * (log_f_event + log_S_cond_cens) + (1.0 - d) * (log_f_cens + log_S_cond_event);
            loglik += li;

            // Derivatives — Mills ratio via inline log-pdf minus log-survival
            double mill_event = std::exp(fast_log_dnorm(w_event) - log_S_cond_cens);
            double mill_cens  = std::exp(fast_log_dnorm(w_cens)  - log_S_cond_event);

            // d/d_ze w_event = rho / sd_cond
            // d/d_zc w_event = -1 / sd_cond
            // d/d_rho w_event = (ze * sd_cond - (rho*ze - zc) * (-rho/sd_cond)) / (1-rho^2)
            //                = (ze * (1-rho^2) + rho*(rho*ze - zc)) / (1-rho^2)^(3/2)
            //                = (ze - rho^2*ze + rho^2*ze - rho*zc) / (1-rho^2)^(3/2)
            //                = (ze - rho*zc) / (1-rho^2)^(3/2)
            
            double d_w_event_d_rho = (ze - rho * zc) / omrs_1p5;
            double d_w_cens_d_rho  = (zc - rho * ze) / omrs_1p5;

            if (d > 0.5) {
                // d * log_f_event
                m_d_ll_d_mu_event[i] += ze / sigma_event;
                d_ll_d_log_sigma_event += (ze * ze - 1.0);

                // d * log_S_cond_cens
                m_d_ll_d_mu_event[i] += mill_event * (rho / sd_cond) * (-1.0 / sigma_event);
                m_d_ll_d_mu_cens[i]  += mill_event * (-1.0 / sd_cond) * (-1.0 / sigma_cens);
                d_ll_d_log_sigma_event += mill_event * (rho / sd_cond) * (-ze);
                d_ll_d_log_sigma_cens  += mill_event * (-1.0 / sd_cond) * (-zc);
                d_ll_d_atanh_rho += mill_event * d_w_event_d_rho * one_minus_rho_sq;
            } else {
                // (1-d) * log_f_cens
                m_d_ll_d_mu_cens[i] += zc / sigma_cens;
                d_ll_d_log_sigma_cens += (zc * zc - 1.0);

                // (1-d) * log_S_cond_event
                m_d_ll_d_mu_cens[i]  += mill_cens * (rho / sd_cond) * (-1.0 / sigma_cens);
                m_d_ll_d_mu_event[i] += mill_cens * (-1.0 / sd_cond) * (-1.0 / sigma_event);
                d_ll_d_log_sigma_cens  += mill_cens * (rho / sd_cond) * (-zc);
                d_ll_d_log_sigma_event += mill_cens * (-1.0 / sd_cond) * (-ze);
                d_ll_d_atanh_rho += mill_cens * d_w_cens_d_rho * one_minus_rho_sq;
            }
        }

        grad.head(m_p) = - m_X.transpose() * m_d_ll_d_mu_event;
        grad.segment(m_p, m_p) = - m_X.transpose() * m_d_ll_d_mu_cens;
        grad[2 * m_p] = - d_ll_d_log_sigma_event;
        grad[2 * m_p + 1] = - d_ll_d_log_sigma_cens;
        grad[2 * m_p + 2] = - d_ll_d_atanh_rho;

        return -loglik;
    }

    Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params) {
        {
            int total_p = params.size();
            Eigen::MatrixXd H(total_p, total_p);
            H.setZero();
            double h = 1e-6;
            Eigen::VectorXd grad_at_params(total_p);
            operator()(params, grad_at_params);

            for (int i = 0; i < total_p; ++i) {
                Eigen::VectorXd p_plus = params;
                p_plus[i] += h;
                Eigen::VectorXd g_plus(total_p);
                operator()(p_plus, g_plus);
                H.col(i) = (g_plus - grad_at_params) / h;
            }
            H = (H + H.transpose()) / 2.0;
            return H;
        }
        const int total_p = params.size();
        const Eigen::VectorXd beta_event = params.head(m_p);
        const Eigen::VectorXd beta_cens = params.segment(m_p, m_p);
        const double log_sigma_event = params[2 * m_p];
        const double log_sigma_cens = params[2 * m_p + 1];
        double atanh_rho = params[2 * m_p + 2];
        atanh_rho = std::min(std::max(atanh_rho, -3.0), 3.0);

        const double sigma_e = std::exp(log_sigma_event);
        const double sigma_c = std::exp(log_sigma_cens);
        const double rho = std::tanh(atanh_rho);
        const double omr2 = std::max(1.0 - rho * rho, 1e-12);
        const double s_cond = std::sqrt(omr2);

        const Eigen::VectorXd mu_e = m_X * beta_event;
        const Eigen::VectorXd mu_c = m_X * beta_cens;
        const Eigen::VectorXd z_e = (m_log_y - mu_e) / sigma_e;
        const Eigen::VectorXd z_c = (m_log_y - mu_c) / sigma_c;

        Eigen::MatrixXd Hess = Eigen::MatrixXd::Zero(total_p, total_p);

        const double omr2_sqrt = std::sqrt(omr2);
        const double omr2_1p5  = omr2 * omr2_sqrt;
        const double omr2_2p5  = omr2 * omr2 * omr2_sqrt;

        for (int i = 0; i < m_n; ++i) {
            const double ze = z_e[i], zc = z_c[i], d = m_dead[i];
            const double w_e = (rho * ze - zc) / s_cond;
            const double w_c = (rho * zc - ze) / s_cond;
            const double log_S_cond_c = fast_log_pnorm(w_e);
            const double log_S_cond_e = fast_log_pnorm(w_c);
            const double mill_e = std::exp(fast_log_dnorm(w_e) - log_S_cond_c);
            const double mill_c = std::exp(fast_log_dnorm(w_c) - log_S_cond_e);

            const double dm_e = -mill_e * (w_e + mill_e);
            const double dm_c = -mill_c * (w_c + mill_c);

            const double d_we_d_ze = rho / s_cond, d_we_d_zc = -1.0 / s_cond;
            const double d_wc_d_zc = rho / s_cond, d_wc_d_ze = -1.0 / s_cond;
            const double d_we_d_rho = (ze - rho * zc) / omr2_1p5;
            const double d_wc_d_rho = (zc - rho * ze) / omr2_1p5;
            const double d_rho_d_ath = omr2;
            const double d2_rho_d_ath2 = -2.0 * rho * omr2;

            const Eigen::VectorXd xi = m_X.row(i).transpose();
            auto add_h = [&](int r, int c, double v) { Hess(r, c) += v; };

            if (d > 0.5) {
                add_h(2*m_p, 2*m_p, -2.0 * ze * ze); 
                for(int j=0; j<m_p; j++) {
                    add_h(j, j, -1.0/(sigma_e*sigma_e) * xi[j] * xi[j]);
                    add_h(j, 2*m_p, -2.0*ze/sigma_e * xi[j]);
                }
                
                const double h_zeze = dm_e * d_we_d_ze * d_we_d_ze;
                const double h_zczc = dm_e * d_we_d_zc * d_we_d_zc;
                const double h_zezc = dm_e * d_we_d_ze * d_we_d_zc;
                const double d2_we_d_ze_d_rho = 1.0 / omr2_1p5;
                const double d2_we_d_zc_d_rho = -rho / omr2_1p5;
                const double d2_we_d_rho2 = (3.0 * rho * ze - (1.0 + 2.0*rho*rho)*zc) / omr2_2p5;

                const double dze_de = -1.0/sigma_e, dze_dls = -ze;
                const double dzc_dc = -1.0/sigma_c, dzc_dlsc = -zc;
                
                for(int j=0; j<m_p; j++) {
                    add_h(j, j, h_zeze * dze_de * dze_de * xi[j] * xi[j]);
                    add_h(m_p+j, m_p+j, h_zczc * dzc_dc * dzc_dc * xi[j] * xi[j]);
                    add_h(j, m_p+j, h_zezc * dze_de * dzc_dc * xi[j] * xi[j]);
                    add_h(j, 2*m_p, (h_zeze*dze_dls + mill_e*d_we_d_ze)*dze_de * xi[j]);
                    add_h(j, 2*m_p+1, h_zezc*dzc_dlsc*dze_de * xi[j]);
                    add_h(m_p+j, 2*m_p, h_zezc*dze_dls*dzc_dc * xi[j]);
                    add_h(m_p+j, 2*m_p+1, (h_zczc*dzc_dlsc + mill_e*d_we_d_zc)*dzc_dc * xi[j]);
                    add_h(j, total_p-1, (dm_e*d_we_d_ze*d_we_d_rho + mill_e*d2_we_d_ze_d_rho)*d_rho_d_ath*dze_de * xi[j]);
                    add_h(m_p+j, total_p-1, (dm_e*d_we_d_zc*d_we_d_rho + mill_e*d2_we_d_zc_d_rho)*d_rho_d_ath*dzc_dc * xi[j]);
                }
                add_h(2*m_p, 2*m_p, h_zeze*dze_dls*dze_dls + mill_e*d_we_d_ze*ze);
                add_h(2*m_p+1, 2*m_p+1, h_zczc*dzc_dlsc*dzc_dlsc + mill_e*d_we_d_zc*zc);
                add_h(2*m_p, 2*m_p+1, h_zezc*dze_dls*dzc_dlsc);
                add_h(2*m_p, total_p-1, (dm_e*d_we_d_ze*d_we_d_rho + mill_e*d2_we_d_ze_d_rho)*d_rho_d_ath*dze_dls);
                add_h(2*m_p+1, total_p-1, (dm_e*d_we_d_zc*d_we_d_rho + mill_e*d2_we_d_zc_d_rho)*d_rho_d_ath*dzc_dlsc);
                add_h(total_p-1, total_p-1, (dm_e*d_we_d_rho*d_we_d_rho + mill_e*d2_we_d_rho2)*d_rho_d_ath*d_rho_d_ath + mill_e*d_we_d_rho*d2_rho_d_ath2);

            } else {
                add_h(2*m_p+1, 2*m_p+1, -2.0 * zc * zc);
                for(int j=0; j<m_p; j++) {
                    add_h(m_p+j, m_p+j, -1.0/(sigma_c*sigma_c) * xi[j] * xi[j]);
                    add_h(m_p+j, 2*m_p+1, -2.0*zc/sigma_c * xi[j]);
                }

                const double d2_wc_d_zc_d_rho = 1.0 / omr2_1p5;
                const double d2_wc_d_ze_d_rho = -rho / omr2_1p5;
                const double d2_wc_d_rho2 = (3.0 * rho * zc - (1.0 + 2.0*rho*rho)*ze) / omr2_2p5;

                const double dze_de = -1.0/sigma_e, dze_dls = -ze;
                const double dzc_dc = -1.0/sigma_c, dzc_dlsc = -zc;
                
                for(int j=0; j<m_p; j++) {
                    add_h(m_p+j, m_p+j, dm_c*d_wc_d_zc*d_wc_d_zc * dzc_dc * dzc_dc * xi[j] * xi[j]);
                    add_h(j, j, dm_c*d_wc_d_ze*d_wc_d_ze * dze_de * dze_de * xi[j] * xi[j]);
                    add_h(j, m_p+j, dm_c*d_wc_d_zc*d_wc_d_ze * dzc_dc * dze_de * xi[j] * xi[j]);
                    add_h(m_p+j, 2*m_p+1, (dm_c*d_wc_d_zc*d_wc_d_zc*dzc_dlsc + mill_c*d_wc_d_zc)*dzc_dc * xi[j]);
                    add_h(m_p+j, 2*m_p, dm_c*d_wc_d_zc*d_wc_d_ze*dze_dls*dzc_dc * xi[j]);
                    add_h(j, 2*m_p+1, dm_c*d_wc_d_zc*d_wc_d_ze*dzc_dlsc*dze_de * xi[j]);
                    add_h(j, 2*m_p, (dm_c*d_wc_d_ze*d_wc_d_ze*dze_dls + mill_c*d_wc_d_ze)*dze_de * xi[j]);
                    add_h(m_p+j, total_p-1, (dm_c*d_wc_d_zc*d_wc_d_rho + mill_c*d2_wc_d_zc_d_rho)*d_rho_d_ath*dzc_dc * xi[j]);
                    add_h(j, total_p-1, (dm_c*d_wc_d_ze*d_wc_d_rho + mill_c*d2_wc_d_ze_d_rho)*d_rho_d_ath*dze_de * xi[j]);
                }
                add_h(2*m_p+1, 2*m_p+1, dm_c*d_wc_d_zc*d_wc_d_zc*dzc_dlsc*dzc_dlsc + mill_c*d_wc_d_zc*zc);
                add_h(2*m_p, 2*m_p, dm_c*d_wc_d_ze*d_wc_d_ze*dze_dls*dze_dls + mill_c*d_wc_d_ze*ze);
                add_h(2*m_p, 2*m_p+1, dm_c*d_wc_d_zc*d_wc_d_ze*dze_dls*dzc_dlsc);
                add_h(2*m_p+1, total_p-1, (dm_c*d_wc_d_zc*d_wc_d_rho + mill_c*d2_wc_d_zc_d_rho)*d_rho_d_ath*dzc_dlsc);
                add_h(2*m_p, total_p-1, (dm_c*d_wc_d_ze*d_wc_d_rho + mill_c*d2_wc_d_ze_d_rho)*d_rho_d_ath*dze_dls);
                add_h(total_p-1, total_p-1, (dm_c*d_wc_d_rho*d_wc_d_rho + mill_c*d2_wc_d_rho2)*d_rho_d_ath*d_rho_d_ath + mill_c*d_wc_d_rho*d2_rho_d_ath2);
            }
        }
        for(int r1=0; r1<total_p; r1++) for(int c1=0; r1>c1; c1++) Hess(r1, c1) = Hess(c1, r1);
        return -Hess;
    }
};

} // namespace

// -----------------------------------------------------------------------------
// R-exposed functions
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP get_clayton_weibull_aft_score_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    SEXP pair_idx_sexp,
    SEXP singleton_rows_sexp,
    SEXP params_sexp
) {
    NumericMatrix X_mat(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector dead_vec(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());
    IntegerMatrix pair_idx_mat(pair_idx_sexp);
    Eigen::Map<const Eigen::MatrixXi> pair_idx(pair_idx_mat.begin(), pair_idx_mat.nrow(), pair_idx_mat.ncol());
    IntegerVector singleton_rows_vec(singleton_rows_sexp);
    Eigen::Map<const Eigen::VectorXi> singleton_rows(singleton_rows_vec.begin(), singleton_rows_vec.size());
    NumericVector params_vec(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());

    ClaytonWeibullLikelihood fun(y, dead, X, pair_idx, singleton_rows);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return wrap(-grad);
}

// [[Rcpp::export]]
SEXP get_clayton_weibull_aft_hessian_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    SEXP pair_idx_sexp,
    SEXP singleton_rows_sexp,
    SEXP params_sexp
) {
    NumericMatrix X_mat(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector dead_vec(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());
    IntegerMatrix pair_idx_mat(pair_idx_sexp);
    Eigen::Map<const Eigen::MatrixXi> pair_idx(pair_idx_mat.begin(), pair_idx_mat.nrow(), pair_idx_mat.ncol());
    IntegerVector singleton_rows_vec(singleton_rows_sexp);
    Eigen::Map<const Eigen::VectorXi> singleton_rows(singleton_rows_vec.begin(), singleton_rows_vec.size());
    NumericVector params_vec(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());

    ClaytonWeibullLikelihood fun(y, dead, X, pair_idx, singleton_rows);
    return wrap(-fun.hessian(params));
}

// [[Rcpp::export]]
SEXP get_dep_cens_transform_score_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    SEXP params_sexp
) {
    NumericMatrix X_mat(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector dead_vec(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());
    NumericVector params_vec(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());

    DepCensTransformLikelihood fun(y, dead, X);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return wrap(-grad);
}

// [[Rcpp::export]]
SEXP get_dep_cens_transform_hessian_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    SEXP params_sexp
) {
    NumericMatrix X_mat(X_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    NumericVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector dead_vec(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());
    NumericVector params_vec(params_sexp);
    Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());

    DepCensTransformLikelihood fun(y, dead, X);
    return wrap(-fun.hessian(params));
}

// [[Rcpp::export]]
List fast_clayton_weibull_aft_optim_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    SEXP pair_idx_sexp,
    SEXP singleton_rows_sexp,
    SEXP warm_start_params_sexp,
    bool estimate_only = false,
    int maxit = 2000,
    double reltol = 1e-9,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    NumericMatrix X_mat(X_sexp);
    NumericVector y_vec(y_sexp);
    NumericVector dead_vec(dead_sexp);
    IntegerMatrix pair_idx_int_mat(pair_idx_sexp);
    IntegerVector singleton_rows_int_vec(singleton_rows_sexp);
    NumericVector warm_start_params_vec(warm_start_params_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());
    Eigen::Map<const Eigen::MatrixXi> pair_idx(pair_idx_int_mat.begin(), pair_idx_int_mat.nrow(), pair_idx_int_mat.ncol());
    Eigen::Map<const Eigen::VectorXi> singleton_rows(singleton_rows_int_vec.begin(), singleton_rows_int_vec.size());
    Eigen::Map<const Eigen::VectorXd> warm_start_params(warm_start_params_vec.begin(), warm_start_params_vec.size());

    ClaytonWeibullLikelihood fun(y, dead, X, pair_idx, singleton_rows);
    Eigen::VectorXd params = warm_start_params;
    FixedParamSpec fixed_spec = make_fixed_param_spec(params.size(), fixed_idx, fixed_values);
    params = apply_fixed_values(params, fixed_spec);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, reltol, optimization_alg, "lbfgs", 0, h_ptr);
    } catch (const std::exception& e) {
        Rcout << "Optimization failed: " << e.what() << std::endl;
        const Eigen::VectorXd p_fixed = apply_fixed_values(params, fixed_spec);
        double last_val = NA_REAL;
        try {
            Eigen::VectorXd g(p_fixed.size());
            last_val = fun(p_fixed, g);
        } catch (...) {}
        return List::create(
            Named("converged") = false, Named("error") = e.what(),
            Named("par") = p_fixed, Named("params") = p_fixed, Named("b") = p_fixed,
            Named("value") = last_val, Named("neg_loglik") = last_val, Named("neg_ll") = last_val
        );
    } catch (...) {
        const Eigen::VectorXd p_fixed = apply_fixed_values(params, fixed_spec);
        double last_val = NA_REAL;
        try {
            Eigen::VectorXd g(p_fixed.size());
            last_val = fun(p_fixed, g);
        } catch (...) {}
        return List::create(
            Named("converged") = false, Named("error") = "unknown",
            Named("par") = p_fixed, Named("params") = p_fixed, Named("b") = p_fixed,
            Named("value") = last_val, Named("neg_loglik") = last_val, Named("neg_ll") = last_val
        );
    }
    params = fit.params;
    
    List out = List::create(
        Named("par") = params,
        Named("params") = params,
        Named("b") = params,
        Named("value") = fit.value,
        Named("neg_loglik") = fit.value,
        Named("neg_ll") = fit.value,
        Named("loglik") = R_finite(fit.value) ? -fit.value : NA_REAL,
        Named("niter") = fit.niter,
        Named("converged") = fit.converged
    );

    if (estimate_only) {
        return out;
    }

    Eigen::MatrixXd observed_information = fun.hessian(params);
    Eigen::VectorXd score(params.size());
    fun(params, score);
    score = -score;
    Eigen::MatrixXd vcov = covariance_from_information(observed_information);
    out["score"] = score;
    out["observed_information"] = observed_information;
    out["information"] = observed_information;
    out["information_type"] = "observed";
    out["hessian"] = -observed_information;
    out["fisher_information"] = observed_information;
    out["vcov"] = vcov;

    return out;
}

// [[Rcpp::export]]
List fast_dep_cens_transform_optim_cpp(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP dead_sexp,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
    bool smart_cold_start = true,
    bool estimate_only = false,
    int maxit = 2000,
    double reltol = 1e-9,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    NumericMatrix X_mat(X_sexp);
    NumericVector y_vec(y_sexp);
    NumericVector dead_vec(dead_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXd> dead(dead_vec.begin(), dead_vec.size());

    int p = X.cols();
    int total = 2 * p + 3;
    Eigen::VectorXd params(total);
    if (warm_start_params.isNotNull()) {
        params = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_params));
    } else if (smart_cold_start) {
        params.setZero();
        Eigen::VectorXd log_y = (y.array() + 1e-8).log().matrix();
        Eigen::VectorXd b_ols = safe_ols_solve(X, log_y);
        params.head(p) = b_ols; // beta_event
        params.segment(p, p) = b_ols; // beta_cens
        params[2 * p] = 0.0; // log_sigma_event
        params[2 * p + 1] = 0.0; // log_sigma_cens
        params[2 * p + 2] = 0.0; // log_theta
    } else {
        params.setZero();
    }
    
    DepCensTransformLikelihood fun(y, dead, X);
    FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);
    params = apply_fixed_values(params, fixed_spec);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit;
    try {
        fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, reltol, optimization_alg, "lbfgs", 0, h_ptr);
    } catch (const std::exception& e) {
        Rcout << "Optimization failed: " << e.what() << std::endl;
        return List::create(Named("converged") = false, Named("error") = e.what());
    } catch (...) {
        return List::create(Named("converged") = false, Named("error") = "unknown");
    }
    params = fit.params;
    
    List out = List::create(
        Named("par") = params,
        Named("params") = params,
        Named("b") = params,
        Named("value") = fit.value,
        Named("neg_loglik") = fit.value,
        Named("neg_ll") = fit.value,
        Named("loglik") = R_finite(fit.value) ? -fit.value : NA_REAL,
        Named("niter") = fit.niter,
        Named("converged") = fit.converged
    );

    if (estimate_only) {
        return out;
    }

    Eigen::MatrixXd observed_information = fun.hessian(params);
    Eigen::VectorXd score(params.size());
    fun(params, score);
    score = -score;
    Eigen::MatrixXd vcov = covariance_from_information(observed_information);
    out["score"] = score;
    out["observed_information"] = observed_information;
    out["information"] = observed_information;
    out["information_type"] = "observed";
    out["hessian"] = -observed_information;
    out["fisher_information"] = observed_information;
    out["vcov"] = vcov;

    return out;
}
