#ifndef EDI_GLMM_ENGINE_H
#define EDI_GLMM_ENGINE_H

#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include "_helper_functions.h"

namespace glmm {

struct GHRule {
    Eigen::VectorXd nodes;
    Eigen::VectorXd log_norm_weights;
};

inline GHRule gauss_hermite_rule(int n) {
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n - 1; ++i) {
        const double v = std::sqrt((i + 1.0) / 2.0);
        J(i, i + 1) = v;
        J(i + 1, i) = v;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
    GHRule rule;
    rule.nodes = es.eigenvalues();
    rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
                            - 0.5 * std::log(M_PI);
    return rule;
}

inline double log_sum_exp(const Eigen::VectorXd& x) {
    const double m = x.maxCoeff();
    if (!std::isfinite(m)) return m;
    return m + std::log((x.array() - m).exp().sum());
}

// Fused log_sum_exp + posterior weights: computes lse(x) and fills
// weights[k] = exp(x[k] - lse) in a single pass — avoids recomputing
// exp(x[k] - lse) separately in the gradient accumulation loop.
inline double log_sum_exp_and_weights(const Eigen::VectorXd& x, Eigen::VectorXd& weights) {
    const double m = x.maxCoeff();
    if (!std::isfinite(m)) { weights.setZero(); return m; }
    weights = (x.array() - m).exp();
    const double S = weights.sum();
    weights /= S;
    return m + std::log(S);
}

struct GLMMData {
    Eigen::MatrixXd X_s;
    Eigen::VectorXd y_s;
    std::vector<int> grp_start, grp_size;
    int n, p, G;
    GHRule gh;
    double max_abs_log_sigma;

    GLMMData(const Eigen::MatrixXd& X, 
             const Eigen::VectorXd& y, 
             const std::vector<int>& gid,
             int n_gh, double max_abs_log_sigma_)
        : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule(n_gh)), max_abs_log_sigma(max_abs_log_sigma_) {
        
        std::vector<int> ord(n);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){ return gid[a] < gid[b]; });

        X_s.resize(n, p);
        y_s.resize(n);
        for(int i = 0; i < n; ++i) {
            X_s.row(i) = X.row(ord[i]);
            y_s[i] = y[ord[i]];
        }

        int prev = -1;
        for(int i = 0; i < n; ++i) {
            if(gid[ord[i]] != prev) {
                grp_start.push_back(i);
                grp_size.push_back(1);
                prev = gid[ord[i]];
            } else {
                grp_size.back()++;
            }
        }
        G = static_cast<int>(grp_start.size());
    }
};

inline double sigma_penalty(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return scale * d * d;
}

inline double sigma_penalty_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double sigma_penalty_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
    const double d = std::abs(log_sigma) - center;
    if (d <= 0.0) return 0.0;
    return 2.0 * scale;
}

// A generic GLMM Objective that uses Gauss-Hermite quadrature.
// Model template must provide:
// - int n_model_params() const
// - void fill_alpha(const VectorXd& par, double* buf) const
//     Precompute any per-step model parameters into caller-allocated buf (length n_model_params()).
//     Called once per optimizer step before the inner loops.
// - double log_prob(double y, double eta, const double* buf) const
// - double log_prob_derivs(double y, double eta, const double* buf, double& d_eta, VectorXd& d_par) const
template <typename Model>
class GLMMObjective {
    const GLMMData& dat;
    Model model;
    // Preallocated once at construction — written by fill_alpha before each loop,
    // only read inside the GH quadrature loop (no LICM concern; inner loop never writes it).
    mutable std::vector<double> m_alpha_buf;

public:
    GLMMObjective(const GLMMData& d, const Model& m)
        : dat(d), model(m), m_alpha_buf(m.n_model_params()) {}

    double value(const Eigen::VectorXd& par) const {
        const int nm = model.n_model_params();
        const double log_sigma = par[nm + dat.p];
        if (!std::isfinite(log_sigma) || std::abs(log_sigma) > dat.max_abs_log_sigma) return 1e100;
        const double sigma = std::exp(log_sigma);
        const double pen = sigma_penalty(log_sigma);

        model.fill_alpha(par, m_alpha_buf.data());

        const Eigen::VectorXd beta = par.segment(nm, dat.p);
        const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
        const int nn = static_cast<int>(b_vals.size());

        double total_ll = 0.0;
        for (int gi = 0; gi < dat.G; ++gi) {
            const int start = dat.grp_start[gi];
            const int sz = dat.grp_size[gi];
            const Eigen::VectorXd eta0 = dat.X_s.middleRows(start, sz) * beta;

            Eigen::VectorXd log_terms(nn);
            for (int k = 0; k < nn; ++k) {
                double ll = dat.gh.log_norm_weights[k];
                for (int r = 0; r < sz; ++r) {
                    ll += model.log_prob(dat.y_s[start + r], eta0[r] + b_vals[k], m_alpha_buf.data());
                }
                log_terms[k] = ll;
            }
            const double ll_g = log_sum_exp(log_terms);
            if (!std::isfinite(ll_g)) return 1e100;
            total_ll += ll_g;
        }
        return -total_ll + pen;
    }

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        const int nm = model.n_model_params();
        const double log_sigma = par[nm + dat.p];
        const double sigma = std::exp(log_sigma);
        const double pen = sigma_penalty(log_sigma);

        model.fill_alpha(par, m_alpha_buf.data());

        const Eigen::VectorXd beta = par.segment(nm, dat.p);
        const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
        const int nn = static_cast<int>(b_vals.size());

        double total_nll = pen;
        grad.setZero();
        grad[nm + dat.p] += sigma_penalty_grad(log_sigma);

        for (int gi = 0; gi < dat.G; ++gi) {
            const int start = dat.grp_start[gi];
            const int sz = dat.grp_size[gi];
            const auto Xg = dat.X_s.middleRows(start, sz);
            const Eigen::VectorXd eta0 = Xg * beta;

            Eigen::VectorXd log_terms(nn);
            Eigen::MatrixXd dL_dp_sum_nodes = Eigen::MatrixXd::Zero(nm, nn);
            Eigen::MatrixXd dL_de_nodes(sz, nn);
            Eigen::VectorXd dp(nm);

            for (int k = 0; k < nn; ++k) {
                double ll = dat.gh.log_norm_weights[k];
                for (int r = 0; r < sz; ++r) {
                    double de;
                    double lp = model.log_prob_derivs(dat.y_s[start + r], eta0[r] + b_vals[k], m_alpha_buf.data(), de, dp);
                    ll += lp;
                    dL_dp_sum_nodes.col(k).noalias() += dp;
                    dL_de_nodes(r, k) = de;
                }
                log_terms[k] = ll;
            }

            Eigen::VectorXd pk_vec(nn);
            const double ll_g = log_sum_exp_and_weights(log_terms, pk_vec);
            total_nll -= ll_g;

            for (int k = 0; k < nn; ++k) {
                const double pk = pk_vec[k];
                if (pk < 1e-15) continue;

                const Eigen::VectorXd dLi_db = Xg.transpose() * dL_de_nodes.col(k);
                grad.head(nm) -= pk * dL_dp_sum_nodes.col(k);
                grad.segment(nm, dat.p) -= pk * dLi_db;
                grad[nm + dat.p] -= pk * dL_de_nodes.col(k).sum() * std::sqrt(2.0) * dat.gh.nodes[k] * sigma;
            }
        }
        return total_nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        return numerical_hessian(*this, par); // Fallback to numerical Hessian for simplicity in generic engine
    }
};

} // namespace glmm

#endif
