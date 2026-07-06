#include "_helper_functions.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class StereotypeLogitRegression {
private:
    const Eigen::Ref<const Eigen::MatrixXd> m_X;
    std::vector<int> m_y;
    int m_n;
    int m_p;
    int m_K;

public:
    StereotypeLogitRegression(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXd>& y) :
        m_X(X), m_n(X.rows()), m_p(X.cols()), m_K(0) {
        std::vector<double> levels = init_levels(y);
        m_K = static_cast<int>(levels.size());
        m_y.resize(m_n);
        for (int i = 0; i < m_n; ++i) {
            double yi = y[i];
            int idx = 0;
            while (idx < m_K && levels[idx] != yi) {
                ++idx;
            }
            m_y[i] = idx + 1;
        }
    }

    static std::vector<double> init_levels(const Eigen::Ref<const Eigen::VectorXd>& y) {
        std::vector<double> levels(y.data(), y.data() + y.size());
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
        return levels;
    }

    int num_categories() const { return m_K; }
    int num_alpha() const { return m_K - 1; }
    int num_gamma() const { return std::max(0, m_K - 2); }
    int num_params() const { return num_alpha() + m_p + num_gamma(); }

    VectorXd initialize_params() const {
        VectorXd params(num_params());
        params.setZero();

        std::vector<double> counts(m_K, 0.5);
        for (int i = 0; i < m_n; ++i) {
            counts[m_y[i] - 1] += 1.0;
        }
        for (int j = 1; j < m_K; ++j) {
            params[j - 1] = std::log(counts[j] / counts[0]);
        }
        return params;
    }

    void compute_scores(
        const Eigen::Ref<const VectorXd>& gamma,
        std::vector<double>& score_vals,
        Eigen::Ref<MatrixXd> dscore_dgamma
    ) const {
        score_vals.assign(m_K, 0.0);
        dscore_dgamma.setZero();

        if (num_gamma() == 0) {
            if (m_K >= 2) {
                score_vals[m_K - 1] = 1.0;
            }
            return;
        }

        VectorXd v = gamma.array().exp();
        double denom = 1.0 + v.sum();
        VectorXd cum_v(num_gamma());
        double running = 0.0;
        for (int r = 0; r < num_gamma(); ++r) {
            running += v[r];
            cum_v[r] = running;
        }

        for (int j = 2; j <= m_K - 1; ++j) {
            int interior_idx = j - 2;
            double c_j = cum_v[interior_idx];
            score_vals[j - 1] = c_j / denom;

            for (int r = 0; r < num_gamma(); ++r) {
                double dc = (r <= interior_idx) ? v[r] : 0.0;
                dscore_dgamma(j - 1, r) = (dc * denom - c_j * v[r]) / (denom * denom);
            }
        }
        score_vals[m_K - 1] = 1.0;
    }

    void compute_scores_with_second_derivatives(
        const Eigen::Ref<const VectorXd>& gamma,
        std::vector<double>& score_vals,
        Eigen::Ref<MatrixXd> dscore_dgamma,
        std::vector<MatrixXd>& d2score_dgamma2
    ) const {
        const int n_gamma = num_gamma();
        compute_scores(gamma, score_vals, dscore_dgamma);
        d2score_dgamma2.assign(m_K, MatrixXd::Zero(n_gamma, n_gamma));

        if (n_gamma == 0) {
            return;
        }

        VectorXd v = gamma.array().exp();
        const double denom = 1.0 + v.sum();
        VectorXd cum_v(n_gamma);
        double running = 0.0;
        for (int r = 0; r < n_gamma; ++r) {
            running += v[r];
            cum_v[r] = running;
        }

        for (int j = 2; j <= m_K - 1; ++j) {
            const int interior_idx = j - 2;
            const double c_j = cum_v[interior_idx];
            MatrixXd& Hs = d2score_dgamma2[j - 1];

            for (int r = 0; r < n_gamma; ++r) {
                const double A_r = (r <= interior_idx) ? 1.0 : 0.0;
                const double numerator_r = A_r * denom - c_j;
                const double first_r = v[r] * numerator_r / (denom * denom);

                for (int t = 0; t < n_gamma; ++t) {
                    const double A_t = (t <= interior_idx) ? 1.0 : 0.0;
                    const double delta_rt = (r == t) ? 1.0 : 0.0;
                    Hs(r, t) =
                        delta_rt * first_r +
                        v[r] * v[t] * (
                            (A_r - A_t) / (denom * denom) -
                            2.0 * numerator_r / (denom * denom * denom)
                        );
                }
            }
        }
    }

    double loglik_grad(
        const Eigen::Ref<const VectorXd>& params,
        VectorXd* grad = NULL
    ) const {
        const int n_alpha = num_alpha();
        const int n_gamma = num_gamma();

        VectorXd alpha = params.head(n_alpha);
        VectorXd beta = params.segment(n_alpha, m_p);
        VectorXd gamma = (n_gamma > 0) ? params.tail(n_gamma) : VectorXd(0);

        if (grad != NULL) {
            grad->setZero();
        }

        std::vector<double> score_vals;
        MatrixXd dscore_dgamma(m_K, n_gamma);
        compute_scores(gamma, score_vals, dscore_dgamma);
        const Map<const VectorXd> score_vec(score_vals.data(), m_K);

        VectorXd logits = VectorXd::Zero(m_K);
        VectorXd probs = VectorXd::Zero(m_K);
        double ll = 0.0;

        for (int i = 0; i < m_n; ++i) {
            const double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;
            logits.setZero();
            for (int j = 2; j <= m_K; ++j) {
                logits[j - 1] = alpha[j - 2] + score_vec[j - 1] * eta;
            }

            const double max_logit = logits.maxCoeff();
            probs = (logits.array() - max_logit).exp().matrix();
            const double denom = probs.sum();
            probs /= denom;

            const int yi = m_y[i] - 1;
            ll += logits[yi] - max_logit - std::log(denom);

            if (grad != NULL) {
                for (int j = 2; j <= m_K; ++j) {
                    (*grad)[j - 2] += ((yi == (j - 1)) ? 1.0 : 0.0) - probs[j - 1];
                }

                if (m_p > 0) {
                    const double observed_score = score_vec[yi];
                    const double expected_score = probs.dot(score_vec);
                    grad->segment(n_alpha, m_p).noalias() += m_X.row(i).transpose() * (observed_score - expected_score);
                }

                if (n_gamma > 0) {
                    const VectorXd expected_dscore = dscore_dgamma.transpose() * probs;
                    grad->tail(n_gamma).noalias() += eta * (dscore_dgamma.row(yi).transpose() - expected_dscore);
                }
            }
        }

        return ll;
    }

    MatrixXd loglik_hessian(const Eigen::Ref<const VectorXd>& params) const {
        const int n_alpha = num_alpha();
        const int n_gamma = num_gamma();
        const int d = num_params();

        VectorXd beta = params.segment(n_alpha, m_p);
        VectorXd gamma = (n_gamma > 0) ? params.tail(n_gamma) : VectorXd(0);

        std::vector<double> score_vals;
        MatrixXd dscore_dgamma(m_K, n_gamma);
        std::vector<MatrixXd> d2score_dgamma2;
        compute_scores_with_second_derivatives(gamma, score_vals, dscore_dgamma, d2score_dgamma2);
        const Map<const VectorXd> score_vec(score_vals.data(), m_K);

        MatrixXd H = MatrixXd::Zero(d, d);
        VectorXd logits = VectorXd::Zero(m_K);
        VectorXd probs = VectorXd::Zero(m_K);
        std::vector<VectorXd> logit_grad(m_K, VectorXd::Zero(d));
        std::vector<MatrixXd> logit_hess(m_K, MatrixXd::Zero(d, d));

        for (int i = 0; i < m_n; ++i) {
            const double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;
            logits.setZero();
            logit_grad[0].setZero();
            logit_hess[0].setZero();

            for (int j = 2; j <= m_K; ++j) {
                const int cat = j - 1;
                logits[cat] = params[j - 2] + score_vec[cat] * eta;

                VectorXd& zj = logit_grad[cat];
                MatrixXd& Bj = logit_hess[cat];
                zj.setZero();
                Bj.setZero();

                zj[j - 2] = 1.0;
                if (m_p > 0) {
                    zj.segment(n_alpha, m_p).noalias() = score_vec[cat] * m_X.row(i).transpose();
                }
                if (n_gamma > 0) {
                    zj.tail(n_gamma).noalias() = eta * dscore_dgamma.row(cat).transpose();
                    for (int r = 0; r < n_gamma; ++r) {
                        const double ds = dscore_dgamma(cat, r);
                        for (int b = 0; b < m_p; ++b) {
                            const int beta_idx = n_alpha + b;
                            const int gamma_idx = n_alpha + m_p + r;
                            const double cross = m_X(i, b) * ds;
                            Bj(beta_idx, gamma_idx) = cross;
                            Bj(gamma_idx, beta_idx) = cross;
                        }
                    }
                    Bj.bottomRightCorner(n_gamma, n_gamma).noalias() =
                        eta * d2score_dgamma2[cat];
                }
            }

            const double max_logit = logits.maxCoeff();
            probs = (logits.array() - max_logit).exp().matrix();
            probs /= probs.sum();

            VectorXd mean_grad = VectorXd::Zero(d);
            MatrixXd mean_hess = MatrixXd::Zero(d, d);
            MatrixXd mean_outer = MatrixXd::Zero(d, d);
            double* mo = mean_outer.data();
            for (int j = 0; j < m_K; ++j) {
                mean_grad.noalias() += probs[j] * logit_grad[j];
                mean_hess.noalias() += probs[j] * logit_hess[j];
                // Raw pointer accumulation for mean_outer (upper triangle only)
                const double* gj = logit_grad[j].data();
                const double pj = probs[j];
                for (int c = 0; c < d; ++c) {
                    const double s = pj * gj[c];
                    for (int r = 0; r <= c; ++r)
                        mo[r + c * d] += s * gj[r];
                }
            }
            // Reflect mean_outer upper triangle to lower
            for (int c = 0; c < d; ++c)
                for (int r = 0; r < c; ++r)
                    mo[c + r * d] = mo[r + c * d];

            const int yi = m_y[i] - 1;
            // Build delta = logit_hess[yi] - mean_hess - mean_outer, then add mean_grad outer product
            Eigen::MatrixXd delta = logit_hess[yi] - mean_hess - mean_outer;
            double* d_data = delta.data();
            const double* mg = mean_grad.data();
            // Add mean_grad * mean_grad.T (upper triangle only)
            for (int c = 0; c < d; ++c) {
                const double mg_c = mg[c];
                for (int r = 0; r <= c; ++r)
                    d_data[r + c * d] += mg_c * mg[r];
            }
            // Reflect delta upper triangle to lower
            for (int c = 0; c < d; ++c)
                for (int r = 0; r < c; ++r)
                    d_data[c + r * d] = d_data[r + c * d];
            H.noalias() += delta;
        }

        return 0.5 * (H + H.transpose());
    }

    // Fisher information: E[-d2LL/dθ2] = Σ_i [Σ_j p_j z_j z_j^T - mean_z mean_z^T]
    // No d2score/dgamma2 needed — cheap and always PSD.
    MatrixXd expected_hessian(const Eigen::Ref<const VectorXd>& params) const {
        const int n_alpha = num_alpha();
        const int n_gamma = num_gamma();
        const int d = num_params();

        const VectorXd beta  = params.segment(n_alpha, m_p);
        const VectorXd gamma = (n_gamma > 0) ? params.tail(n_gamma) : VectorXd(0);

        std::vector<double> score_vals;
        MatrixXd dscore_dgamma(m_K, n_gamma);
        compute_scores(gamma, score_vals, dscore_dgamma);

        MatrixXd I = MatrixXd::Zero(d, d);
        double* I_data = I.data();
        std::vector<double> z(d), mean_grad(d), logits(m_K);

        for (int i = 0; i < m_n; ++i) {
            const double eta = (m_p > 0) ? m_X.row(i).dot(beta) : 0.0;
            const double* xi = m_X.data() + i; // column-major: xi[j*m_n] = X(i,j)

            // Compute logits and softmax probabilities
            logits[0] = 0.0;
            for (int cat = 1; cat < m_K; ++cat)
                logits[cat] = params[cat - 1] + score_vals[cat] * eta;
            double max_l = *std::max_element(logits.begin(), logits.end());
            std::vector<double> probs(m_K);
            double sp = 0.0;
            for (int cat = 0; cat < m_K; ++cat) { probs[cat] = std::exp(logits[cat] - max_l); sp += probs[cat]; }
            for (int cat = 0; cat < m_K; ++cat) probs[cat] /= sp;

            // Accumulate mean_grad = Σ_j p_j * z_j
            for (int jj = 0; jj < d; ++jj) mean_grad[jj] = 0.0;
            for (int cat = 1; cat < m_K; ++cat) {
                const double pj = probs[cat];
                const double sv = score_vals[cat];
                mean_grad[cat - 1] += pj;
                for (int b = 0; b < m_p; ++b)
                    mean_grad[n_alpha + b] += pj * sv * xi[b * m_n];
                for (int r = 0; r < n_gamma; ++r)
                    mean_grad[n_alpha + m_p + r] += pj * eta * dscore_dgamma(cat, r);
            }

            // Accumulate Σ_j p_j * z_j * z_j^T (upper triangle)
            for (int cat = 1; cat < m_K; ++cat) {
                const double pj = probs[cat];
                const double sv = score_vals[cat];
                for (int jj = 0; jj < d; ++jj) z[jj] = 0.0;
                z[cat - 1] = 1.0;
                for (int b = 0; b < m_p; ++b) z[n_alpha + b] = sv * xi[b * m_n];
                for (int r = 0; r < n_gamma; ++r) z[n_alpha + m_p + r] = eta * dscore_dgamma(cat, r);
                for (int c = 0; c < d; ++c) {
                    if (z[c] == 0.0) continue;
                    const double pz = pj * z[c];
                    for (int r = 0; r <= c; ++r)
                        I_data[r + c * d] += pz * z[r];
                }
            }

            // Subtract mean_grad * mean_grad^T (upper triangle)
            for (int c = 0; c < d; ++c) {
                if (mean_grad[c] == 0.0) continue;
                const double mg = mean_grad[c];
                for (int r = 0; r <= c; ++r)
                    I_data[r + c * d] -= mg * mean_grad[r];
            }
        }

        // Reflect upper triangle to lower
        for (int c = 0; c < d; ++c)
            for (int r = 0; r < c; ++r)
                I_data[c + r * d] = I_data[r + c * d];
        return I;
    }
};

static MatrixXd numeric_hessian_from_gradient(
    const StereotypeLogitRegression& model,
    const Eigen::Ref<const VectorXd>& params,
    double h = 1e-5
) {
    (void)h;
    return model.loglik_hessian(params);
}

static MatrixXd pseudo_inverse_symmetric(const Eigen::Ref<const MatrixXd>& A, double tol = 1e-8) {
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    VectorXd sing = svd.singularValues();
    MatrixXd Dinv = MatrixXd::Zero(sing.size(), sing.size());
    for (int i = 0; i < sing.size(); ++i) {
        if (sing[i] > tol) {
            Dinv(i, i) = 1.0 / sing[i];
        }
    }
    return svd.matrixV() * Dinv * svd.matrixU().transpose();
}

static VectorXd set_beta_and_pack_nuisance(
    const Eigen::Ref<const VectorXd>& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    const int d = n_alpha + p + n_gamma;
    VectorXd params(d);
    int cursor = 0;

    if (n_alpha > 0) {
        params.head(n_alpha) = nuisance.segment(cursor, n_alpha);
        cursor += n_alpha;
    }

    for (int j = 0; j < p; ++j) {
        if (j == beta_index) {
            params[n_alpha + j] = beta_fixed;
        } else {
            params[n_alpha + j] = nuisance[cursor];
            ++cursor;
        }
    }

    if (n_gamma > 0) {
        params.tail(n_gamma) = nuisance.tail(n_gamma);
    }

    return params;
}

static VectorXd extract_nuisance_params(
    const Eigen::Ref<const VectorXd>& params,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    VectorXd nuisance(params.size() - 1);
    int cursor = 0;

    if (n_alpha > 0) {
        nuisance.segment(cursor, n_alpha) = params.head(n_alpha);
        cursor += n_alpha;
    }

    for (int j = 0; j < p; ++j) {
        if (j == beta_index) {
            continue;
        }
        nuisance[cursor] = params[n_alpha + j];
        ++cursor;
    }

    if (n_gamma > 0) {
        nuisance.tail(n_gamma) = params.tail(n_gamma);
    }

    return nuisance;
}

static double nuisance_loglik_grad(
    const StereotypeLogitRegression& model,
    const Eigen::Ref<const VectorXd>& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    VectorXd* grad = NULL,
    int beta_index = 0
) {
    VectorXd params = set_beta_and_pack_nuisance(nuisance, beta_fixed, n_alpha, p, n_gamma, beta_index);
    VectorXd full_grad(params.size());
    double ll = model.loglik_grad(params, grad == NULL ? NULL : &full_grad);
    if (grad != NULL) {
        *grad = extract_nuisance_params(full_grad, n_alpha, p, n_gamma, beta_index);
    }
    return ll;
}

static MatrixXd nuisance_numeric_hessian_from_gradient(
    const StereotypeLogitRegression& model,
    const Eigen::Ref<const VectorXd>& nuisance,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    double h = 1e-5,
    int beta_index = 0
) {
    (void)h;
    VectorXd params = set_beta_and_pack_nuisance(nuisance, beta_fixed, n_alpha, p, n_gamma, beta_index);
    MatrixXd H_full = model.loglik_hessian(params);
    const int d = nuisance.size();
    MatrixXd H(d, d);
    std::vector<int> keep;
    keep.reserve(d);
    for (int j = 0; j < params.size(); ++j) {
        if (j != n_alpha + beta_index) {
            keep.push_back(j);
        }
    }
    for (int r = 0; r < d; ++r) {
        for (int c = 0; c < d; ++c) {
            H(r, c) = H_full(keep[r], keep[c]);
        }
    }
    return H;
}

static VectorXd optimize_nuisance_given_beta(
    const StereotypeLogitRegression& model,
    const Eigen::Ref<const VectorXd>& full_params_start,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int maxit = 50,
    double tol = 1e-8,
    int beta_index = 0
) {
    VectorXd nuisance = extract_nuisance_params(full_params_start, n_alpha, p, n_gamma, beta_index);
    const int d = nuisance.size();
    if (d == 0) {
        return nuisance;
    }

    VectorXd grad(d);
    for (int iter = 0; iter < maxit; ++iter) {
        double current_ll = nuisance_loglik_grad(model, nuisance, beta_fixed, n_alpha, p, n_gamma, &grad, beta_index);
        if (grad.norm() < tol) {
            break;
        }

        MatrixXd H = nuisance_numeric_hessian_from_gradient(
            model, nuisance, beta_fixed, n_alpha, p, n_gamma, 1e-5, beta_index
        );
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) {
            break;
        }

        VectorXd step = lu.solve(grad);
        double scale = 1.0;
        bool accepted = false;
        while (scale > 1e-8) {
            VectorXd next_nuisance = nuisance - scale * step;
            double next_ll = nuisance_loglik_grad(
                model, next_nuisance, beta_fixed, n_alpha, p, n_gamma, NULL, beta_index
            );
            if (R_finite(next_ll) && next_ll > current_ll) {
                nuisance = next_nuisance;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }

        if (!accepted || (scale * step).norm() < tol) {
            break;
        }
    }

    return nuisance;
}

static double profile_loglik_for_beta(
    const StereotypeLogitRegression& model,
    const Eigen::Ref<const VectorXd>& full_params_start,
    double beta_fixed,
    int n_alpha,
    int p,
    int n_gamma,
    int beta_index = 0
) {
    VectorXd nuisance = optimize_nuisance_given_beta(
        model, full_params_start, beta_fixed, n_alpha, p, n_gamma, 50, 1e-8, beta_index
    );
    VectorXd params = set_beta_and_pack_nuisance(nuisance, beta_fixed, n_alpha, p, n_gamma, beta_index);
    return model.loglik_grad(params, NULL);
}

static VectorXd stereotype_newton_fit(
    const StereotypeLogitRegression& model,
    int maxit,
    double tol,
    bool* converged = NULL
) {
    VectorXd params = model.initialize_params();
    const int d = params.size();
    VectorXd grad(d);
    bool did_converge = false;
    const double score_tol = std::max(tol, std::sqrt(std::max(tol, 0.0)));

    auto score_is_small = [&](const VectorXd& g) {
        return g.allFinite() && (g.norm() / std::sqrt((double)std::max(1, d))) <= score_tol;
    };

    if (converged != NULL) {
        *converged = false;
    }

    for (int iter = 0; iter < maxit; ++iter) {
        double current_ll = model.loglik_grad(params, &grad);
        if (!R_finite(current_ll) || !grad.allFinite()) {
            break;
        }
        if (score_is_small(grad)) {
            did_converge = true;
            break;
        }

        MatrixXd H = numeric_hessian_from_gradient(model, params);
        FullPivLU<MatrixXd> lu(H);
        if (!lu.isInvertible()) {
            break;
        }

        VectorXd step = lu.solve(grad);
        double scale = 1.0;
        bool accepted = false;

        while (scale > 1e-8) {
            VectorXd next_params = params - scale * step;
            double next_ll = model.loglik_grad(next_params, NULL);
            if (R_finite(next_ll) && next_ll > current_ll) {
                params = next_params;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }

        if (!accepted) {
            model.loglik_grad(params, &grad);
            did_converge = score_is_small(grad);
            break;
        }

        model.loglik_grad(params, &grad);
        if (score_is_small(grad)) {
            did_converge = true;
            break;
        }

        if ((scale * step).norm() < tol) {
            break;
        }
    }

    if (!did_converge) {
        model.loglik_grad(params, &grad);
        did_converge = score_is_small(grad);
    }
    if (converged != NULL) {
        *converged = did_converge;
    }

    return params;
}

struct StereotypeObjective {
    const StereotypeLogitRegression& model;
    StereotypeObjective(const StereotypeLogitRegression& m) : model(m) {}

    double operator()(const Eigen::Ref<const VectorXd>& params, Eigen::Ref<VectorXd> grad) const {
        VectorXd g = VectorXd::Zero(params.size());
        double val = -model.loglik_grad(params, &g);
        grad = -g;
        return val;
    }

    MatrixXd hessian(const Eigen::Ref<const VectorXd>& params) const {
        return -model.loglik_hessian(params);
    }

    MatrixXd expected_hessian(const Eigen::Ref<const VectorXd>& params) const {
        return model.expected_hessian(params);
    }
};

//' @title Compute Stereotype Logit Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a stereotype logit model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters.
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
SEXP get_stereotype_logit_score_cpp(const Rcpp::NumericMatrix& X,
											   const Rcpp::NumericVector& y,
											   const Rcpp::NumericVector& params) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());

	StereotypeLogitRegression model(map_X, map_y);
	Eigen::VectorXd grad(map_params.size());
	model.loglik_grad(map_params, &grad);
	return wrap(grad);
}

//' @title Compute Stereotype Logit Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a stereotype logit model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters.
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
SEXP get_stereotype_logit_hessian_cpp(const Rcpp::NumericMatrix& X,
												 const Rcpp::NumericVector& y,
												 const Rcpp::NumericVector& params) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());

	StereotypeLogitRegression model(map_X, map_y);
	return wrap(model.loglik_hessian(map_params));
}

//' @title Fast Stereotype Logit Regression (C++)
//' @description High-performance stereotype logit regression fitting using Newton-Raphson.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, thresholds, scores, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_stereotype_logit_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, int maxit = 100, double tol = 1e-8,
                                bool smart_cold_start = true,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "newton_raphson",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                bool estimate_only = false) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    StereotypeLogitRegression model(map_X, map_y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    VectorXd params = model.initialize_params();
    int n_par = params.size();
    int n_alpha = model.num_alpha();
    int p = map_X.cols();
    
    if (warm_start_params.isNotNull()) {
        params = as<VectorXd>(warm_start_params);
        if (params.size() != n_par) stop("warm_start_params size mismatch");
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == p) {
            params.segment(n_alpha, p) = sb;
        }
    }
    // smart_cold_start: alpha from initialize_params() empirical log-ratios; beta/gamma = 0

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);
    StereotypeObjective obj(model);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, info_start_ptr);
    params = fit.params;

    if (estimate_only) {
        return List::create(
            Named("b") = params.segment(n_alpha, p),
            Named("alpha") = params.head(n_alpha),
            Named("scores_raw") = (model.num_gamma() > 0) ? params.tail(model.num_gamma()) : VectorXd(0),
            Named("params") = params,
            Named("converged") = fit.converged
        );
    }
    return List::create(
        Named("b") = params.segment(n_alpha, p),
        Named("alpha") = params.head(n_alpha),
        Named("scores_raw") = (model.num_gamma() > 0) ? params.tail(model.num_gamma()) : VectorXd(0),
        Named("params") = params,
        Named("converged") = fit.converged,
        Named("fisher_information") = (-model.loglik_hessian(params))
    );
}

//' @title Fast Stereotype Logit Regression with Variance (C++)
//' @description Stereotype logit regression fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @return A list containing coefficients, variance estimates, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_stereotype_logit_with_var_cpp(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& y, int maxit = 100, double tol = 1e-8,
                                         bool smart_cold_start = true,
                                         Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                         Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                         std::string optimization_alg = "newton_raphson",
                                         Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                         Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
                                         Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                         bool estimate_only = false) {

    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    StereotypeLogitRegression model(map_X, map_y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    VectorXd params = model.initialize_params();
    int n_par = params.size();
    
    if (warm_start_params.isNotNull()) {
        params = as<VectorXd>(warm_start_params);
        if (params.size() != n_par) stop("warm_start_params size mismatch");
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == map_X.cols()) {
            params.segment(model.num_alpha(), map_X.cols()) = sb;
        }
    }
    // smart_cold_start: alpha from initialize_params() empirical log-ratios; beta/gamma = 0

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);
    StereotypeObjective obj(model);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, info_start_ptr);
    params = fit.params;
    MatrixXd H = model.loglik_hessian(params);
    MatrixXd info = -H;

    int n_alpha = model.num_alpha();
    int p = map_X.cols();

    MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
    int free_j = -1;
    for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
        if (fixed_spec.free_idx[jj] == n_alpha) { free_j = jj + 1; break; }
    double ssq_b_1 = (p >= 1 && free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;
    
    // Fallback profiling logic for variance if still not finite
    if (!R_finite(ssq_b_1) && p >= 1 && !fixed_spec.has_fixed) {
        double beta_hat = params[n_alpha];
        double h = std::max(1e-4, 1e-3 * (std::abs(beta_hat) + 1.0));
        double ll_0 = profile_loglik_for_beta(model, params, beta_hat, n_alpha, p, model.num_gamma(), 0);
        double ll_p = profile_loglik_for_beta(model, params, beta_hat + h, n_alpha, p, model.num_gamma(), 0);
        double ll_m = profile_loglik_for_beta(model, params, beta_hat - h, n_alpha, p, model.num_gamma(), 0);
        double info_beta = -(ll_p - 2.0 * ll_0 + ll_m) / (h * h);
        if (R_finite(info_beta) && info_beta > 0) {
            ssq_b_1 = 1.0 / info_beta;
        }
    }

    return List::create(
        Named("b") = params.segment(n_alpha, p),
        Named("alpha") = params.head(n_alpha),
        Named("params") = params,
        Named("ssq_b_1") = ssq_b_1,
        Named("ssq_b_j") = ssq_b_1,
        Named("vcov") = R_NilValue,
        Named("converged") = fit.converged,
        Named("fisher_information") = info
    );
}

//' @title Compute Stereotype Profile Log-Likelihood (C++)
//' @description Calculates the profile log-likelihood for a fixed beta in a stereotype logit model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param beta_fixed The fixed value for the first beta parameter.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @return The profile log-likelihood value.
//' @export
//' @keywords internal
// [[Rcpp::export]]
double fast_stereotype_profile_loglik_cpp(
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& y,
    double beta_fixed,
    int maxit = 100,
    double tol = 1e-8,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue
) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    StereotypeLogitRegression model(map_X, map_y);
    if (model.num_categories() < 2) {
        stop("Stereotype logistic regression requires at least two observed outcome categories.");
    }

    int n_alpha = model.num_alpha();
    int p = map_X.cols();
    int n_gamma = model.num_gamma();
    int total = n_alpha + p + n_gamma;

    VectorXd params = model.initialize_params();
    if (warm_start_params.isNotNull()) {
        params = as<VectorXd>(warm_start_params);
        if (params.size() != total) stop("warm_start_params size mismatch");
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == total) {
            params = sb;
        } else if (sb.size() == p) {
            params.segment(n_alpha, p) = sb;
        }
    }

    return profile_loglik_for_beta(model, params, beta_fixed, n_alpha, p, n_gamma, 0);
}

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

//' Parallel Stereotype Logit Randomization Distribution
//'
//' @param X Matrix of covariates (without intercept or treatment).
//' @param y Numeric vector of response values (pre-null-shifted for treated).
//' @param w_mat Integer matrix of permuted treatment assignments (n x nsim).
//' @param delta Null treatment effect (additive shift).
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length nsim with treatment coefficients.
// [[Rcpp::export]]
NumericVector compute_stereotype_logit_distr_parallel_cpp(
	const Rcpp::NumericMatrix& X,
	const Rcpp::NumericVector& y,
	const Rcpp::IntegerMatrix& w_mat,
	double delta,
	int num_cores
) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

	int nsim = w_mat.cols();
	int n = map_y.size();
	int p_covars = map_X.cols();
	int p_full = p_covars + 1; // treatment + covars (no intercept — thresholds handle location)

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
				X_full(i, 1 + k) = map_X(i, k);
			}
			y_shifted[i] = (w_col[i] == 1) ? map_y[i] + delta : map_y[i];
		}


		StereotypeLogitRegression model(X_full, y_shifted);
		if (model.num_categories() < 2) continue;

		bool converged = false;
		Eigen::VectorXd params = stereotype_newton_fit(model, 100, 1e-8, &converged);

		int n_alpha = model.num_alpha();
		// accept result even if not formally converged, matching R generate_mod behaviour
		if ((int)params.size() >= n_alpha + 1 && std::isfinite(params[n_alpha])) {
			results[b] = params[n_alpha];
		}
	}

	return wrap(results);
}
