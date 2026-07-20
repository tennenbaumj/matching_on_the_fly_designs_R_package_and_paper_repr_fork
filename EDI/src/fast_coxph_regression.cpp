#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct CoxData {
    Eigen::VectorXd y;
    Eigen::VectorXd dead;
    RowMajorMatrixXd X;
    std::vector<double> unique_event_times;
    std::vector<int> event_counts;
    std::vector<int> input_idx;
    int n;
    int p;

    template<typename Derived>
    CoxData(const Eigen::VectorXd& y_in, const Eigen::VectorXd& dead_in, const Eigen::MatrixBase<Derived>& X_in) :
        n((int)y_in.size()), p((int)X_in.cols()) {

        input_idx.resize(n);
        std::iota(input_idx.begin(), input_idx.end(), 0);
        std::sort(input_idx.begin(), input_idx.end(), [&](int i, int j) {
            if (y_in[i] == y_in[j]) return dead_in[i] > dead_in[j];
            return y_in[i] < y_in[j];
        });

        y.resize(n);
        dead.resize(n);
        X.resize(n, p);

        for (int i = 0; i < n; ++i) {
            int id = input_idx[i];
            y[i] = y_in[id];
            dead[i] = dead_in[id];
            X.row(i) = X_in.row(id);
        }

        for (int i = 0; i < n; ++i) {
            if (dead[i] > 0.5) {
                if (unique_event_times.empty() || y[i] > unique_event_times.back()) {
                    unique_event_times.push_back(y[i]);
                    event_counts.push_back(0);
                }
            }
        }

        int k = -1;
        for (int i = 0; i < n; ++i) {
            if (k + 1 < (int)unique_event_times.size() && y[i] == unique_event_times[k+1])
                k++;
            if (k >= 0 && y[i] == unique_event_times[k] && dead[i] > 0.5)
                event_counts[k]++;
        }
    }

    inline const double* row(int i) const { return X.row(i).data(); }
    inline Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix_map() const {
        return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(X.data(), n, p);
    }
};

struct CoxWorkspace {
    Eigen::VectorXd eta;
    Eigen::VectorXd exp_eta;
    Eigen::VectorXd S1;
    Eigen::MatrixXd S2;
    Eigen::VectorXd sum_x_dk;
    Eigen::VectorXd e_z;
    Eigen::VectorXd grad;
    Eigen::MatrixXd hess;

    CoxWorkspace() = default;

    CoxWorkspace(int n, int p) :
        eta(n),
        exp_eta(n),
        S1(p),
        S2(p, p),
        sum_x_dk(p),
        e_z(p),
        grad(p),
        hess(p, p) {}
};

double compute_cox_neg_ll_only(const CoxData& data, const std::vector<double>& beta, CoxWorkspace& workspace) {
    const int n = data.n;
    const int p = data.p;
    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);
    workspace.eta.noalias() = data.matrix_map() * beta_map;
    const double max_eta = workspace.eta.maxCoeff();
    workspace.exp_eta.array() = (workspace.eta.array() - max_eta).exp();
    double r_exp = workspace.exp_eta.sum();
    double neg_ll = 0.0;
    int j = 0;
    for (size_t k = 0; k < data.unique_event_times.size(); ++k) {
        const double tk = data.unique_event_times[k];
        while (j < n && data.y[j] < tk) {
            r_exp -= workspace.exp_eta[j];
            ++j;
        }
        double sum_eta_dk = 0.0;
        const int count_dk = data.event_counts[k];
        for (int m = j; m < n && data.y[m] == tk; ++m) {
            if (data.dead[m] > 0.5) sum_eta_dk += workspace.eta[m];
        }
        if (count_dk > 0) {
            const double safe_r = std::max(r_exp, 1e-100);
            neg_ll -= (sum_eta_dk - count_dk * (std::log(safe_r) + max_eta));
        }
    }
    return neg_ll;
}

std::vector<CoxWorkspace> make_cox_workspaces(const std::vector<CoxData>& strata_data) {
    std::vector<CoxWorkspace> workspaces;
    workspaces.reserve(strata_data.size());
    for (const CoxData& sd : strata_data) {
        workspaces.emplace_back(sd.n, sd.p);
    }
    return workspaces;
};

double compute_cox_ll_grad_hess_fast(
        const CoxData& data,
        const std::vector<double>& beta,
        Eigen::VectorXd& grad,
        Eigen::MatrixXd& hess,
        bool estimate_only,
        CoxWorkspace& workspace
) {
    const int n = data.n;
    const int p = data.p;
    
    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);
    workspace.eta.noalias() = data.matrix_map() * beta_map;
    const double max_eta = workspace.eta.maxCoeff();
    workspace.exp_eta.array() = (workspace.eta.array() - max_eta).exp();

    grad.setZero();
    if (!estimate_only) hess.setZero();

    double neg_ll = 0.0;
    double S0 = 0.0;
    workspace.S1.setZero(); 
    if (!estimate_only) workspace.S2.setZero(); 

    int data_idx = n - 1; 

    for (int k = (int)data.unique_event_times.size() - 1; k >= 0; --k) {
        double tk = data.unique_event_times[k];
        int dk = data.event_counts[k];

        int start_idx_for_tk = data_idx;
        while (data_idx >= 0 && data.y[data_idx] >= tk) {
            int id = data_idx;
            const double wi = workspace.exp_eta[id];
            const double* xi = data.row(id);
            
            S0 += wi;
            double* r_x_exp_ptr = workspace.S1.data();
            for (int q = 0; q < p; ++q) {
                r_x_exp_ptr[q] += wi * xi[q];
            }
            if (!estimate_only) {
                double* r_xx_exp_ptr = workspace.S2.data();
                for (int q1 = 0; q1 < p; ++q1) {
                    double w_xi_q1 = wi * xi[q1];
                    for (int q2 = q1; q2 < p; ++q2) {
                        r_xx_exp_ptr[q2 + q1 * p] += w_xi_q1 * xi[q2];
                    }
                }
            }
            --data_idx;
        }
        
        double sum_eta_dk = 0.0;
        workspace.sum_x_dk.setZero();
        double* sum_x_dk_ptr = workspace.sum_x_dk.data();
        
        for (int id = data_idx + 1; id <= start_idx_for_tk; ++id) {
            if (data.dead[id] > 0.5) {
                sum_eta_dk += workspace.eta[id];
                const double* xi = data.row(id);
                for (int q = 0; q < p; ++q) {
                    sum_x_dk_ptr[q] += xi[q];
                }
            }
        }

        if (dk > 0 && S0 > 0) {
            double safe_r = std::max(S0, 1e-100);
            double inv_r = 1.0 / safe_r;
            neg_ll -= (sum_eta_dk - dk * (std::log(safe_r) + max_eta));
            
            workspace.e_z.noalias() = workspace.S1 * inv_r;
            grad.noalias() -= (workspace.sum_x_dk - dk * workspace.e_z);
            if (!estimate_only) {
                for (int q1 = 0; q1 < p; ++q1) {
                    const double ez_q1 = workspace.e_z[q1];
                    for (int q2 = q1; q2 < p; ++q2) {
                        const double contribution = dk * (
                            workspace.S2(q2, q1) * inv_r - workspace.e_z[q2] * ez_q1
                        );
                        hess(q2, q1) += contribution;
                        if (q2 != q1) hess(q1, q2) += contribution;
                    }
                }
            }
        }
    }

    return neg_ll;
}

struct CoxFitResult {
    std::vector<double> beta;
    Eigen::MatrixXd vcov;
    Eigen::MatrixXd hess_mat;
    double neg_ll;
    bool converged;
    int iterations;

    CoxFitResult() : neg_ll(NA_REAL), converged(false), iterations(0) {}
};

CoxFitResult cox_newton_raphson(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    if (strata_data.empty()) return CoxFitResult();
    const int p = (int)strata_data[0].p;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        NumericVector sb(warm_start_beta);
        for (int q = 0; q < p; ++q) beta[q] = sb[q];
    } else if (smart_cold_start) {
        int total_n = 0;
        for (const auto& sd : strata_data) total_n += sd.n;
        Eigen::MatrixXd X_full(total_n, p);
        Eigen::VectorXd log_y(total_n);
        int offset = 0;
        for (const auto& sd : strata_data) {
            X_full.block(offset, 0, sd.n, p) = sd.X;
            for (int i = 0; i < sd.n; ++i) log_y[offset + i] = std::log(std::max(sd.y[i], 1e-8));
            offset += sd.n;
        }
        beta = -safe_ols_solve(X_full, log_y);
    }
    beta = apply_fixed_values(beta, fixed_spec);

    std::vector<CoxWorkspace> workspaces = make_cox_workspaces(strata_data);
    Eigen::VectorXd beta_candidate(p);
    std::vector<double> beta_cand_vec(p);
    std::vector<double> beta_vec(p);
    
    Eigen::VectorXd total_grad(p);
    Eigen::MatrixXd total_hess(p, p);

    double old_ll = 1e300;
    int iter = 0;

    for (iter = 0; iter < maxit; ++iter) {
        total_grad.setZero();
        total_hess.setZero();
        double ll = 0.0;

        for (int q = 0; q < p; ++q) beta_vec[q] = beta[q];

        for (std::size_t s = 0; s < strata_data.size(); ++s) {
            const CoxData& sd = strata_data[s];
            CoxWorkspace& ws = workspaces[s];
            ll += compute_cox_ll_grad_hess_fast(sd, beta_vec, ws.grad, ws.hess, false, ws);
            total_grad.noalias() += ws.grad;
            total_hess.noalias() += ws.hess;
        }

        if (std::abs(old_ll - ll) < tol) break;
        old_ll = ll;

        Eigen::MatrixXd H;
        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            H = as<Eigen::MatrixXd>(warm_start_fisher_info);
            if (H.rows() != p || H.cols() != p) {
                H = total_hess;
            }
        } else {
            H = total_hess;
        }

        if (!H.allFinite() || !total_grad.allFinite()) break;

        for (int i = 0; i < (int)fixed_spec.fixed_idx.size(); ++i) {
            int idx = fixed_spec.fixed_idx[i];
            total_grad[idx] = 0.0;
            H.row(idx).setZero();
            H.col(idx).setZero();
            H(idx, idx) = 1.0;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) break;
        const Eigen::VectorXd delta = ldlt.solve(total_grad);
        if (!delta.allFinite()) break;

        // Step-halving line search: ensure the neg log-likelihood decreases.
        static const int kMaxHalvings = 10;
        auto eval_ll_candidate = [&](double step) -> double {
            for (int q = 0; q < p; ++q) beta_candidate[q] = beta[q] - step * delta[q];
            for (int fi = 0; fi < (int)fixed_spec.fixed_idx.size(); ++fi)
                beta_candidate[fixed_spec.fixed_idx[fi]] = fixed_spec.fixed_values[fi];
            for (int q = 0; q < p; ++q) beta_cand_vec[q] = beta_candidate[q];
            double ll_c = 0.0;
            for (std::size_t s = 0; s < strata_data.size(); ++s)
                ll_c += compute_cox_neg_ll_only(strata_data[s], beta_cand_vec, workspaces[s]);
            return ll_c;
        };

        double step = 1.0;
        double ll_candidate = eval_ll_candidate(step);
        for (int h = 0; h < kMaxHalvings && ll_candidate >= ll; ++h) {
            step *= 0.5;
            ll_candidate = eval_ll_candidate(step);
        }
        beta = beta_candidate;
    }

    CoxFitResult res;
    res.beta.assign(beta.data(), beta.data() + p);
    res.neg_ll = old_ll;
    res.converged = (iter < maxit);
    res.iterations = iter;
    res.hess_mat = total_hess;

    if (!estimate_only) {
        Eigen::MatrixXd H_free = subset_matrix(res.hess_mat, fixed_spec.free_idx, fixed_spec.free_idx);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            Eigen::MatrixXd vcov_free = lu.inverse();
            res.vcov = expand_free_covariance(p, fixed_spec, vcov_free, true);
        } else {
            res.vcov = Eigen::MatrixXd::Constant(p, p, NA_REAL);
        }
    }

    return res;
}

class StratifiedCoxObjective {
    const std::vector<CoxData>& m_strata;
    const int m_p;
    mutable std::vector<CoxWorkspace> m_workspaces;
public:
    StratifiedCoxObjective(const std::vector<CoxData>& strata, int p)
        : m_strata(strata), m_p(p), m_workspaces(make_cox_workspaces(strata)) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        double nll = 0.0;
        grad.setZero();
        for (std::size_t s = 0; s < m_strata.size(); ++s) {
            const CoxData& sd = m_strata[s];
            CoxWorkspace& ws = m_workspaces[s];
            nll += compute_cox_ll_grad_hess_fast(sd, beta, ws.grad, ws.hess, true, ws);
            grad.noalias() += ws.grad;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        Eigen::MatrixXd h(m_p, m_p);
        h.setZero();
        for (std::size_t s = 0; s < m_strata.size(); ++s) {
            const CoxData& sd = m_strata[s];
            CoxWorkspace& ws = m_workspaces[s];
            compute_cox_ll_grad_hess_fast(sd, beta, ws.grad, ws.hess, false, ws);
            h.noalias() += ws.hess;
        }
        return h;
    }
};

CoxFitResult cox_lbfgs(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    if (strata_data.empty()) return CoxFitResult();
    const int p = (int)strata_data[0].p;
    Eigen::VectorXd par = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        NumericVector sb(warm_start_beta);
        for (int q = 0; q < p; ++q) par[q] = sb[q];
    } else if (smart_cold_start) {
        int total_n = 0;
        for (const auto& sd : strata_data) total_n += sd.n;
        Eigen::MatrixXd X_full(total_n, p);
        Eigen::VectorXd log_y(total_n);
        int offset = 0;
        for (const auto& sd : strata_data) {
            X_full.block(offset, 0, sd.n, p) = sd.X;
            for (int i = 0; i < sd.n; ++i) log_y[offset + i] = std::log(std::max(sd.y[i], 1e-8));
            offset += sd.n;
        }
        par = -safe_ols_solve(X_full, log_y);
    }
    par = apply_fixed_values(par, fixed_spec);

    StratifiedCoxObjective obj(strata_data, p);
    Eigen::MatrixXd info_start;
    Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, tol, "lbfgs", "lbfgs", 0, info_start_ptr);

    CoxFitResult res;
    res.beta.assign(fit.params.data(), fit.params.data() + p);
    res.neg_ll    = fit.value;
    res.converged = fit.converged;
    res.iterations = fit.niter;

    if (!estimate_only && (fit.converged || true)) { // Always try to get Hessian if not estimate_only
        res.hess_mat = obj.hessian(fit.params);
        Eigen::MatrixXd H_free = subset_matrix(res.hess_mat, fixed_spec.free_idx, fixed_spec.free_idx);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            Eigen::MatrixXd vcov_free = lu.inverse();
            res.vcov = expand_free_covariance(p, fixed_spec, vcov_free, true);
        } else {
            res.vcov = Eigen::MatrixXd::Constant(p, p, NA_REAL);
        }
    }
    return res;
}

CoxFitResult cox_fit(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    const std::string& optimization_alg,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    if (optimization_alg == "lbfgs")
        return cox_lbfgs(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol);
    return cox_newton_raphson(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, warm_start_fisher_info);
}

Eigen::MatrixXd compute_robust_vcov(
    const std::vector<CoxData>& strata_data,
    const std::vector<double>& beta,
    const Eigen::MatrixXd& H_inv,
    const std::vector<int>& cluster)
{
    int n_total = 0;
    for (const CoxData& sd : strata_data) n_total += sd.n;
    const int p = (int)beta.size();
    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);

    Eigen::MatrixXd U(n_total, p);
    U.setZero();

    int row_offset = 0;
    for (const CoxData& sd : strata_data) {
        const int ns = sd.n;
        Eigen::VectorXd eta = sd.X * beta_map;
        Eigen::VectorXd exp_eta = (eta.array() - eta.maxCoeff()).exp().matrix();

        double r_exp = 0.0;
        Eigen::VectorXd r_x_exp = Eigen::VectorXd::Zero(p);
        for (int i = 0; i < ns; ++i) {
            r_exp += exp_eta[i];
            for (int q = 0; q < p; ++q) r_x_exp[q] += sd.X(i, q) * exp_eta[i];
        }

        const int n_events = (int)sd.unique_event_times.size();
        Eigen::VectorXd dk_over_Rk = Eigen::VectorXd::Zero(n_events);
        Eigen::MatrixXd ek = Eigen::MatrixXd::Zero(p, n_events);
        Eigen::MatrixXd dk_ek_over_Rk = Eigen::MatrixXd::Zero(p, n_events);

        int j = 0;
        for (int k = 0; k < n_events; ++k) {
            double tk = sd.unique_event_times[k];
            while (j < ns && sd.y[j] < tk) {
                int id = j;
                r_exp -= exp_eta[id];
                for (int q = 0; q < p; ++q) r_x_exp[q] -= sd.X(id, q) * exp_eta[id];
                ++j;
            }
            int dk = sd.event_counts[k];
            if (dk > 0) {
                double safe_r = std::max(r_exp, 1e-100);
                dk_over_Rk[k] = (double)dk / safe_r;
                for (int q = 0; q < p; ++q) {
                    ek(q, k) = r_x_exp[q] / safe_r;
                    dk_ek_over_Rk(q, k) = (double)dk * ek(q, k) / safe_r;
                }
            }
        }

        Eigen::VectorXd cum_A = Eigen::VectorXd::Zero(n_events);
        Eigen::MatrixXd cum_B = Eigen::MatrixXd::Zero(p, n_events);
        if (n_events > 0) {
            cum_A[0] = dk_over_Rk[0];
            cum_B.col(0) = dk_ek_over_Rk.col(0);
            for (int k = 1; k < n_events; ++k) {
                cum_A[k] = cum_A[k-1] + dk_over_Rk[k];
                cum_B.col(k) = cum_B.col(k - 1) + dk_ek_over_Rk.col(k);
            }
        }

        for (int i = 0; i < ns; ++i) {
            double yi = sd.y[i];
            double ei = exp_eta[i];
            double di = sd.dead[i];

            int k_last = -1;
            if (n_events > 0) {
                int lo = 0, hi = n_events - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) / 2;
                    if (sd.unique_event_times[mid] <= yi) { k_last = mid; lo = mid + 1; }
                    else hi = mid - 1;
                }
            }
            double A_i = (k_last >= 0) ? cum_A[k_last] : 0.0;
            for (int q = 0; q < p; ++q) {
                double B_iq    = (k_last >= 0) ? cum_B(q, k_last) : 0.0;
                double e_yi_q  = (di > 0.5 && k_last >= 0) ? ek(q, k_last) : 0.0;
                U(row_offset + sd.input_idx[i], q) =
                    (sd.X(i, q) - e_yi_q) * di - ei * (sd.X(i, q) * A_i - B_iq);
            }
        }
        row_offset += ns;
    }

    std::unordered_map<int, Eigen::VectorXd> cluster_scores;
    cluster_scores.reserve(cluster.size());
    for (int i = 0; i < n_total; ++i) {
        int c = cluster[i];
        auto [it, inserted] = cluster_scores.try_emplace(c, U.row(i).transpose());
        if (!inserted) it->second += U.row(i).transpose();
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(p, p);
    for (auto& kv : cluster_scores) B += kv.second * kv.second.transpose();
    return H_inv * B * H_inv;
}

} // namespace

// [[Rcpp::export]]
SEXP build_cox_data_cache_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead) {
    auto* data = new std::vector<CoxData>();
    data->emplace_back(y, dead, X);
    return Rcpp::XPtr<std::vector<CoxData>>(data, true);
}

// [[Rcpp::export]]
SEXP build_stratified_cox_data_cache_cpp(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Rcpp::IntegerVector& strata)
{
    const int n = (int)y.size();
    const int p = (int)X.cols();
    std::map<int, std::vector<int>> strata_map;
    for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);
    RowMajorMatrixXd X_rm = X;
    auto* data = new std::vector<CoxData>();
    data->reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        RowMajorMatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            int id = idx[ii];
            y_s[ii] = y[id]; dead_s[ii] = dead[id]; X_s.row(ii) = X_rm.row(id);
        }
        data->emplace_back(y_s, dead_s, X_s);
    }
    return Rcpp::XPtr<std::vector<CoxData>>(data, true);
}

// [[Rcpp::export]]
List fast_coxph_regression_prebuilt_cpp(
    SEXP cox_data_xptr,
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    bool estimate_only = false,
    int maxit = 20,
    double tol = 1e-9,
    Nullable<IntegerVector> fixed_idx = R_NilValue,
    Nullable<NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "newton_raphson",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    Rcpp::XPtr<std::vector<CoxData>> data_ptr(cox_data_xptr);
    const std::vector<CoxData>& strata_data = *data_ptr;
    if (strata_data.empty()) return List::create();
    const int p = strata_data[0].p;
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    return List::create(_["coefficients"] = coef_r, _["vcov"] = fit.vcov, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
List fast_coxph_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead,
                               Nullable<NumericVector> warm_start_beta = R_NilValue,
                               bool smart_cold_start = true,
                               bool estimate_only = false,
                               int maxit = 20,
                               double tol = 1e-9,
                               Nullable<IntegerVector> cluster = R_NilValue,
                               Nullable<IntegerVector> fixed_idx = R_NilValue,
                               Nullable<NumericVector> fixed_values = R_NilValue,
                               std::string optimization_alg = "newton_raphson",
                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int p = (int)X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::vector<CoxData> strata_data;
    strata_data.emplace_back(y, dead, X);
    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    Eigen::MatrixXd vcov_mat = (cluster.isNotNull()) ? compute_robust_vcov(strata_data, fit.beta, fit.vcov, std::vector<int>(IntegerVector(cluster).begin(), IntegerVector(cluster).end())) : fit.vcov;
    return List::create(_["coefficients"] = coef_r, _["vcov"] = vcov_mat, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
List fast_stratified_coxph_regression_cpp(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Rcpp::IntegerVector& strata,
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    bool estimate_only = false,
    int maxit = 20,
    double tol = 1e-9,
    Nullable<IntegerVector> fixed_idx = R_NilValue,
    Nullable<NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "newton_raphson",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    const int n = (int)y.size();
    const int p = (int)X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    // Efficiently group by strata
    std::map<int, std::vector<int>> strata_map;
    for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);

    RowMajorMatrixXd X_rm = X;
    std::vector<CoxData> strata_data;
    strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        RowMajorMatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            int id = idx[ii];
            y_s[ii] = y[id]; dead_s[ii] = dead[id]; X_s.row(ii) = X_rm.row(id);
        }
        strata_data.emplace_back(y_s, dead_s, X_s);
    }

    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = (p > 0 && !fit.beta.empty()) ? fit.beta[q] : NA_REAL;
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    return List::create(_["coefficients"] = coef_r, _["vcov"] = fit.vcov, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
Eigen::VectorXd get_coxph_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data; strata_data.emplace_back(y, dead, X);
    CoxWorkspace ws((int)y.size(), (int)X.cols());
    Eigen::VectorXd grad(beta.size());
    Eigen::MatrixXd hess(beta.size(), beta.size());
    std::vector<double> beta_vec(beta.data(), beta.data() + beta.size());
    compute_cox_ll_grad_hess_fast(strata_data[0], beta_vec, grad, hess, true, ws);
    return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_coxph_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data; strata_data.emplace_back(y, dead, X);
    CoxWorkspace ws((int)y.size(), (int)X.cols());
    Eigen::VectorXd grad(beta.size());
    Eigen::MatrixXd hess(beta.size(), beta.size());
    std::vector<double> beta_vec(beta.data(), beta.data() + beta.size());
    compute_cox_ll_grad_hess_fast(strata_data[0], beta_vec, grad, hess, false, ws);
    return -hess;
}

// [[Rcpp::export]]
Eigen::VectorXd get_stratified_coxph_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Rcpp::IntegerVector& strata, const Eigen::VectorXd& beta) {
    const int n = (int)y.size(); const int p = (int)X.cols();
    std::map<int, std::vector<int>> strata_map; for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);
    std::vector<CoxData> strata_data; strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        int ns = (int)idx.size(); Eigen::VectorXd ys(ns), ds(ns); Eigen::MatrixXd Xs(ns, p);
        for (int ii = 0; ii < ns; ++ii) { int id = idx[ii]; ys[ii] = y[id]; ds[ii] = dead[id]; Xs.row(ii) = X.row(id); }
        strata_data.emplace_back(ys, ds, Xs);
    }
    StratifiedCoxObjective obj(strata_data, p); Eigen::VectorXd grad(beta.size()); obj(beta, grad); return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_stratified_coxph_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Rcpp::IntegerVector& strata, const Eigen::VectorXd& beta) {
    const int n = (int)y.size(); const int p = (int)X.cols();
    std::map<int, std::vector<int>> strata_map; for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);
    std::vector<CoxData> strata_data; strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        int ns = (int)idx.size(); Eigen::VectorXd ys(ns), ds(ns); Eigen::MatrixXd Xs(ns, p);
        for (int ii = 0; ii < ns; ++ii) { int id = idx[ii]; ys[ii] = y[id]; ds[ii] = dead[id]; Xs.row(ii) = X.row(id); }
        strata_data.emplace_back(ys, ds, Xs);
    }
    StratifiedCoxObjective obj(strata_data, p); return -obj.hessian(beta);
}

//' @title Bootstrap Randomization Test distribution for Cox PH (treatment-only model)
//' @description Computes the BRT null distribution of the treatment log-hazard-ratio from a
//'   treatment-only Cox PH model across B pre-generated (i_mat, w_mat) draws. Each draw
//'   resamples rows and applies a fresh treatment assignment; the sharp-null shift is
//'   multiplicative on treated survival times (exp(delta) on the log-time scale).
//' @param y0 Numeric vector of original survival times (length n).
//' @param dead Numeric vector of event indicators (length n).
//' @param i_mat Integer matrix (n x B) of 1-based row indices.
//' @param w_mat Integer matrix (n x B) of treatment assignments (0/1).
//' @param delta Sharp-null shift on log-time scale.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length B with treatment log-HR per draw (NA on non-convergence).
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_coxph_rand_bootstrap_cpp(
    const Eigen::VectorXd& y0,
    const Eigen::VectorXd& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    int num_cores)
{
    const int n  = (int)y0.size();
    const int B  = i_mat.cols();
    NumericVector out(B, NA_REAL);
    const double mult = std::exp(delta);
    const int* i_ptr = i_mat.begin();
    const int* w_ptr = w_mat.begin();

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, num_cores));
#endif

#pragma omp parallel for schedule(dynamic) if(num_cores > 1)
    for (int b = 0; b < B; ++b) {
        Eigen::VectorXd y_b(n), dead_b(n);
        RowMajorMatrixXd X_b(n, 1);
        for (int i = 0; i < n; ++i) {
            int row0 = i_ptr[b * n + i] - 1;
            int  w_i = w_ptr[b * n + i];
            double yi = y0[row0];
            if (delta != 0.0 && w_i == 1) yi *= mult;
            y_b[i]    = yi;
            dead_b[i] = dead[row0];
            X_b(i, 0) = (double)w_i;
        }
        std::vector<CoxData> strata;
        strata.emplace_back(y_b, dead_b, X_b);
        CoxFitResult fit = cox_fit(strata, R_NilValue, false, FixedParamSpec(),
                                   true, 20, 1e-9, "newton_raphson");
        if (fit.converged && !fit.beta.empty() && std::isfinite(fit.beta[0]))
            out[b] = fit.beta[0];
    }
    return out;
}

// [[Rcpp::export]]
NumericVector compute_coxph_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores)
{
    const int n      = i_mat.nrow();
    const int nsim   = i_mat.ncol();
    const int n_full = y0.size();
    const int p_cov  = Xc.ncol();
    const int p      = 1 + p_cov;  // treatment + covariates; no intercept in Cox PH

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(p);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 1) = xc_ptr[(size_t)j * n_full + r];
                double yv = y0_ptr[r];
                if (has_noise) yv += noise_ptr[(size_t)b * n + i];
                y_b(i)    = (wt && mult != 1.0) ? yv * mult : yv;
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            std::vector<CoxData> strata;
            strata.emplace_back(y_b, dead_b, X_b);
            CoxFitResult fit = cox_newton_raphson(strata, R_NilValue, true, fspec,
                                                  true, 20, 1e-6);
            if (fit.converged && !fit.beta.empty() && std::isfinite(fit.beta[0]))
                res_ptr[b] = fit.beta[0];
        }
    }
    return wrap(results);
}
