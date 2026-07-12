// Gaussian Linear Mixed Model for KK designs
// Model: y_ij = X_ij β + u_i + ε_ij
//   u_i ~ N(0, σ²_b)  (random intercept per matched pair / singleton)
//   ε_ij ~ N(0, σ²_e) (residual)
//
// Parameter vector: par = (β[0..p-1], log_σ_e, log_σ_b)
//   log_σ_e = log(σ_e),  log_σ_b = log(σ_b)
//
// Groups: matched pairs have size 2; reservoir subjects have size 1.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

namespace {

// ── Pre-computed, design-fixed quantities per group ─────────────────────────
struct GroupInfo {
    int    size;       // m_i  (1 or 2)
    int    warm_start_params;      // index into sorted-obs arrays
    double w;          // case weight for this group (1.0 = unweighted)
};

// ── Core data object (passed by const-ref to all functions) ──────────────────
struct LMMData {
    // Observations sorted so that obs within each group are contiguous
    Eigen::VectorXd y_s;         // n sorted responses
    Eigen::MatrixXd X_s;         // n × p sorted design (includes intercept at col 0)
    std::vector<GroupInfo> grps; // one entry per group
    int n;
    int p;
    int G;  // number of groups

    LMMData(const Eigen::Ref<const VectorXd>& y,
            const Eigen::Ref<const MatrixXd>& X,
            const std::vector<int>& gid,              // 0-based group ids, length n
            const std::vector<double>& obs_weights = {}) // per-obs case weights (empty = uniform 1)
        : n(y.size()), p(X.cols())
    {
        // Sort observations by group id
        std::vector<int> ord(n);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(),
                         [&](int a, int b){ return gid[a] < gid[b]; });

        y_s.resize(n);
        X_s.resize(n, p);
        for (int i = 0; i < n; ++i) {
            y_s[i]   = y[ord[i]];
            X_s.row(i) = X.row(ord[i]);
        }

        // Build group list; all obs in a matched pair share the same case weight
        const bool have_weights = ((int)obs_weights.size() == n);
        int prev = -1, gi = 0;
        for (int i = 0; i < n; ++i) {
            int g = gid[ord[i]];
            if (g != prev) {
                GroupInfo info;
                info.warm_start_params = i;
                info.size  = 1;
                info.w     = have_weights ? obs_weights[ord[i]] : 1.0;
                grps.push_back(info);
                prev = g;
                gi++;
            } else {
                grps.back().size++;
            }
        }
        G = (int)grps.size();
    }
};

// ── Neg log-likelihood with analytical gradient ──────────────────────────────
//
// neg_ll = (n/2)log(2π)
//        + Σ_g [ (m_g-1)/2 · log(v_e) + 1/2 · log(a_g) ]
//        + 1/(2 v_e) · Σ_g [ Q_g - (v_b/a_g)·S_g² ]
//
// where a_g = v_e + m_g · v_b,  S_g = Σ r,  Q_g = Σ r²,  r = y - Xβ
//
// Analytical gradients:
//   ∂/∂β       = -(1/v_e)·X^T·r + (v_b/v_e) Σ_g (S_g/a_g)·Σ_{j∈g} x_j
//   ∂/∂log_σ_e = Σ_g(m_g-1) + v_e·Σ_g(1/a_g) - (1/v_e)·Σ_g W_g
//   ∂/∂log_σ_b = Σ_g(m_g·v_b/a_g) - v_b·Σ_g(S_g²/a_g²)
//
// where W_g = Q_g - (v_b/a_g)·S_g²

// neg_ll_and_grad with caller-supplied scratch vectors (avoids per-call heap alloc).
double neg_ll_and_grad(const LMMData& dat,
                       const Eigen::Ref<const Eigen::VectorXd>& par,
                       Eigen::Ref<VectorXd> grad,
                       Eigen::VectorXd& r,          // scratch n (preallocated)
                       Eigen::VectorXd& grad_beta)   // scratch p  (preallocated)
{
    const int p = dat.p, n = dat.n;
    const double lse = par[p];
    const double lsb = par[p + 1];

    const double v_e = std::exp(2.0 * lse);
    const double v_b = std::exp(2.0 * lsb);

    if (!std::isfinite(v_e) || !std::isfinite(v_b) || v_e < 1e-300)
        return 1e300;

    r.noalias() = dat.y_s - dat.X_s * par.head(p);

    double d_lse = 0.0;
    double d_lsb = 0.0;

    grad_beta.setZero();

    double neg_ll = (n * 0.5) * LOG2PI;

    for (int gi = 0; gi < dat.G; ++gi) {
        const GroupInfo& g = dat.grps[gi];
        const int m = g.size;
        const int s = g.warm_start_params;
        const double a = v_e + m * v_b;

        double S = 0.0, Q = 0.0;
        for (int j = s; j < s + m; ++j) {
            S += r[j];
            Q += r[j] * r[j];
        }

        const double vb_over_a = v_b / a;   // ∈ [0, 1]: avoids v_e*a underflow below
        const double W      = Q - vb_over_a * S * S;
        const double inv_a  = 1.0 / a;
        const double gw     = g.w;  // case weight for this group

        // --- neg_ll contribution (weighted) ---
        neg_ll += gw * (0.5 * ((m - 1) * 2.0 * lse + std::log(a)) + W / (2.0 * v_e));

        // --- gradient for β (weighted) ---
        // Use (v_b/a)/v_e instead of v_b/(v_e*a) to avoid v_e*a underflowing to 0
        // which would produce Inf, and then Inf*0 = NaN when S≈0.
        const double coeff_S = gw * (vb_over_a / v_e) * S;
        for (int j = s; j < s + m; ++j) {
            const double w_j = gw * (-r[j] / v_e) + coeff_S;
            // grad_beta -= w_j * x_j  (remember: neg_ll, so ∂neg_ll/∂β = -score)
            grad_beta -= w_j * dat.X_s.row(j).transpose();
        }

        // --- gradient for log σ_e (weighted) ---
        d_lse += gw * ((m - 1) + v_e * inv_a + vb_over_a * (S * S) * inv_a / v_e - W / v_e);

        // --- gradient for log σ_b (weighted) ---
        d_lsb += gw * ((m * v_b * inv_a) - vb_over_a * (S * S) * inv_a);
    }

    grad.head(p) = -grad_beta;
    grad[p]     = d_lse;
    grad[p + 1] = d_lsb;

    return neg_ll;
}

// Scratch-allocating wrapper for non-hot callers (score export, fisher export).
double neg_ll_and_grad(const LMMData& dat,
                       const Eigen::Ref<const Eigen::VectorXd>& par,
                       Eigen::Ref<VectorXd> grad)
{
    Eigen::VectorXd r(dat.n), gb(dat.p);
    return neg_ll_and_grad(dat, par, grad, r, gb);
}

// Stack-allocated 2×2 matrix/vector — no heap alloc for group-level temporaries.
// MaxRows=2, MaxCols=2 tells Eigen to use internal fixed storage even for dynamic size.
using SM2 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 2, 2>;
using SV2 = Eigen::Matrix<double, Eigen::Dynamic, 1,              Eigen::ColMajor, 2, 1>;

// Analytical inverse of V = v_e*I + v_b*J (m×m, m ∈ {1,2}).
// For m=1: P = 1/(v_e+v_b).
// For m=2: det = v_e*(v_e+2*v_b), P_diag=(v_e+v_b)/det, P_off=-v_b/det.
inline SM2 lmm_P_analytic(int m, double v_e, double v_b) {
    SM2 P(m, m);
    if (m == 1) {
        P(0, 0) = 1.0 / (v_e + v_b);
    } else {
        const double det = v_e * (v_e + 2.0 * v_b);
        P(0,0) = P(1,1) = (v_e + v_b) / det;
        P(0,1) = P(1,0) = -v_b / det;
    }
    return P;
}

// lmm_analytic_hessian: uses preallocated scratch (buf_r: n, buf_XtA: p×max_m).
// All per-group temporaries are SM2/SV2 (stack-allocated ≤2×2) — no per-group heap alloc.
Eigen::MatrixXd lmm_analytic_hessian(const LMMData& dat,
                                     const Eigen::Ref<const Eigen::VectorXd>& par,
                                     Eigen::VectorXd& buf_r,    // scratch n
                                     Eigen::MatrixXd& buf_XtA)  // scratch p × max_m
{
    const int p = dat.p;
    const int k = (int)par.size();
    const double v_e = std::exp(2.0 * par[p]);
    const double v_b = std::exp(2.0 * par[p + 1]);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(k, k);
    if (!std::isfinite(v_e) || !std::isfinite(v_b) || v_e < 1e-300) {
        H.setConstant(NA_REAL);
        return H;
    }

    buf_r.noalias() = dat.y_s - dat.X_s * par.head(p);

    for (int gi = 0; gi < dat.G; ++gi) {
        const int m = dat.grps[gi].size;
        const int s = dat.grps[gi].warm_start_params;
        const double gw = dat.grps[gi].w;  // case weight for this group

        // All per-group matrices are SM2 (stack-allocated ≤2×2)
        const SM2 P = lmm_P_analytic(m, v_e, v_b);

        // Derivative matrices w.r.t. log σ_e and log σ_b
        SM2 dV_e(m, m); dV_e  = 2.0 * v_e * SM2::Identity(m, m);
        SM2 dV_b(m, m); dV_b.setConstant(2.0 * v_b);  // 2*v_b*J (J = all-ones)
        SM2 d2V_ee(m, m); d2V_ee = 4.0 * v_e * SM2::Identity(m, m);
        SM2 d2V_bb(m, m); d2V_bb.setConstant(4.0 * v_b);
        SM2 d2V_zero(m, m); d2V_zero.setZero();

        // dP/d(log σ) = -P * dV * P
        // dP_e = -2*v_e * P²
        SM2 tmp(m, m), dP_e(m, m), dP_b(m, m);
        tmp.noalias() = P * P;
        dP_e = (-2.0 * v_e) * tmp;

        // dP_b = -2*v_b * P * J * P  (P*J: row i = row_sum_of_P[i] broadcast across cols)
        SM2 PJ(m, m);
        for (int i = 0; i < m; ++i) PJ.row(i).setConstant(P.row(i).sum());
        tmp.noalias() = PJ * P;
        dP_b = (-2.0 * v_b) * tmp;

        // β-β block: gw * Xg^T * P * Xg  (buf_XtA reused, no Xg copy)
        buf_XtA.leftCols(m).noalias() = dat.X_s.middleRows(s, m).transpose() * P;
        H.topLeftCorner(p, p).noalias() += gw * (buf_XtA.leftCols(m) * dat.X_s.middleRows(s, m));

        // β-σ blocks: -gw * Xg^T * dP_a * rg  (SV2 on stack, no rg copy)
        SV2 sv(m);
        sv.noalias() = dP_e * buf_r.segment(s, m);
        H.topRightCorner(p, 1).noalias() -= gw * (dat.X_s.middleRows(s, m).transpose() * sv);
        sv.noalias() = dP_b * buf_r.segment(s, m);
        H.block(0, p + 1, p, 1).noalias() -= gw * (dat.X_s.middleRows(s, m).transpose() * sv);

        // σ-σ block: loop (a,b) ∈ {(e,e),(e,b),(b,b)}
        const SM2* dV_arr[2]  = {&dV_e,   &dV_b};
        const SM2* dP_arr[2]  = {&dP_e,   &dP_b};
        const SM2* d2V_arr[2][2] = {{&d2V_ee, &d2V_zero}, {&d2V_zero, &d2V_bb}};

        SM2 term(m, m), tmp2(m, m), tmp3(m, m);
        for (int a = 0; a < 2; ++a) {
            for (int b = a; b < 2; ++b) {
                tmp2.noalias() = *dP_arr[a] * *dV_arr[b];    // dP[a]*dV[b]
                tmp3.noalias() = P * *d2V_arr[a][b];          // P*d2V[a][b]
                const double h_tr = 0.5 * (tmp2 + tmp3).trace();

                // term = dP[a]*dV[b]*P + P*d2V[a][b]*P + P*dV[b]*dP[a]
                term.noalias()  = tmp2 * P;
                term.noalias() += tmp3 * P;
                tmp.noalias()   = P * *dV_arr[b];
                term.noalias() += tmp * *dP_arr[a];

                sv.noalias() = term * buf_r.segment(s, m);
                const double h_ab = h_tr - 0.5 * buf_r.segment(s, m).dot(sv);
                H(p + a, p + b) += gw * h_ab;
                if (a != b) H(p + b, p + a) += gw * h_ab;
            }
        }
    }

    H.bottomLeftCorner(2, p) = H.topRightCorner(p, 2).transpose();
    return (H + H.transpose()) / 2.0;
}

// ── Wrapper satisfying LBFGSpp operator() signature ─────────────────────────
class GaussianLMMObjective {
public:
    const LMMData& dat;

private:
    static int max_grp_size(const LMMData& d) {
        int mm = 1;
        for (const auto& g : d.grps) mm = std::max(mm, g.size);
        return mm;
    }

    Eigen::VectorXd m_r;        // scratch n  — reused across operator() and hessian()
    Eigen::VectorXd m_grad_beta; // scratch p  — reused across operator() calls
    Eigen::MatrixXd m_buf_XtA;  // scratch p × max_m — reused in hessian()

public:
    explicit GaussianLMMObjective(const LMMData& d)
        : dat(d),
          m_r(d.n),
          m_grad_beta(d.p),
          m_buf_XtA(d.p, max_grp_size(d)) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        if ((int)grad.size() != dat.p + 2) grad.resize(dat.p + 2);
        double nll = neg_ll_and_grad(dat, par, grad, m_r, m_grad_beta);
        if (!std::isfinite(nll) || nll >= 1e299) {
            grad.setZero();
            return 1e300;
        }
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        return lmm_analytic_hessian(dat, par, m_r, m_buf_XtA);
    }
};

// ── Hessian of neg_ll (for Fisher info / vcov) ─────────────────────────────
Eigen::MatrixXd lmm_fisher_hessian(const LMMData& dat,
                                   const Eigen::Ref<const Eigen::VectorXd>& par,
                                   double h_rel = 1e-4)
{
    (void)h_rel;
    // Non-hot path: allocate scratch locally
    Eigen::VectorXd buf_r(dat.n);
    int mm = 1;
    for (const auto& g : dat.grps) mm = std::max(mm, g.size);
    Eigen::MatrixXd buf_XtA(dat.p, mm);
    return lmm_analytic_hessian(dat, par, buf_r, buf_XtA);
}

// ── Starting values: OLS β, residual-based σ_e, σ_b = σ_e/2 ─────────────────
Eigen::VectorXd make_start(const LMMData& dat)
{
    const int p = dat.p;
    // OLS: β = (X^T X)^{-1} X^T y
    Eigen::MatrixXd XtX = dat.X_s.transpose() * dat.X_s;
    Eigen::VectorXd Xty = dat.X_s.transpose() * dat.y_s;
    Eigen::VectorXd beta = XtX.ldlt().solve(Xty);

    // Residual SD
    Eigen::VectorXd res = dat.y_s - dat.X_s * beta;
    double sigma2_e = res.squaredNorm() / std::max(dat.n - p, 1);
    double sigma_e  = std::sqrt(std::max(sigma2_e, 1e-8));

    Eigen::VectorXd warm_start_params(p + 2);
    warm_start_params.head(p)   = beta;
    warm_start_params[p]     = std::log(sigma_e);        // log σ_e
    warm_start_params[p + 1] = std::log(sigma_e * 0.5);  // log σ_b ≈ σ_e/2
    return warm_start_params;
}

} // anonymous namespace

// ── R-exported: fit Gaussian LMM ─────────────────────────────────────────────
// [[Rcpp::export]]
List fast_gaussian_lmm_cpp(
    SEXP X_r,       // n × p, intercept in col 0, treatment in col 1
    SEXP y_r,
    SEXP group_id_r, // 1-based group IDs (length n)
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
    bool  estimate_only = false,
    int   maxit  = 300,
    double eps_g = 1e-6,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "lbfgs",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue  // per-obs case weights (NULL = uniform)
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());

    const int n = y.size(), p = X.cols();

    // Convert group_id to 0-based sorted integer IDs
    std::vector<int> gid(n);
    {
        // Map R group ids (any positive integers) to 0-based consecutive ints
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }

    std::vector<double> w_vec;
    if (weights.isNotNull()) {
        NumericVector wts(weights);
        w_vec.assign(wts.begin(), wts.end());
    }
    LMMData dat(y, X, gid, w_vec);
    GaussianLMMObjective obj(dat);

    // Starting point
    Eigen::VectorXd par = make_start(dat);
    if (warm_start_params.isNotNull()) {
        NumericVector sp(warm_start_params);
        if (sp.size() == p + 2)
            for (int i = 0; i < p + 2; ++i) par[i] = sp[i];
    } else if (warm_start_beta.isNotNull()) {
        VectorXd sb = as<VectorXd>(warm_start_beta);
        if (sb.size() == p + 2) {
            par = sb;
        } else if (sb.size() == p) {
            par.head(p) = sb;
        }
    }
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 2, fixed_idx, fixed_values);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    double neg_ll = 1e300;
    int niter = maxit;
    bool converged = false;
    try {
        LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
        par = fit.params;
        neg_ll = fit.value;
        niter = fit.niter;
        converged = std::isfinite(neg_ll) && fit.converged;
    } catch (...) {
        converged = false;
    }

    // Return par names: β[0..p-1], log_sigma_e, log_sigma_b
    NumericVector b_r(p + 2);
    for (int i = 0; i < p + 2; ++i) b_r[i] = par[i];
    CharacterVector b_names(p + 2);
    for (int i = 0; i < p; ++i) b_names[i] = "b" + std::to_string(i);
    b_names[p]   = "log_sigma_e";
    b_names[p+1] = "log_sigma_b";
    b_r.names() = b_names;

    if (estimate_only) {
        return List::create(
            Named("b")         = b_r,
            Named("ssq_b_T")   = NA_REAL,
            Named("neg_loglik")= neg_ll,
            Named("converged") = converged,
            Named("niter")     = niter
        );
    }

    // Variance-covariance via Hessian of neg_ll (reuses obj's scratch buffers)
    Eigen::MatrixXd H = obj.hessian(par);
    Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);

    double ssq_b_T = NA_REAL;
    NumericMatrix vcov_r(p + 2, p + 2);
    std::fill(vcov_r.begin(), vcov_r.end(), NA_REAL);

    if (ldlt.info() == Eigen::Success) {
        Eigen::MatrixXd V_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
        Eigen::MatrixXd V = expand_free_covariance(p + 2, fixed_spec, V_free, true);
        if (V.allFinite()) {
            ssq_b_T = V(1, 1);   // treatment is index 1 (after intercept at 0)
            for (int i = 0; i < p + 2; ++i)
                for (int j = 0; j < p + 2; ++j)
                    vcov_r(i, j) = V(i, j);
        }
    }

    return List::create(
        Named("b")         = b_r,
        Named("params")    = b_r,
        Named("ssq_b_T")   = ssq_b_T,
        Named("vcov")      = vcov_r,
        Named("neg_loglik")= neg_ll,
        Named("converged") = converged,
        Named("niter")     = niter,
        Named("fisher_information") = H
    );
}

// ── R-exported: GLS estimator with fixed variance components ─────────────────
// Given fixed log_sigma_e and log_sigma_b (from a prior full fit), solve for
// beta via GLS without any L-BFGS optimization.  ~100x faster than full MLE.
// Valid for permutation tests (VC fixed at null-fit MLE — exact test by
// exchangeability) and for non-studentised bootstrap (VC approximated).
//
// GLS normal equations with V_g = v_e*I + v_b*J:
//   A     = sum_g gw * [ X_g'X_g - c_g * sx_g sx_g' ]   (v_e cancels)
//   b_rhs = sum_g gw * [ X_g'y_g - c_g * sx_g sy_g  ]
// where a_g = v_e + m_g*v_b,  c_g = v_b/(v_e*a_g),
//       sx_g = X_g'1_m,  sy_g = 1_m'y_g.
// [[Rcpp::export]]
Rcpp::NumericVector fast_gaussian_lmm_gls_cpp(
    SEXP X_r,
    SEXP y_r,
    SEXP group_id_r,
    double log_sigma_e,
    double log_sigma_b,
    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());

    const int n = y.size(), p = X.cols();

    std::vector<int> gid(n);
    {
        const int* gid_ptr = group_id_int.begin();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }

    std::vector<double> w_vec;
    if (weights.isNotNull()) {
        NumericVector wts(weights);
        w_vec.assign(wts.begin(), wts.end());
    }
    LMMData dat(y, X, gid, w_vec);

    const double v_e = std::exp(2.0 * log_sigma_e);
    const double v_b = std::exp(2.0 * log_sigma_b);
    if (!std::isfinite(v_e) || !std::isfinite(v_b) || v_e < 1e-300)
        return Rcpp::NumericVector(p, NA_REAL);

    Eigen::MatrixXd A     = Eigen::MatrixXd::Zero(p, p);
    Eigen::VectorXd b_rhs = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd sx(p);

    for (int gi = 0; gi < dat.G; ++gi) {
        const GroupInfo& g = dat.grps[gi];
        const int m = g.size, s = g.warm_start_params;
        const double gw = g.w;
        const double a_g = v_e + m * v_b;
        const double c_g = v_b / (v_e * a_g);

        sx.noalias() = dat.X_s.middleRows(s, m).colwise().sum().transpose();
        const double sy = dat.y_s.segment(s, m).sum();

        A.noalias()     += gw * (dat.X_s.middleRows(s, m).transpose() * dat.X_s.middleRows(s, m)
                                 - c_g * (sx * sx.transpose()));
        b_rhs.noalias() += gw * (dat.X_s.middleRows(s, m).transpose() * dat.y_s.segment(s, m)
                                 - c_g * (sx * sy));
    }

    Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
    if (ldlt.info() != Eigen::Success)
        return Rcpp::NumericVector(p, NA_REAL);
    Eigen::VectorXd beta = ldlt.solve(b_rhs);

    NumericVector result(p);
    CharacterVector names(p);
    for (int i = 0; i < p; ++i) {
        if (!std::isfinite(beta[i])) return Rcpp::NumericVector(p, NA_REAL);
        result[i] = beta[i];
        names[i]  = "b" + std::to_string(i);
    }
    result.names() = names;
    return result;
}

// ── R-exported: score (gradient of log_lik) at arbitrary par ─────────────────
// [[Rcpp::export]]
NumericVector get_gaussian_lmm_score_cpp(
    SEXP X_r,
    SEXP y_r,
    SEXP group_id_r,
    SEXP par_sexp
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    NumericVector par_r(par_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
    Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

    const int n = y.size();
    std::vector<int> gid(n);
    {
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    const int k = X.cols() + 2;  // p betas + log_sigma_e + log_sigma_b
    if (par.size() != k) {
        return NumericVector(k, NA_REAL);
    }
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(k);
    neg_ll_and_grad(dat, par, grad);
    // score = -grad of neg_ll
    Eigen::VectorXd score = -grad;
    // Safety net: if any component is non-finite, return all-NA.
    for (int i = 0; i < k; ++i) {
        if (!std::isfinite(score[i])) {
            return NumericVector(k, NA_REAL);
        }
    }
    return wrap(score);
}

// ── R-exported: observed Fisher information (Hessian of neg_ll) at par ───────
// [[Rcpp::export]]
NumericMatrix get_gaussian_lmm_fisher_cpp(
    SEXP X_r,
    SEXP y_r,
    SEXP group_id_r,
    SEXP par_sexp,
    double h_rel = 1e-4
) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    IntegerVector group_id_int(group_id_r);
    NumericVector par_r(par_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
    Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

    const int n = y.size();
    std::vector<int> gid(n);
    {
        const int* gid_ptr = group_id.data();
        std::vector<int> gid_r(gid_ptr, gid_ptr + n);
        std::vector<int> uniq = gid_r;
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        for (int i = 0; i < n; ++i)
            gid[i] = (int)(std::lower_bound(uniq.begin(), uniq.end(), gid_r[i]) - uniq.begin());
    }
    LMMData dat(y, X, gid);
    Eigen::MatrixXd H = lmm_fisher_hessian(dat, par, h_rel);
    return wrap(H);
}
