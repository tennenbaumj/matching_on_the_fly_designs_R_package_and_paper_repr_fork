#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

namespace {

inline double log_sum_exp_wf(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline Eigen::ArrayXd clamp_weibull_wf(const Eigen::ArrayXd& w) {
	return w.min(700.0);
}

struct GHRuleWF {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleWF gauss_hermite_rule_wf(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleWF rule;
	rule.nodes = es.eigenvalues();
	Eigen::VectorXd weights = std::sqrt(M_PI) * es.eigenvectors().row(0).array().square().matrix();
	rule.log_norm_weights = weights.array().log() - 0.5 * std::log(M_PI);
	return rule;
}

// Weibull AFT GLMM with log-Normal random intercept per cluster.
// params: [beta(p), log_sigma_eps(1), log_sigma_u(1)]
// log_sigma_eps: log of Weibull scale (error SD on log scale)
// log_sigma_u: log of random-intercept SD
class WeibullFrailtyLikelihood {
private:
	Eigen::VectorXd m_y;
	Eigen::VectorXd m_dead;
	Eigen::MatrixXd m_X;
	Eigen::VectorXd m_log_y;
	const int m_n;
	const int m_p;
	const GHRuleWF m_gh;
	const double m_max_abs_log_sigma;

	std::vector<int> m_group_start;
	std::vector<int> m_group_end;
	Eigen::VectorXd m_group_dead_sum;
	Eigen::VectorXd m_group_dead_log_y_sum;
	int m_n_groups;
	int m_max_group_size;
	Eigen::VectorXd m_eta_all;
	Eigen::VectorXd m_u_vals;
	Eigen::MatrixXd m_log_terms_mat;
	Eigen::MatrixXd m_r_all_k_mat;
	Eigen::MatrixXd m_w_all_k_mat;
	Eigen::VectorXd m_ll_g_vec;
	Eigen::VectorXd m_weighted_r;
	Eigen::VectorXd m_base_wi;    // (log_y - eta) / sigma_eps per obs — precomputed once per optimizer step
	Eigen::VectorXd m_base_exp;   // exp(m_base_wi) per obs — n exp calls, contiguous for SIMD
	Eigen::MatrixXd m_post_weights; // G x K posterior weights — preallocated to avoid per-call heap alloc

	void build_group_structure(const Eigen::Ref<const Eigen::VectorXi>& group_id) {
		std::vector<int> order(m_n);
		std::iota(order.begin(), order.end(), 0);
		std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
			return group_id[a] < group_id[b];
		});

		Eigen::VectorXd y_s(m_n), dead_s(m_n), log_y_s(m_n);
		Eigen::MatrixXd X_s(m_n, m_p);
		for (int i = 0; i < m_n; ++i) {
			y_s[i]     = m_y[order[i]];
			dead_s[i]  = m_dead[order[i]];
			log_y_s[i] = m_log_y[order[i]];
			X_s.row(i) = m_X.row(order[i]);
		}
		m_y     = y_s;
		m_dead  = dead_s;
		m_log_y = log_y_s;
		m_X     = X_s;

			Eigen::VectorXi sorted_gid(m_n);
			for (int i = 0; i < m_n; ++i) sorted_gid[i] = group_id[order[i]];

			auto layout = build_contiguous_group_layout(m_n, [&](int i) { return sorted_gid[i]; });
			m_group_start = layout.start;
			m_group_end.resize(layout.G);
			for (int gi = 0; gi < layout.G; ++gi) {
				m_group_end[gi] = layout.start[gi] + layout.size[gi];
			}
		m_n_groups      = layout.G;
		m_max_group_size = layout.max_size;

		m_group_dead_sum.resize(m_n_groups);
		m_group_dead_log_y_sum.resize(m_n_groups);
		for (int g = 0; g < m_n_groups; ++g) {
			m_group_dead_sum[g]       = m_dead.segment(m_group_start[g], m_group_end[g] - m_group_start[g]).sum();
			m_group_dead_log_y_sum[g] = (m_dead.segment(m_group_start[g], m_group_end[g] - m_group_start[g]).array() *
			                            m_log_y.segment(m_group_start[g], m_group_end[g] - m_group_start[g]).array()).sum();
		}
	}

public:
	WeibullFrailtyLikelihood(const Eigen::Ref<const Eigen::VectorXd>& y, 
							 const Eigen::Ref<const Eigen::VectorXd>& dead, 
							 const Eigen::Ref<const Eigen::MatrixXd>& X, 
							 const Eigen::Ref<const Eigen::VectorXi>& group_id, 
							 int n_gh, 
							 double max_abs_log_sigma) :
		m_y(y), m_dead(dead), m_X(X), m_log_y(y.array().log().matrix()),
		m_n(y.size()), m_p(X.cols()), m_gh(gauss_hermite_rule_wf(n_gh)),
		m_max_abs_log_sigma(max_abs_log_sigma) {
		
		build_group_structure(group_id);

		m_eta_all.resize(m_n);
		m_u_vals = std::sqrt(2.0) * m_gh.nodes;
		m_log_terms_mat.resize(m_n_groups, m_gh.nodes.size());
		m_r_all_k_mat.resize(m_n, m_gh.nodes.size());
		m_w_all_k_mat.resize(m_n, m_gh.nodes.size());
		m_ll_g_vec.resize(m_n_groups);
		m_weighted_r.resize(m_n);
		m_base_wi.resize(m_n);
		m_base_exp.resize(m_n);
		m_post_weights.resize(m_n_groups, m_gh.nodes.size());
	}

	double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
		const Eigen::VectorXd beta = params.head(m_p);
		const double log_sigma_eps  = std::max(-m_max_abs_log_sigma, std::min(m_max_abs_log_sigma, params[m_p]));
		const double sigma_eps      = std::exp(log_sigma_eps);
		const double log_sigma_u    = std::max(-m_max_abs_log_sigma, std::min(m_max_abs_log_sigma, params[m_p + 1]));
		const double sigma_u        = std::exp(log_sigma_u);
		const double inv_eps        = 1.0 / sigma_eps;

		m_eta_all.noalias() = m_X * beta;
		const int K = (int)m_gh.nodes.size();

		// Precompute per-obs base_wi = (log_y - eta) / sigma_eps and base_exp = exp(base_wi).
		// n exp calls over a contiguous array — vectorized by the compiler into batched SIMD exp.
		m_base_wi  = (m_log_y - m_eta_all) * inv_eps;
		m_base_exp = m_base_wi.array().exp();

		for (int k = 0; k < K; ++k) {
			const double uk          = sigma_u * m_u_vals[k];
			const double delta_k     = -uk * inv_eps;          // node-only shift: wik = base_wi[i] + delta_k
			const double exp_delta_k = std::exp(delta_k);      // one exp per node (K=20 total)

			for (int g = 0; g < m_n_groups; ++g) {
				const int start = m_group_start[g];
				const int sz    = m_group_end[g] - start;
				// Per-obs conditional log-lik: dead*(w - log_sigma_eps - log_y) - exp(w)
				double log_lik_g_k = m_group_dead_sum[g] * (-log_sigma_eps) - m_group_dead_log_y_sum[g];

				for (int i = 0; i < sz; ++i) {
					const int idx = start + i;
					const double wik = m_base_wi[idx] + delta_k;
					// exp(wik) = exp(base_wi) * exp(delta_k); fallback for rare float overflow.
					double exp_wik = m_base_exp[idx] * exp_delta_k;
					if (!std::isfinite(exp_wik)) exp_wik = std::exp(std::min(wik, 700.0));
					log_lik_g_k += m_dead[idx] * wik - exp_wik;
					m_r_all_k_mat(idx, k) = exp_wik;
					m_w_all_k_mat(idx, k) = wik;
				}
				m_log_terms_mat(g, k) = log_lik_g_k + m_gh.log_norm_weights[k];
			}
		}

		// Fused LSE + posterior weights: compute exp(row - max) once, reuse as weights.
		// Saves G*K=1000 exp calls vs. separate log_sum_exp + post_weights passes.
		double total_neg_ll = 0.0;
		grad.setZero();
		for (int g = 0; g < m_n_groups; ++g) {
			const auto row = m_log_terms_mat.row(g);
			const double mx = row.maxCoeff();
			m_post_weights.row(g) = (row.array() - mx).exp();
			const double S = m_post_weights.row(g).sum();
			m_post_weights.row(g) /= S;
			const double ll_g = mx + std::log(S);
			total_neg_ll -= ll_g;
		}

		// grad accumulates the gradient of the NEGATIVE log-likelihood
		// (matching the convention of the other likelihood functors, e.g.
		// fast_weibull_regression.cpp: return -loglik, grad = d(-loglik)/d params).
		for (int k = 0; k < K; ++k) {
			const double uk = sigma_u * m_u_vals[k];
			for (int g = 0; g < m_n_groups; ++g) {
				const double pw_gk = m_post_weights(g, k);
				const int start    = m_group_start[g];
				const int sz       = m_group_end[g] - start;
				double sum_s_gk = 0.0;
				for (int i = 0; i < sz; ++i) {
					const int idx = start + i;
					const double d_log_f_d_eta = (m_r_all_k_mat(idx, k) - m_dead[idx]) * inv_eps;
					sum_s_gk += d_log_f_d_eta;
					grad.head(m_p).noalias() -= pw_gk * d_log_f_d_eta * m_X.row(idx).transpose();

					const double wik = m_w_all_k_mat(idx, k);
					const double d_log_f_d_log_sigma_eps = m_r_all_k_mat(idx, k) * wik - m_dead[idx] * (wik + 1.0);
					grad[m_p] -= pw_gk * d_log_f_d_log_sigma_eps;
				}
				// d loglik_g / d log_sigma_u flows only through u_k = sigma_u*sqrt(2)*t_k:
				// d w_ik / d log_sigma_u = -u_k/sigma_eps, so d l_ik/d log_sigma_u = s_ik*u_k.
				grad[m_p + 1] -= pw_gk * uk * sum_s_gk;
			}
		}

		return total_neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
		return numerical_hessian(*this, params);
	}
};

} // namespace

// [[Rcpp::export]]
double get_weibull_frailty_neg_loglik_cpp(
	SEXP X_sexp,
	SEXP y_sexp,
	SEXP dead_sexp,
	SEXP group_id_sexp,
	SEXP params_sexp,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	NumericMatrix X_r(X_sexp);
	NumericVector y_r(y_sexp);
	NumericVector dead_r(dead_sexp);
	IntegerVector group_id_r(group_id_sexp);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	Eigen::VectorXd grad(params.size());
	return obj(params, grad);
}

// [[Rcpp::export]]
Eigen::VectorXd get_weibull_frailty_score_cpp(
	SEXP X_sexp,
	SEXP y_sexp,
	SEXP dead_sexp,
	SEXP group_id_sexp,
	SEXP params_sexp,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	NumericMatrix X_r(X_sexp);
	NumericVector y_r(y_sexp);
	NumericVector dead_r(dead_sexp);
	IntegerVector group_id_r(group_id_sexp);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_weibull_frailty_hessian_cpp(
	SEXP X_sexp,
	SEXP y_sexp,
	SEXP dead_sexp,
	SEXP group_id_sexp,
	SEXP params_sexp,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	NumericMatrix X_r(X_sexp);
	NumericVector y_r(y_sexp);
	NumericVector dead_r(dead_sexp);
	IntegerVector group_id_r(group_id_sexp);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	return -obj.hessian(params);
}

// [[Rcpp::export]]
List fast_weibull_frailty_cpp(
	SEXP X_sexp,
	SEXP y_sexp,
	SEXP dead_sexp,
	SEXP group_id_sexp,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	NumericMatrix X_r(X_sexp);
	NumericVector y_r(y_sexp);
	NumericVector dead_r(dead_sexp);
	IntegerVector group_id_r(group_id_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());

	const int p     = X.cols();
	const int n_par = p + 2;

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);

	Eigen::VectorXd par(n_par);
	if (warm_start_params.isNotNull()) {
		par = Rcpp::as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_params));
	} else if (warm_start_beta.isNotNull()) {
		VectorXd sb = as<VectorXd>(warm_start_beta);
		if (sb.size() == n_par) {
			par = sb;
		} else if (sb.size() == p) {
			par.head(p) = sb;
			par[p] = 0.0;
			par[p+1] = -3.0;
		}
	} else {
		Eigen::VectorXd log_y = y.array().log().matrix();
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
		Eigen::VectorXd beta_init = qr.solve(log_y);
		par.head(p) = beta_init;
		Eigen::VectorXd resid = log_y - X * par.head(p);
		const double std_resid = std::sqrt(resid.squaredNorm() / std::max(1, (int)y.size() - p));
		par[p]     = std::log(std_resid * 0.7797);
		par[p + 1] = 0.0;
	}

	FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll  = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par       = fit.params;
		neg_ll    = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {}
	// If neg_ll is still NA (optimizer threw or returned non-finite value), evaluate at last iterate
	if (!std::isfinite(neg_ll) && par.allFinite()) {
		Eigen::VectorXd grad_tmp(n_par);
		double val = obj(par, grad_tmp);
		if (std::isfinite(val)) neg_ll = val;
	}

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(n_par, n_par, NA_REAL);
	Eigen::VectorXd score = Eigen::VectorXd::Constant(n_par, NA_REAL);
	Eigen::MatrixXd information = Eigen::MatrixXd::Constant(n_par, n_par, NA_REAL);
	if (!estimate_only && converged) {
		Eigen::VectorXd grad(n_par);
		obj(par, grad);
		score = -grad;
		information = obj.hessian(par);
		Eigen::MatrixXd info_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(info_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(info_free.rows(), info_free.cols()));
			vcov = expand_free_covariance(n_par, fixed_spec, inv_free, true);
			if (vcov.allFinite()) ssq_b_T = vcov(0, 0);  // j_T = 0 usually
		}
	}

	return List::create(
		Named("params")        = par,
		Named("b")             = par.head(p),
		Named("log_sigma_eps") = par[p],
		Named("log_sigma_u")   = par[p + 1],
		Named("ssq_b_T")       = ssq_b_T,
		Named("vcov")          = vcov,
		Named("score")         = score,
		Named("observed_information") = information,
		Named("information")   = information,
		Named("information_type") = "observed",
		Named("hessian")       = -information,
		Named("converged")     = converged,
		Named("neg_loglik")    = neg_ll,
		Named("neg_ll")        = neg_ll,
		Named("loglik")        = R_finite(neg_ll) ? -neg_ll : NA_REAL
	);
}
