// Poisson GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  log E[Y_ij | u_i] = X_ij' beta + u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in {0, 1, 2, ...} count responses
//   X includes an intercept column
//
// Parameter vector: par = [beta_0, beta_1(treatment), ..., beta_{p-1}, log_sigma]
//   Total length: p + 1

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

namespace {

struct GHRuleP {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleP gauss_hermite_rule_poisson(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleP rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
		                    - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_p(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline double soft_barrier(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double soft_barrier_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double soft_barrier_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

inline Eigen::ArrayXd clamp_eta_pg(const Eigen::ArrayXd& eta) {
	return eta.min(700.0);
}

struct PoissonGLMMData {
	Eigen::MatrixXd X_s;
	Eigen::VectorXd y_s;
	Eigen::VectorXd log_fact_y;  // precomputed log(y!) = lgamma(y+1)
	Eigen::VectorXd w_s;         // per-row weights (1.0 if unweighted)
	std::vector<int> grp_start;
	std::vector<int> grp_size;
	int n, p, G;
	GHRuleP gh;

	PoissonGLMMData(
		const Eigen::Ref<const Eigen::MatrixXd>& X,
		const Eigen::Ref<const Eigen::VectorXd>& y,
		const std::vector<int>& group_id,
		int n_gh,
		const Eigen::VectorXd* row_weights = nullptr
	) : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule_poisson(n_gh)) {

		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id[a] < group_id[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		log_fact_y.resize(n);
		w_s.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y[ord[i]];
			log_fact_y[i] = std::lgamma(y_s[i] + 1.0);
			w_s[i] = row_weights ? (*row_weights)[ord[i]] : 1.0;
		}

		int prev = -1;
		for (int i = 0; i < n; ++i) {
			int g = group_id[ord[i]];
			if (g != prev) {
				grp_start.push_back(i);
				grp_size.push_back(1);
				prev = g;
			} else {
				grp_size.back()++;
			}
		}
		G = (int)grp_start.size();
	}
};

class PoissonGLMMObjective {
	const PoissonGLMMData& dat;
	const int m_n_nodes;
	const int m_total;    // p + 1
	const int m_max_grp;  // max group size for per-group scratch buffers

	// Preallocated once at construction — written each call, never heap-allocated per call.
	// Shared by operator() and hessian():
	Eigen::VectorXd m_eta_all;        // n
	Eigen::ArrayXd  m_term_all_k;     // n  (per-node log-likelihood terms)
	Eigen::MatrixXd m_log_terms_mat;  // G × n_nodes
	Eigen::MatrixXd m_mu_all_k_mat;   // n × n_nodes  (replaces std::vector<VectorXd>)
	Eigen::VectorXd m_ll_g_vec;       // G
	// operator()-only:
	Eigen::VectorXd m_post_k_exp;     // n  (posterior weight expanded to subjects)
	Eigen::VectorXd m_wres_all_k;     // n  (weighted residuals for beta gradient)
	// hessian()-only:
	Eigen::MatrixXd m_E_Hik_sum;      // total × total
	Eigen::MatrixXd m_E_GiGiT_sum;    // total × total
	Eigen::MatrixXd m_G_avg_outer;    // total × total
	Eigen::VectorXd m_pk_vec;         // G
	Eigen::VectorXd m_pk_exp;         // n  (pk expanded per subject)
	Eigen::VectorXd m_G_avg_gi;       // total
	Eigen::MatrixXd m_E_GiGiT_gi;     // total × total
	Eigen::VectorXd m_G_ik;           // total
	Eigen::VectorXd m_wres_k_gi;      // max_grp (per-group scratch)
	Eigen::VectorXd m_wmu_k_gi;       // max_grp (per-group scratch)

public:
	explicit PoissonGLMMObjective(const PoissonGLMMData& d)
		: dat(d),
		  m_n_nodes((int)d.gh.nodes.size()),
		  m_total(d.p + 1),
		  m_max_grp(d.grp_size.empty() ? 1 : *std::max_element(d.grp_size.begin(), d.grp_size.end())),
		  m_eta_all(d.n),
		  m_term_all_k(d.n),
		  m_log_terms_mat(d.G, (int)d.gh.nodes.size()),
		  m_mu_all_k_mat(d.n, (int)d.gh.nodes.size()),
		  m_ll_g_vec(d.G),
		  m_post_k_exp(d.n),
		  m_wres_all_k(d.n),
		  m_E_Hik_sum(d.p + 1, d.p + 1),
		  m_E_GiGiT_sum(d.p + 1, d.p + 1),
		  m_G_avg_outer(d.p + 1, d.p + 1),
		  m_pk_vec(d.G),
		  m_pk_exp(d.n),
		  m_G_avg_gi(d.p + 1),
		  m_E_GiGiT_gi(d.p + 1, d.p + 1),
		  m_G_ik(d.p + 1),
		  m_wres_k_gi(d.grp_size.empty() ? 1 : *std::max_element(d.grp_size.begin(), d.grp_size.end())),
		  m_wmu_k_gi(d.grp_size.empty() ? 1 : *std::max_element(d.grp_size.begin(), d.grp_size.end())) {}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& par, Eigen::VectorXd& grad) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const Eigen::ArrayXd  y_all  = dat.y_s.array();

		m_eta_all.noalias() = dat.X_s * par.head(dat.p);

		for (int k = 0; k < m_n_nodes; ++k) {
			m_mu_all_k_mat.col(k) = (m_eta_all.array() + b_vals[k]).min(700.0).exp().matrix();
			m_term_all_k = dat.w_s.array() * (y_all * (m_eta_all.array() + b_vals[k])
			               - m_mu_all_k_mat.col(k).array() - dat.log_fact_y.array());
			for (int gi = 0; gi < dat.G; ++gi) {
				m_log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] +
				    m_term_all_k.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}

		double total_nll = soft_barrier(log_sigma);
		for (int gi = 0; gi < dat.G; ++gi) {
			m_ll_g_vec[gi] = log_sum_exp_p(m_log_terms_mat.row(gi));
			if (!std::isfinite(m_ll_g_vec[gi])) { grad.setZero(m_total); return 1e100; }
			total_nll -= m_ll_g_vec[gi];
		}

		grad.setZero(m_total);
		const double d_pen = std::abs(log_sigma) - 5.0;
		if (d_pen > 0.0) grad[dat.p] += 20.0 * d_pen * (log_sigma > 0 ? 1.0 : -1.0);

		Eigen::VectorXd grad_beta = Eigen::VectorXd::Zero(dat.p);
		double grad_log_sigma = 0.0;

		for (int k = 0; k < m_n_nodes; ++k) {
			double dLL_dlog_sigma_k = 0.0;

			for (int gi = 0; gi < dat.G; ++gi) {
				const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
				const int s  = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				if (pk < 1e-15) {
					m_post_k_exp.segment(s, sz).setZero();
					continue;
				}
				m_post_k_exp.segment(s, sz).setConstant(pk);
				const double wres_sum_k_gi = (dat.w_s.array().segment(s, sz) *
				    (y_all.segment(s, sz) - m_mu_all_k_mat.col(k).array().segment(s, sz))).sum();
				dLL_dlog_sigma_k += pk * wres_sum_k_gi * b_vals[k];
			}

			m_wres_all_k = dat.w_s.cwiseProduct(y_all.matrix() - m_mu_all_k_mat.col(k));
			grad_beta.noalias() -= dat.X_s.transpose() * m_post_k_exp.cwiseProduct(m_wres_all_k);
			grad_log_sigma -= dLL_dlog_sigma_k;
		}

		grad.head(dat.p) = grad_beta;
		grad[dat.p] += grad_log_sigma;
		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& par) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const Eigen::ArrayXd  y_all  = dat.y_s.array();

		m_eta_all.noalias() = dat.X_s * par.head(dat.p);

		for (int k = 0; k < m_n_nodes; ++k) {
			m_mu_all_k_mat.col(k) = (m_eta_all.array() + b_vals[k]).min(700.0).exp().matrix();
			m_term_all_k = dat.w_s.array() * (y_all * (m_eta_all.array() + b_vals[k])
			               - m_mu_all_k_mat.col(k).array() - dat.log_fact_y.array());
			for (int gi = 0; gi < dat.G; ++gi) {
				m_log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] +
				    m_term_all_k.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}
		for (int gi = 0; gi < dat.G; ++gi) m_ll_g_vec[gi] = log_sum_exp_p(m_log_terms_mat.row(gi));

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m_total, m_total);
		H(dat.p, dat.p) = soft_barrier_hessian(log_sigma);

		m_E_Hik_sum.setZero();
		m_E_GiGiT_sum.setZero();
		m_G_avg_outer.setZero();

		for (int k = 0; k < m_n_nodes; k++) {
			for (int gi = 0; gi < dat.G; gi++) m_pk_vec[gi] = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
			for (int gi = 0; gi < dat.G; gi++)
				m_pk_exp.segment(dat.grp_start[gi], dat.grp_size[gi]).setConstant(m_pk_vec[gi]);

			// Beta-Beta block: -X^T * diag(pk_exp * w * mu_k) * X
			m_E_Hik_sum.topLeftCorner(dat.p, dat.p).noalias() -=
				weighted_crossprod(dat.X_s, m_pk_exp.cwiseProduct(dat.w_s).cwiseProduct(m_mu_all_k_mat.col(k)));

			const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
			for (int gi = 0; gi < dat.G; ++gi) {
				const double pk = m_pk_vec[gi];
				if (pk < 1e-15) continue;
				const int s  = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const double sum_wmu  = (dat.w_s.array().segment(s, sz) * m_mu_all_k_mat.col(k).array().segment(s, sz)).sum();
				const double sum_wres = (dat.w_s.array().segment(s, sz) *
				    (y_all.segment(s, sz) - m_mu_all_k_mat.col(k).array().segment(s, sz))).sum();
				m_E_Hik_sum(dat.p, dat.p) += pk * ((-sum_wmu * node_factor * node_factor * sigma + sum_wres * node_factor) * sigma);

				// Beta-Sigma block: -X_g^T * (w_g * mu_k_g) * node_factor * sigma
				m_wmu_k_gi.head(sz) = dat.w_s.segment(s, sz).cwiseProduct(m_mu_all_k_mat.col(k).segment(s, sz));
				m_E_Hik_sum.block(0, dat.p, dat.p, 1).noalias() +=
					pk * (-(dat.X_s.middleRows(s, sz).transpose() * m_wmu_k_gi.head(sz)) * (node_factor * sigma));
			}
		}

		for (int gi = 0; gi < dat.G; ++gi) {
			m_G_avg_gi.setZero();
			m_E_GiGiT_gi.setZero();
			const int s  = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];
			for (int k = 0; k < m_n_nodes; k++) {
				const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
				if (pk < 1e-15) continue;
				m_wres_k_gi.head(sz) = (dat.w_s.array().segment(s, sz) *
				    (y_all.segment(s, sz) - m_mu_all_k_mat.col(k).array().segment(s, sz))).matrix();
				m_G_ik.setZero();
				m_G_ik.head(dat.p).noalias() = dat.X_s.middleRows(s, sz).transpose() * m_wres_k_gi.head(sz);
				const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
				m_G_ik[dat.p] = m_wres_k_gi.head(sz).sum() * node_factor * sigma;
				m_G_avg_gi.noalias() += pk * m_G_ik;
				m_E_GiGiT_gi.noalias() += pk * (m_G_ik * m_G_ik.transpose());
			}
			m_E_GiGiT_sum.noalias() += m_E_GiGiT_gi;
			m_G_avg_outer.noalias()  += m_G_avg_gi * m_G_avg_gi.transpose();
		}

		for (int r = 0; r < dat.p; r++) for (int c = 0; c < r; c++) m_E_Hik_sum(r, c) = m_E_Hik_sum(c, r);
		m_E_Hik_sum.block(dat.p, 0, 1, dat.p) = m_E_Hik_sum.block(0, dat.p, dat.p, 1).transpose();

		H -= (m_E_Hik_sum + m_E_GiGiT_sum - m_G_avg_outer);
		return H;
	}
};

} // namespace

// [[Rcpp::export]]
SEXP fast_poisson_glmm_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	int j_T,
	Nullable<NumericVector> warm_start_params = R_NilValue,
	bool smart_cold_start = true,
	bool estimate_only = false,
	int n_gh = 20,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> row_weights = R_NilValue
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());

	if (X_r.rows() != y_r.size() || X_r.rows() != group_id_r.size()) {
		Rcpp::stop("Dimension mismatch: X_r has %d rows, y_r has %d elements, group_id_r has %d elements",
		           X_r.rows(), y_r.size(), group_id_r.size());
	}

	const int n = X.rows();
	const int p = X.cols();
	const int total = p + 1;

	std::vector<int> gid_v(n);
	for (int i = 0; i < n; ++i) gid_v[i] = group_id[i];

	Eigen::VectorXd rw_vec;
	const Eigen::VectorXd* rw_ptr = nullptr;
	if (row_weights.isNotNull()) {
		NumericVector rw_r(row_weights);
		if (rw_r.size() != n)
			Rcpp::stop("row_weights length (%d) must equal nrow(X) (%d)", rw_r.size(), n);
		rw_vec = Eigen::Map<const Eigen::VectorXd>(rw_r.begin(), n);
		rw_ptr = &rw_vec;
	}

	PoissonGLMMData dat(X, y, gid_v, n_gh, rw_ptr);

	// Initialize
	Eigen::VectorXd par(total);
	if (warm_start_params.isNotNull()) {
		NumericVector sp(warm_start_params);
		if (sp.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sp[i];
		}
	} else if (smart_cold_start) {
		// Init: beta via OLS on log(y+0.5), log_sigma = -3
		Eigen::VectorXd log_y_safe = (y.array() + 0.5).log().matrix();
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
		Eigen::VectorXd beta_init = qr.solve(log_y_safe);
		par.head(p) = beta_init;
		par[total - 1] = -3.0;
	} else {
		par.head(p).setZero();
		par[total - 1] = -3.0;
	}

	PoissonGLMMObjective obj(dat);
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par       = fit.params;
		neg_ll    = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

		const double pen = soft_barrier(par[total - 1]);
		const double true_neg_ll = neg_ll - pen;
		if (estimate_only) {
			return List::create(
				Named("b")          = par.head(p),
				Named("log_sigma")  = par[total - 1],
				Named("ssq_b_T")    = NA_REAL,
				Named("vcov")       = R_NilValue,
				Named("converged")  = converged,
				Named("neg_loglik") = true_neg_ll,
				Named("fisher_information") = R_NilValue
			);
		}

		Eigen::MatrixXd information = obj.hessian(par);
		information(total - 1, total - 1) -= soft_barrier_hessian(par[total - 1]);

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
		if (converged) {
		Eigen::MatrixXd H_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
			vcov = expand_free_covariance(total, fixed_spec, inv_free, true);
			if (vcov.allFinite() && j_T < p) ssq_b_T = vcov(j_T, j_T);
		}
	}

	return List::create(
		Named("b")          = par.head(p),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("vcov")       = vcov,
		Named("converged")  = converged,
		Named("neg_loglik") = true_neg_ll,
		Named("fisher_information") = information
	);
}

// ── R-exported: score (gradient of log_lik) at arbitrary par ─────────────────
// [[Rcpp::export]]
SEXP get_poisson_glmm_score_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	const NumericVector& par_r,
	int n_gh = 20
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < group_id.size(); ++i) gid_v[i] = group_id[i];
	PoissonGLMMData dat(X, y, gid_v, n_gh);
	PoissonGLMMObjective obj(dat);
	Eigen::VectorXd grad;
	obj(par, grad);
	return wrap(-grad);
}

// ── R-exported: observed information (Hessian of neg_ll) at par ─────────────
// [[Rcpp::export]]
SEXP get_poisson_glmm_hessian_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	const NumericVector& par_r,
	int n_gh = 20
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < group_id.size(); ++i) gid_v[i] = group_id[i];
	PoissonGLMMData dat(X, y, gid_v, n_gh);
	PoissonGLMMObjective obj(dat);
	return wrap(-obj.hessian(par));
}
