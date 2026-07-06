// Logistic GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  logit P(Y_ij = 1 | u_i) = X_ij' beta + u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in [0,1]          (binary 0/1 or continuous proportion)
//   X includes an intercept column
//
// Parameter vector: par = [beta_0, beta_1(treatment), ..., beta_{p-1}, log_sigma]
//   Total length: p + 1
//
// log_sigma is parameterized with a soft-barrier penalty rather than a hard cut.
// This avoids the infinite-gradient issue at the boundary.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

namespace {

struct GHRule {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRule gauss_hermite_rule_log(int n) {
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

inline double log_sum_exp_v(const Eigen::Ref<const Eigen::VectorXd>& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}


// Soft barrier: adds a smooth penalty when |log_sigma| is large,
// avoiding the hard-wall numerical catastrophe.
// penalty(t) = 0 when |t| <= center, smoothly increases beyond.
inline double log_sigma_penalty(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double log_sigma_penalty_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double log_sigma_penalty_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

struct LogisticGLMMData {
	Eigen::MatrixXd X_s;        // n x p (includes intercept)
	std::vector<double> y_s;    // responses in [0,1], length n
	std::vector<int> grp_start; // warm_start_params index of each group
	std::vector<int> grp_size;  // size of each group
	Eigen::VectorXd grp_y_sum;  // sum of responses in each group
	int n, p, G, max_group_size;
	GHRule gh;

	LogisticGLMMData(
		const Eigen::Ref<const Eigen::MatrixXd>& X,
		const Eigen::Ref<const Eigen::VectorXd>& y,
		const std::vector<int>& group_id,
		int n_gh
	) : n((int)X.rows()), p((int)X.cols()), gh(gauss_hermite_rule_log(n_gh)) {

		// Sort by group_id for contiguous group access
		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id[a] < group_id[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y[ord[i]];
		}

		auto layout = build_contiguous_group_layout(n, [&](int i) { return group_id[ord[i]]; });
		grp_start = layout.start;
		grp_size = layout.size;
		G = layout.G;
		max_group_size = layout.max_size;
		grp_y_sum.resize(G);
		for (int gi = 0; gi < G; ++gi) {
			double sum_y = 0.0;
			for (int r = 0; r < grp_size[gi]; ++r) sum_y += y_s[grp_start[gi] + r];
			grp_y_sum[gi] = sum_y;
		}
	}
};

class LogisticGLMMObjective {
	const LogisticGLMMData& dat;
	Eigen::VectorXd m_eta_all;
	Eigen::VectorXd m_group_y_eta0;
	Eigen::VectorXd m_ll_g_vec;
	Eigen::MatrixXd m_log_terms_mat;
	Eigen::MatrixXd m_mu_all_k_mat;
	Eigen::MatrixXd m_mu_group_sum_mat;
	Eigen::VectorXd m_weighted_res;
	const int m_n_nodes;

public:
	explicit LogisticGLMMObjective(const LogisticGLMMData& d)
		: dat(d),
		  m_eta_all(d.n),
		  m_group_y_eta0(d.G),
		  m_ll_g_vec(d.G),
		  m_log_terms_mat(d.G, (int)d.gh.nodes.size()),
		  m_mu_all_k_mat(d.n, (int)d.gh.nodes.size()),
		  m_mu_group_sum_mat(d.G, (int)d.gh.nodes.size()),
		  m_weighted_res(std::max(1, d.max_group_size)),
		  m_n_nodes((int)d.gh.nodes.size()) {}

	double value(const Eigen::Ref<const Eigen::VectorXd>& par) const {
		const double log_sigma = par[dat.p];
		// Soft barrier for |log_sigma| > 5 (instead of hard cut)
		const double pen = log_sigma_penalty(log_sigma);
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		const Eigen::Map<const Eigen::VectorXd> y_all(dat.y_s.data(), dat.n);
		double total_ll = 0.0;
		Eigen::VectorXd log_terms(n_nodes);
		for (int gi = 0; gi < dat.G; ++gi) {
			const int gs = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];
			const Eigen::ArrayXd eta0 = (dat.X_s.middleRows(gs, sz) * beta).array();
			const double y_eta0 = y_all.segment(gs, sz).dot(eta0.matrix());
			const double y_sum  = dat.grp_y_sum[gi];
			for (int k = 0; k < n_nodes; ++k) {
				// vectorized softplus over the group's eta (was a scalar r-loop with log1pexp_s)
				log_terms[k] = dat.gh.log_norm_weights[k] + y_eta0 + b_vals[k] * y_sum
				             - log1pexp_array_fast(eta0 + b_vals[k]).sum();
			}
			const double ll_g = log_sum_exp_v(log_terms);
			if (!std::isfinite(ll_g)) return 1e50;
			total_ll += ll_g;
		}
		return -total_ll + pen;
	}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& par, Eigen::Ref<VectorXd> grad) {
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const Eigen::Map<const Eigen::VectorXd> y_all(dat.y_s.data(), dat.n);

		m_eta_all.noalias() = dat.X_s * beta;
		for (int gi = 0; gi < dat.G; ++gi) {
			m_group_y_eta0[gi] = y_all.segment(dat.grp_start[gi], dat.grp_size[gi]).dot(
				m_eta_all.segment(dat.grp_start[gi], dat.grp_size[gi])
			);
		}

		for (int k = 0; k < m_n_nodes; ++k) {
			for (int gi = 0; gi < dat.G; ++gi) {
				const int warm_start_params = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const Eigen::ArrayXd eta_g_k = m_eta_all.segment(warm_start_params, sz).array() + b_vals[k];
				const Eigen::ArrayXd mu_g_k = plogis_array_safe(eta_g_k);
				m_mu_all_k_mat.block(warm_start_params, k, sz, 1) = mu_g_k.matrix();
				m_mu_group_sum_mat(gi, k) = mu_g_k.sum();
				m_log_terms_mat(gi, k) = dat.gh.log_norm_weights[k]
					+ m_group_y_eta0[gi]
					+ b_vals[k] * dat.grp_y_sum[gi]
					- log1pexp_array_fast(eta_g_k).sum();
			}
		}

		double total_nll = log_sigma_penalty(log_sigma);
		for (int gi = 0; gi < dat.G; ++gi) {
			m_ll_g_vec[gi] = log_sum_exp_v(m_log_terms_mat.row(gi));
			if (!std::isfinite(m_ll_g_vec[gi])) { grad.setZero(); return 1e100; }
			total_nll -= m_ll_g_vec[gi];
		}

		grad.setZero();
		const double center = 5.0, scale = 10.0;
		const double d_pen = std::abs(log_sigma) - center;
		if (d_pen > 0.0) grad[dat.p] += 2.0 * scale * d_pen * (log_sigma > 0 ? 1.0 : -1.0);

		Eigen::VectorXd grad_beta = Eigen::VectorXd::Zero(dat.p);
		double grad_log_sigma = 0.0;

		for (int gi = 0; gi < dat.G; ++gi) {
			const int warm_start_params = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];
			m_weighted_res.head(sz).setZero();

			for (int k = 0; k < m_n_nodes; ++k) {
				const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
				if (pk < 1e-15) continue;
				m_weighted_res.head(sz).array() += pk * (
					y_all.segment(warm_start_params, sz).array() -
					m_mu_all_k_mat.block(warm_start_params, k, sz, 1).array()
				);
				const double res_sum = dat.grp_y_sum[gi] - m_mu_group_sum_mat(gi, k);
				grad_log_sigma -= pk * res_sum * b_vals[k];
			}

			grad_beta.noalias() -= dat.X_s.middleRows(warm_start_params, sz).transpose() * m_weighted_res.head(sz);
		}

		grad.head(dat.p) = grad_beta;
		grad[dat.p] += grad_log_sigma;

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& par) {
		const int total = dat.p + 1;
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		const Eigen::VectorXd eta_all = dat.X_s * beta;
		const Eigen::Map<const Eigen::VectorXd> y_all(dat.y_s.data(), dat.n);
		Eigen::VectorXd group_y_eta0(dat.G);
		for (int gi = 0; gi < dat.G; ++gi) {
			group_y_eta0[gi] = y_all.segment(dat.grp_start[gi], dat.grp_size[gi]).dot(
				eta_all.segment(dat.grp_start[gi], dat.grp_size[gi])
			);
		}

		Eigen::MatrixXd log_terms_mat(dat.G, n_nodes);
		Eigen::MatrixXd mu_all_k_mat(dat.n, n_nodes);
		Eigen::MatrixXd w_all_k_mat(dat.n, n_nodes);
		Eigen::MatrixXd mu_group_sum_mat(dat.G, n_nodes);
		Eigen::MatrixXd w_group_sum_mat(dat.G, n_nodes);
		Eigen::MatrixXd res_group_sum_mat(dat.G, n_nodes);
		for (int k = 0; k < n_nodes; ++k) {
			const Eigen::ArrayXd eta_all_k = eta_all.array() + b_vals[k];
			const Eigen::VectorXd mu_all_k = plogis_array_safe(eta_all_k).matrix();
			mu_all_k_mat.col(k) = mu_all_k;
			w_all_k_mat.col(k) = (mu_all_k.array() * (1.0 - mu_all_k.array())).matrix();
			for (int gi = 0; gi < dat.G; ++gi) {
				const int warm_start_params = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const Eigen::ArrayXd eta_g_k = eta_all.segment(warm_start_params, sz).array() + b_vals[k];
				const double mu_sum = mu_all_k_mat.block(warm_start_params, k, sz, 1).sum();
				mu_group_sum_mat(gi, k) = mu_sum;
				w_group_sum_mat(gi, k) = w_all_k_mat.block(warm_start_params, k, sz, 1).sum();
				res_group_sum_mat(gi, k) = dat.grp_y_sum[gi] - mu_sum;
				log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] +
				                       group_y_eta0[gi] +
				                       b_vals[k] * dat.grp_y_sum[gi] -
				                       log1pexp_array_fast(eta_g_k).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(dat.G);
		for (int gi = 0; gi < dat.G; ++gi) ll_g_vec[gi] = log_sum_exp_v(log_terms_mat.row(gi));

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);
		H(dat.p, dat.p) = log_sigma_penalty_hessian(log_sigma);

		Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(total, total);

		for (int k = 0; k < n_nodes; k++) {
			Eigen::VectorXd pk_vec(dat.G);
			for (int gi = 0; gi < dat.G; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);

			const double node_factor = b_vals[k];
			for (int gi = 0; gi < dat.G; gi++) {
				const double pk = pk_vec[gi];
				if (pk < 1e-15) continue;
				const int warm_start_params = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const auto Xg = dat.X_s.middleRows(warm_start_params, sz);
				const auto w_seg = w_all_k_mat.col(k).segment(warm_start_params, sz);
				const double sum_w_k_gi = w_group_sum_mat(gi, k);
				const double sum_res_k_gi = res_group_sum_mat(gi, k);
				E_Hik_sum.topLeftCorner(dat.p, dat.p).noalias() -= pk * weighted_crossprod(Xg, w_seg);
				
				E_Hik_sum(dat.p, dat.p) += pk * ((-sum_w_k_gi * node_factor * node_factor * sigma + sum_res_k_gi * node_factor) * sigma);
				
				Eigen::VectorXd d2L_db_dlogsigma_k_gi = -(Xg.transpose() * w_seg) * (node_factor * sigma);
				E_Hik_sum.block(0, dat.p, dat.p, 1).noalias() += pk * d2L_db_dlogsigma_k_gi;
			}
		}

			for (int gi = 0; gi < dat.G; gi++) {
				Eigen::VectorXd beta_avg_gi = Eigen::VectorXd::Zero(dat.p);
				double sigma_avg_gi = 0.0;
				Eigen::MatrixXd beta_beta_gi = Eigen::MatrixXd::Zero(dat.p, dat.p);
				Eigen::VectorXd beta_sigma_gi = Eigen::VectorXd::Zero(dat.p);
				double sigma_sigma_gi = 0.0;
				const int warm_start_params = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const Eigen::MatrixXd Xg = dat.X_s.middleRows(warm_start_params, sz);

				for (int k = 0; k < n_nodes; k++) {
					const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
					if (pk < 1e-15) continue;
					const Eigen::VectorXd g_beta = Xg.transpose() * (y_all.segment(warm_start_params, sz) - mu_all_k_mat.col(k).segment(warm_start_params, sz));
					const double g_sigma = res_group_sum_mat(gi, k) * b_vals[k];

					beta_avg_gi.noalias() += pk * g_beta;
					sigma_avg_gi += pk * g_sigma;
					beta_beta_gi.noalias() += pk * (g_beta * g_beta.transpose());
					beta_sigma_gi.noalias() += pk * (g_beta * g_sigma);
					sigma_sigma_gi += pk * g_sigma * g_sigma;
				}
				E_GiGiT_sum.topLeftCorner(dat.p, dat.p).noalias() += beta_beta_gi;
				E_GiGiT_sum.block(0, dat.p, dat.p, 1).noalias() += beta_sigma_gi;
				E_GiGiT_sum.block(dat.p, 0, 1, dat.p).noalias() += beta_sigma_gi.transpose();
				E_GiGiT_sum(dat.p, dat.p) += sigma_sigma_gi;
				G_avg_outer_sum.topLeftCorner(dat.p, dat.p).noalias() += beta_avg_gi * beta_avg_gi.transpose();
				G_avg_outer_sum.block(0, dat.p, dat.p, 1).noalias() += beta_avg_gi * sigma_avg_gi;
				G_avg_outer_sum.block(dat.p, 0, 1, dat.p).noalias() += (beta_avg_gi * sigma_avg_gi).transpose();
				G_avg_outer_sum(dat.p, dat.p) += sigma_avg_gi * sigma_avg_gi;
			}

		for (int r = 0; r < dat.p; r++) for (int c = 0; c < r; c++) E_Hik_sum(r, c) = E_Hik_sum(c, r);
		E_Hik_sum.block(dat.p, 0, 1, dat.p) = E_Hik_sum.block(0, dat.p, dat.p, 1).transpose();

		H -= (E_Hik_sum + E_GiGiT_sum - G_avg_outer_sum);
		return H;
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_logistic_glmm_score_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 20
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> gid_v(group_id.data(), group_id.data() + group_id.size());

	LogisticGLMMData dat(X, y, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);

	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	grad[X.cols()] -= log_sigma_penalty_grad(params[X.cols()]);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_glmm_hessian_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 20
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> gid_v(group_id.data(), group_id.data() + group_id.size());

	LogisticGLMMData dat(X, y, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);

	Eigen::MatrixXd information = obj.hessian(params);
	information(X.cols(), X.cols()) -= log_sigma_penalty_hessian(params[X.cols()]);
	return -information;
}

// [[Rcpp::export]]
double get_logistic_glmm_neg_loglik_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 20
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> gid_v(group_id.data(), group_id.data() + group_id.size());

	LogisticGLMMData dat(X, y, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);
	return likelihood_value(obj, params) - log_sigma_penalty(params[X.cols()]);
}

// [[Rcpp::export]]
List fast_logistic_glmm_cpp(
	SEXP X_r,       // n x p, includes intercept; treatment at col j_T (0-based)
	SEXP y_r,       // responses in [0,1], length n
	SEXP group_id_r,// group IDs, sorted internally
	int j_T,                        // 0-based treatment column index in X
	Nullable<NumericVector> warm_start_params = R_NilValue,
	bool smart_cold_start = true,
	bool estimate_only = false,
	int n_gh = 20,
	int maxit = 300,
	double eps_g = 1e-6,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	const int n = (int)X.rows();
	const int p = (int)X.cols();
	const int total = p + 1; // betas + log_sigma

	std::vector<int> gid_v(group_id.data(), group_id.data() + n);

	LogisticGLMMData dat(X, y, gid_v, n_gh);

	// Initialize
	Eigen::VectorXd par(total);
	if (warm_start_params.isNotNull()) {
		NumericVector sp(warm_start_params);
		if (sp.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sp[i];
		}
	} else if (smart_cold_start) {
		// Logistic smart warm_start_params: OLS on y
		par.head(p) = ols_smart_cold_start_beta(X, y);
		par[total - 1] = -3.0;
	} else {
		par.head(p).setZero();
		par[total - 1] = -3.0;
	}

	LogisticGLMMObjective obj(dat);
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll = NA_REAL;
	int niter = maxit;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par = fit.params;
		neg_ll = fit.value;
		niter = fit.niter;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	// Remove soft-barrier penalty from neg_ll so it reflects the true neg log-likelihood
	double true_neg_ll = neg_ll;
	if (estimate_only) {
		const double pen = log_sigma_penalty(par[total - 1]);
		true_neg_ll = neg_ll - pen;
		return List::create(
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("vcov")       = R_NilValue,
			Named("score")      = R_NilValue,
			Named("observed_information") = R_NilValue,
			Named("information") = R_NilValue,
			Named("information_type") = "observed",
			Named("hessian")    = R_NilValue,
			Named("converged")  = converged,
			Named("neg_loglik") = true_neg_ll,
			Named("neg_ll")     = true_neg_ll,
			Named("loglik")     = R_finite(true_neg_ll) ? -true_neg_ll : NA_REAL,
			Named("fisher_information") = R_NilValue
		);
	}

	Eigen::VectorXd score = Eigen::VectorXd::Constant(total, NA_REAL);
	Eigen::MatrixXd information = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	try {
		const double pen = log_sigma_penalty(par[total - 1]);
		true_neg_ll = neg_ll - pen;
		obj(par, score);
		score[total - 1] -= log_sigma_penalty_grad(par[total - 1]);
		score = -score;
		information = obj.hessian(par);
		information(total - 1, total - 1) -= log_sigma_penalty_hessian(par[total - 1]);
	} catch (...) {
		converged = false;
	}

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	if (converged) {
		Eigen::MatrixXd information_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::MatrixXd cov_free = covariance_from_information(information_free);
		vcov = expand_free_covariance(total, fixed_spec, cov_free, true);
		if (j_T < p) {
			ssq_b_T = vcov(j_T, j_T);
		}
	}

	return List::create(
		Named("params")     = par,
		Named("b")          = par.head(p),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("vcov")       = vcov,
		Named("score")      = score,
		Named("observed_information") = information,
		Named("information") = information,
		Named("information_type") = "observed",
		Named("hessian")    = -information,
		Named("converged")  = converged,
		Named("neg_loglik") = true_neg_ll,
		Named("neg_ll")     = true_neg_ll,
		Named("loglik")     = R_finite(true_neg_ll) ? -true_neg_ll : NA_REAL,
		Named("fisher_information") = information
	);
}
