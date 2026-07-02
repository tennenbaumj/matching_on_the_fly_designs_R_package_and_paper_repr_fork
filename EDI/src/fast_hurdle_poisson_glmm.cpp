// Zero-truncated Poisson GLMM for the count component of a hurdle Poisson GLMM.
//
// The hurdle likelihood factors into independent zi (logistic) and count
// (zero-truncated Poisson) components with separate random intercepts.
// This file handles the count component via Gauss-Hermite quadrature.
//
// Model (count component, positive observations only):
//   y_ij | u_i ~ TruncPoisson(exp(X_ij' beta + u_i)),  y_ij >= 1
//   u_i ~ N(0, sigma^2)
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

struct GHRuleHP {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleHP gauss_hermite_rule_hp(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleHP rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
		                    - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_hp(const Eigen::Ref<const Eigen::VectorXd>& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline double soft_barrier_hp(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double soft_barrier_hp_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double soft_barrier_hp_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

template<typename XDerived, typename WDerived, typename OutDerived>
inline void crossprod_rhs_assign_hp(const Eigen::MatrixBase<XDerived>& X,
                                    const Eigen::MatrixBase<WDerived>& w,
                                    Eigen::MatrixBase<OutDerived>& out) {
	const int n = X.rows();
	const int p = X.cols();
	out.derived().setZero();
	for (int j = 0; j < p; ++j) {
		double acc = 0.0;
		for (int i = 0; i < n; ++i) acc += X(i, j) * w(i);
		out(j) = acc;
	}
}

template<typename XDerived, typename WDerived, typename OutDerived>
inline void weighted_crossprod_assign_hp(const Eigen::MatrixBase<XDerived>& X,
                                         const Eigen::MatrixBase<WDerived>& w,
                                         Eigen::MatrixBase<OutDerived>& out) {
	const int n = X.rows();
	const int p = X.cols();
	out.derived().setZero();
	for (int j = 0; j < p; ++j) {
		for (int k = j; k < p; ++k) {
			double acc = 0.0;
			for (int i = 0; i < n; ++i) acc += X(i, j) * w(i) * X(i, k);
			out(j, k) = acc;
			if (k != j) out(k, j) = acc;
		}
	}
}

template<typename XDerived, typename WDerived>
inline void crossprod_rhs_to_col_hp(const Eigen::MatrixBase<XDerived>& X,
                                    const Eigen::MatrixBase<WDerived>& w,
                                    Eigen::MatrixXd& out,
                                    int col,
                                    double scale) {
	const int n = X.rows();
	const int p = X.cols();
	for (int j = 0; j < p; ++j) {
		double acc = 0.0;
		for (int i = 0; i < n; ++i) acc += X(i, j) * w(i);
		out(j, col) = acc * scale;
	}
}

struct HurdlePoissonGLMMData {
	Eigen::MatrixXd X_s;
	Eigen::VectorXd y_s;
	Eigen::VectorXd log_fact_y;
	std::vector<int> grp_start;
	std::vector<int> grp_size;
	int n, p, G, max_grp_sz;
	GHRuleHP gh;

	HurdlePoissonGLMMData(
		const Eigen::Ref<const Eigen::MatrixXd>& X,
		const Eigen::Ref<const Eigen::VectorXd>& y,
		const std::vector<int>& group_id,
		int n_gh
	) : n((int)X.rows()), p((int)X.cols()), gh(gauss_hermite_rule_hp(n_gh)) {

		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id[a] < group_id[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		log_fact_y.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y[ord[i]];
			log_fact_y[i] = std::lgamma(y_s[i] + 1.0);
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
		max_grp_sz = (G > 0) ? *std::max_element(grp_size.begin(), grp_size.end()) : 1;
	}
};

class HurdlePoissonGLMMObjective {
	const HurdlePoissonGLMMData& dat;
	// Preallocated column-major buffers: element [i + k*max_grp_sz] = obs i at node k.
	// Sized once at construction; reused across every operator() and hessian() call
	// to eliminate G*K heap allocations per optimizer step.
	std::vector<double> m_lam;   // exp(eta0[i] + b_vals[k])
	std::vector<double> m_eneg;  // exp(-lam[i,k]) — shared by ll and gradient
	std::vector<double> m_w;     // per-obs gradient weight accumulator
	std::vector<double> m_exp_bvals;
	Eigen::VectorXd m_b_vals;
	Eigen::VectorXd m_log_terms;
	Eigen::VectorXd m_eta0;
	Eigen::VectorXd m_res_k;
	Eigen::VectorXd m_d2e_k;
	Eigen::VectorXd m_G_ik;
	Eigen::MatrixXd m_H_ik;
	Eigen::MatrixXd m_E_Hik;
	Eigen::MatrixXd m_E_GiGiT;
	Eigen::VectorXd m_G_avg;

public:
	explicit HurdlePoissonGLMMObjective(const HurdlePoissonGLMMData& d) : dat(d) {
		const int K = (int)d.gh.nodes.size();
		const int total = d.p + 1;
		m_lam.resize(d.max_grp_sz * K);
		m_eneg.resize(d.max_grp_sz * K);
		m_w.resize(d.max_grp_sz);
		m_exp_bvals.resize(K);
		m_b_vals.resize(K);
		m_log_terms.resize(K);
		m_eta0.resize(d.max_grp_sz);
		m_res_k.resize(d.max_grp_sz);
		m_d2e_k.resize(d.max_grp_sz);
		m_G_ik.resize(total);
		m_H_ik.resize(total, total);
		m_E_Hik.resize(total, total);
		m_E_GiGiT.resize(total, total);
		m_G_avg.resize(total);
	}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& par, Eigen::Ref<Eigen::VectorXd> grad) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const int p = dat.p;
		const int K = (int)dat.gh.nodes.size();
		m_b_vals.noalias() = (std::sqrt(2.0) * sigma) * dat.gh.nodes;

		// Precompute exp(b_vals[k]) once per optimizer step (K exps amortized over G groups).
		for (int k = 0; k < K; ++k) m_exp_bvals[k] = std::exp(m_b_vals[k]);

		double total_nll = soft_barrier_hp(log_sigma);
		grad.setZero();
		const double d = std::abs(log_sigma) - 5.0;
		if (d > 0.0) grad[p] += 20.0 * d * (log_sigma > 0 ? 1.0 : -1.0);

		for (int gi = 0; gi < dat.G; ++gi) {
			const int gs = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];
			const auto Xg = dat.X_s.middleRows(gs, sz);
			m_eta0.head(sz).noalias() = Xg * par.head(p);
			const auto eta0 = m_eta0.head(sz);
			const double* y_g = dat.y_s.data() + gs;
			const double* lfg = dat.log_fact_y.data() + gs;

			// Forward pass: multiplicative exp decomposition.
			// lam[i,k] = exp(eta0[i]) * exp_bvals[k]  (no exp call per pair).
			// eneg[i,k] = exp(-lam[i,k])  — computed once, reused in ll and gradient.
			for (int k = 0; k < K; ++k) {
				const double eb    = m_exp_bvals[k];
				double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;
				for (int i = 0; i < sz; ++i) {
					const double lam = std::exp(std::min(eta0[i], 700.0)) * eb;
					lam_k[i]  = lam;
					eneg_k[i] = std::exp(-lam);
				}
			}

			// Log-likelihood contribution per node using shared lam/eneg.
			for (int k = 0; k < K; ++k) {
				const double bk     = m_b_vals[k];
				const double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				const double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;
				double ll_k = dat.gh.log_norm_weights[k];
				for (int i = 0; i < sz; ++i) {
					const double lam    = lam_k[i];
					const double eneg   = eneg_k[i];
					const double eta_ki = eta0[i] + bk;
					// log(1-exp(-lam)): for tiny lam use log(lam)=eta_ki (avoids cancellation).
					const double lne = (lam < 1e-10) ? eta_ki : std::log1p(-eneg);
					ll_k += y_g[i] * eta_ki - lam - lfg[i] - lne;
				}
				m_log_terms[k] = ll_k;
			}

			const double ll_g = log_sum_exp_hp(m_log_terms);
			if (!std::isfinite(ll_g)) { grad.setZero(); return 1e100; }
			total_nll -= ll_g;

			// Gradient: accumulate per-obs weights across all K nodes, then one GEMV per group.
			// Previously: K separate Xg' * res_k GEMV calls.
			double* w_g = m_w.data();
			std::fill(w_g, w_g + sz, 0.0);
			double w_sigma = 0.0;

			for (int k = 0; k < K; ++k) {
				const double post_k = std::exp(m_log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;
				const double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				const double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;
				const double pbk = post_k * m_b_vals[k];
				double score_sum = 0.0;
				for (int i = 0; i < sz; ++i) {
					const double lam  = lam_k[i];
					const double eneg = eneg_k[i];
					const double si   = (lam > 30.0) ? (y_g[i] - lam)
					                  : (lam < 1e-8)  ? (y_g[i] - 1.0)
					                  :                  (y_g[i] - lam / (1.0 - eneg));
					w_g[i]    += post_k * si;
					score_sum += si;
				}
				w_sigma += pbk * score_sum;
			}

			const Eigen::Map<const Eigen::VectorXd> w_vec(w_g, sz);
			grad.head(p).noalias() -= dat.X_s.middleRows(gs, sz).transpose() * w_vec;
			grad[p] -= w_sigma;
		}
		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& par) {
		const int total = dat.p + 1;
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const int p = dat.p;
		const int K = (int)dat.gh.nodes.size();
		m_b_vals.noalias() = (std::sqrt(2.0) * sigma) * dat.gh.nodes;

		for (int k = 0; k < K; ++k) m_exp_bvals[k] = std::exp(m_b_vals[k]);

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);
		H(dat.p, dat.p) = soft_barrier_hp_hessian(log_sigma);

		for (int gi = 0; gi < dat.G; gi++) {
			const int gs = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];
			const auto Xg = dat.X_s.middleRows(gs, sz);
			m_eta0.head(sz).noalias() = Xg * par.head(p);
			const auto eta0 = m_eta0.head(sz);
			const double* y_g = dat.y_s.data() + gs;
			const double* lfg = dat.log_fact_y.data() + gs;

			for (int k = 0; k < K; ++k) {
				const double eb = m_exp_bvals[k];
				double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;
				for (int i = 0; i < sz; ++i) {
					const double lam = std::exp(std::min(eta0[i], 700.0)) * eb;
					lam_k[i]  = lam;
					eneg_k[i] = std::exp(-lam);
				}
			}
			for (int k = 0; k < K; ++k) {
				const double bk     = m_b_vals[k];
				const double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				const double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;
				double ll_k = dat.gh.log_norm_weights[k];
				for (int i = 0; i < sz; ++i) {
					const double lam    = lam_k[i];
					const double eneg   = eneg_k[i];
					const double eta_ki = eta0[i] + bk;
					const double lne    = (lam < 1e-10) ? eta_ki : std::log1p(-eneg);
					ll_k += y_g[i] * eta_ki - lam - lfg[i] - lne;
				}
				m_log_terms[k] = ll_k;
			}
			const double ll_g = log_sum_exp_hp(m_log_terms);

			m_E_Hik.setZero();
			m_E_GiGiT.setZero();
			m_G_avg.setZero();

			for (int k = 0; k < K; k++) {
				const double pk = std::exp(m_log_terms[k] - ll_g);
				if (pk < 1e-15) continue;

				const double* lam_k  = m_lam.data()  + k * dat.max_grp_sz;
				const double* eneg_k = m_eneg.data() + k * dat.max_grp_sz;

				auto res_k = m_res_k.head(sz);
				auto d2e_k = m_d2e_k.head(sz);
				for (int i = 0; i < sz; ++i) {
					const double lam  = lam_k[i];
					const double eneg = eneg_k[i];
					if (lam > 30.0) {
						res_k[i] = y_g[i] - lam;
						d2e_k[i] = -lam;
					} else if (lam < 1e-8) {
						res_k[i] = y_g[i] - 1.0;
						d2e_k[i] = -lam / 2.0;
					} else {
						const double one_minus = 1.0 - eneg;
						res_k[i] = y_g[i] - lam / one_minus;
						d2e_k[i] = -lam * (one_minus - lam * eneg) / (one_minus * one_minus);
					}
				}

				m_G_ik.setZero();
				m_H_ik.setZero();
				const double sum_res = res_k.sum();
				const double sum_d2e = d2e_k.sum();

				auto G_beta = m_G_ik.head(p);
				crossprod_rhs_assign_hp(Xg, res_k, G_beta);
				auto H_beta = m_H_ik.topLeftCorner(p, p);
				weighted_crossprod_assign_hp(Xg, d2e_k, H_beta);

				const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
				m_G_ik[dat.p] = sum_res * node_factor * sigma;

				m_H_ik(dat.p, dat.p) = (sum_d2e * node_factor * node_factor * sigma + sum_res * node_factor) * sigma;

				crossprod_rhs_to_col_hp(Xg, d2e_k, m_H_ik, dat.p, node_factor * sigma);

				for (int r1 = 0; r1 < total; r1++) for (int c1 = 0; r1 > c1; c1++) m_H_ik(r1, c1) = m_H_ik(c1, r1);

				m_E_Hik.noalias() += pk * m_H_ik;
				m_G_avg.noalias() += pk * m_G_ik;
				m_E_GiGiT.noalias() += pk * (m_G_ik * m_G_ik.transpose());
			}
			H.noalias() -= m_E_Hik;
			H.noalias() -= m_E_GiGiT;
			H.noalias() += m_G_avg * m_G_avg.transpose();
		}
		return H;
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_hurdle_poisson_glmm_score_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 7
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	grad[X.cols()] -= soft_barrier_hp_grad(params[X.cols()]);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_hurdle_poisson_glmm_hessian_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 7
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	Eigen::MatrixXd information = obj.hessian(params);
	information(X.cols(), X.cols()) -= soft_barrier_hp_hessian(params[X.cols()]);
	return -information;
}

// [[Rcpp::export]]
double get_hurdle_poisson_glmm_neg_loglik_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	SEXP params_sexp,
	int n_gh = 7
) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	IntegerVector group_id_int(group_id_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_int.begin(), group_id_int.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	return likelihood_value(obj, params) - soft_barrier_hp(params[X.cols()]);
}

// [[Rcpp::export]]
List fast_hurdle_poisson_glmm_cpp(
	SEXP X_r,
	SEXP y_r,
	SEXP group_id_r,
	int j_T,
	Nullable<NumericVector> warm_start_params = R_NilValue,
	bool smart_cold_start = true,
	bool estimate_only = false,
	int n_gh = 7,
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
	const int n_all = (int)X.rows();
	const int p     = (int)X.cols();
	const int total = p + 1;

	// Filter to positive-count observations (truncated Poisson component)
	std::vector<int> pos_idx;
	pos_idx.reserve(n_all);
	for (int i = 0; i < n_all; ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = (int)pos_idx.size();

	if (n_pos <= p) {
		return List::create(
			Named("params")     = Eigen::VectorXd::Constant(total, NA_REAL),
			Named("b")          = Eigen::VectorXd::Constant(p, NA_REAL),
			Named("log_sigma")  = NA_REAL,
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	Eigen::MatrixXd X_pos(n_pos, p);
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k]     = y[i];
		gid_pos[k]   = group_id[i];
	}

	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);

	Eigen::VectorXd par(total);
	if (warm_start_params.isNotNull()) {
		par = as<Eigen::VectorXd>(NumericVector(warm_start_params));
	} else if (smart_cold_start) {
		par.setZero();
		Eigen::VectorXd log_y = y_pos.array().log().matrix();
		Eigen::VectorXd legacy_beta = Eigen::VectorXd::Zero(p);
		if (try_safe_ols_solve(X_pos, log_y, legacy_beta)) {
			par.head(p) = legacy_beta;
		}
		par[total - 1] = -3.0;
	} else {
		par.head(p).setZero();
		par[total - 1] = -3.0;
	}

	HurdlePoissonGLMMObjective obj(dat);
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
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const double pen         = soft_barrier_hp(par[total - 1]);
	const double true_neg_ll = neg_ll - pen;

	// Early return: skip score/Hessian computation when only point estimates are needed.
	if (estimate_only) {
		return List::create(
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("vcov")       = Eigen::MatrixXd::Constant(total, total, NA_REAL),
			Named("converged")  = converged,
			Named("neg_loglik") = true_neg_ll,
			Named("neg_ll")     = true_neg_ll,
			Named("loglik")     = R_finite(true_neg_ll) ? -true_neg_ll : NA_REAL
		);
	}

	Eigen::VectorXd score(total);
	obj(par, score);
	score[total - 1] -= soft_barrier_hp_grad(par[total - 1]);
	score = -score;
	Eigen::MatrixXd information = obj.hessian(par);
	information(total - 1, total - 1) -= soft_barrier_hp_hessian(par[total - 1]);

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	if (converged) {
		Eigen::MatrixXd information_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::MatrixXd cov_free = covariance_from_information(information_free);
		vcov = expand_free_covariance(total, fixed_spec, cov_free, true);
		if (j_T < p) ssq_b_T = vcov(j_T, j_T);
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
