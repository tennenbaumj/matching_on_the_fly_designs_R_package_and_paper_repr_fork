#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {


inline double log_sum_exp_cpp(const Eigen::Ref<const Eigen::VectorXd>& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

struct GHRule {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRule gauss_hermite_rule(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	Eigen::VectorXd nodes = es.eigenvalues();
	Eigen::VectorXd weights = std::sqrt(M_PI) * es.eigenvectors().row(0).array().square().matrix();
	GHRule rule;
	rule.nodes = nodes;
	rule.log_norm_weights = weights.array().log() - 0.5 * std::log(M_PI);
	return rule;
}

class ClogitPlusGLMMObjective {
private:
	const Eigen::Map<const Eigen::MatrixXd> X_disc;
	const Eigen::Map<const Eigen::VectorXd> y_disc;
	const Eigen::Map<const Eigen::MatrixXd> X_conc;
	const Eigen::Map<const Eigen::VectorXd> y_conc;
	const Eigen::Map<const Eigen::VectorXi> group_conc;
	const int q;
	const bool has_discordant;
	const bool has_concordant;
	const GHRule gh;
	const double max_abs_log_sigma;
	const std::vector<int> m_grp_start_conc;
	const std::vector<int> m_grp_size_conc;
	const int m_G_conc;
	const int m_max_group_size_conc;
	const Eigen::VectorXd m_grp_y_sum_conc;
	mutable Eigen::VectorXd m_eta_conc_all;
	mutable Eigen::VectorXd m_grp_y_eta0_conc;
	Eigen::VectorXd m_ll_g_vec;
	mutable Eigen::MatrixXd m_log_terms_mat;
	Eigen::MatrixXd m_mu_conc_all_k_mat;
	Eigen::VectorXd m_weighted_res;
	mutable Eigen::ArrayXd m_eta_k_g;    // per-group scratch — eliminates per-call inner-loop allocs
	Eigen::VectorXd m_grad_beta_conc;    // beta gradient accumulator — eliminates per-call alloc
	// hessian() scratch — promoted from per-call locals
	mutable Eigen::MatrixXd m_w_conc_all_k_mat;  // n_conc × n_gh
	mutable Eigen::MatrixXd m_mu_group_sum_mat;  // G × n_gh
	mutable Eigen::MatrixXd m_w_group_sum_mat;   // G × n_gh
	mutable Eigen::MatrixXd m_res_group_sum_mat; // G × n_gh
	mutable Eigen::VectorXd m_hess_pk_vec;       // G
	mutable Eigen::VectorXd m_hess_d2L;          // n_beta
	mutable Eigen::VectorXd m_hess_beta_avg;     // n_beta
	mutable Eigen::MatrixXd m_hess_beta_beta;    // n_beta × n_beta
	mutable Eigen::VectorXd m_hess_beta_sigma;   // n_beta
	mutable Eigen::VectorXd m_hess_g_beta;       // n_beta

	static ContiguousGroupLayout build_concordant_layout(const Eigen::Ref<const Eigen::VectorXi>& group_conc_) {
		return build_contiguous_group_layout(group_conc_.size(), [&](int idx) { return group_conc_[idx]; });
	}

	static Eigen::VectorXd build_group_y_sum(const Eigen::Ref<const Eigen::VectorXd>& y_conc_,
	                                         const std::vector<int>& grp_start,
	                                         const std::vector<int>& grp_size) {
		Eigen::VectorXd out(grp_start.size());
		for (int gi = 0; gi < static_cast<int>(grp_start.size()); ++gi) {
			out[gi] = y_conc_.segment(grp_start[gi], grp_size[gi]).sum();
		}
		return out;
	}

public:
	ClogitPlusGLMMObjective(
		const Eigen::Map<const Eigen::MatrixXd>& X_disc_,
		const Eigen::Map<const Eigen::VectorXd>& y_disc_,
		const Eigen::Map<const Eigen::MatrixXd>& X_conc_,
		const Eigen::Map<const Eigen::VectorXd>& y_conc_,
		const Eigen::Map<const Eigen::VectorXi>& group_conc_,
		const bool has_discordant_,
		const bool has_concordant_,
		const int n_gh = 20,
		const double max_abs_log_sigma_ = 8.0
	) :
		X_disc(X_disc_), y_disc(y_disc_), X_conc(X_conc_), y_conc(y_conc_),
		group_conc(group_conc_), q(has_discordant_ ? (int)X_disc_.cols() : (int)X_conc_.cols()), has_discordant(has_discordant_),
		has_concordant(has_concordant_), gh(gauss_hermite_rule(n_gh)),
		max_abs_log_sigma(max_abs_log_sigma_),
			m_grp_start_conc(build_concordant_layout(group_conc_).start),
		m_grp_size_conc(build_concordant_layout(group_conc_).size),
		m_G_conc(static_cast<int>(m_grp_start_conc.size())),
		m_max_group_size_conc(build_concordant_layout(group_conc_).max_size),
		m_grp_y_sum_conc(build_group_y_sum(y_conc_, m_grp_start_conc, m_grp_size_conc)),
		m_eta_conc_all(std::max(1, static_cast<int>(X_conc_.rows()))),
		m_grp_y_eta0_conc(std::max(1, static_cast<int>(m_grp_start_conc.size()))),
		m_ll_g_vec(std::max(1, static_cast<int>(m_grp_start_conc.size()))),
		m_log_terms_mat(std::max(1, static_cast<int>(m_grp_start_conc.size())), n_gh),
		m_mu_conc_all_k_mat(std::max(1, static_cast<int>(X_conc_.rows())), n_gh),
		m_weighted_res(std::max(1, build_concordant_layout(group_conc_).max_size)),
		m_eta_k_g(std::max(1, m_max_group_size_conc)),
		m_grad_beta_conc(std::max(1, static_cast<int>(X_conc_.cols()))),
		m_w_conc_all_k_mat(std::max(1, static_cast<int>(X_conc_.rows())), n_gh),
		m_mu_group_sum_mat(std::max(1, static_cast<int>(m_grp_start_conc.size())), n_gh),
		m_w_group_sum_mat(std::max(1, static_cast<int>(m_grp_start_conc.size())), n_gh),
		m_res_group_sum_mat(std::max(1, static_cast<int>(m_grp_start_conc.size())), n_gh),
		m_hess_pk_vec(std::max(1, static_cast<int>(m_grp_start_conc.size()))),
		m_hess_d2L(std::max(1, static_cast<int>(X_conc_.cols()))),
		m_hess_beta_avg(std::max(1, static_cast<int>(X_conc_.cols()))),
		m_hess_beta_beta(std::max(1, static_cast<int>(X_conc_.cols())),
		                 std::max(1, static_cast<int>(X_conc_.cols()))),
		m_hess_beta_sigma(std::max(1, static_cast<int>(X_conc_.cols()))),
		m_hess_g_beta(std::max(1, static_cast<int>(X_conc_.cols()))) {}

	double neg_clogit(const Eigen::Ref<const Eigen::VectorXd>& beta_no_intercept) const {
		if (!has_discordant) return 0.0;
		const Eigen::ArrayXd eta = (X_disc * beta_no_intercept).array();
		const Eigen::ArrayXd ll = y_disc.array() * eta - log1pexp_array_safe(eta);  // vectorized softplus
		if (!ll.allFinite()) return 1e100;
		return -ll.sum();
	}

	double neg_glmm(const Eigen::Ref<const Eigen::VectorXd>& par_full) const {
		if (!has_concordant) return 0.0;
		const int p_full = par_full.size();
		const double log_sigma = par_full[p_full - 1];
		if (!std::isfinite(log_sigma) || std::abs(log_sigma) > max_abs_log_sigma) return 1e100;
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
		const int n_nodes = (int)b_vals.size();

		m_eta_conc_all.noalias() = X_conc * par_full.head(p_full - 1);
		for (int gi = 0; gi < m_G_conc; ++gi) {
			const int warm_start_params = m_grp_start_conc[gi];
			const int sz = m_grp_size_conc[gi];
			m_grp_y_eta0_conc[gi] = y_conc.segment(warm_start_params, sz).dot(m_eta_conc_all.segment(warm_start_params, sz));
		}
		for (int k = 0; k < n_nodes; ++k) {
			for (int gi = 0; gi < m_G_conc; ++gi) {
				const int warm_start_params = m_grp_start_conc[gi];
				const int sz = m_grp_size_conc[gi];
				m_eta_k_g.head(sz) = m_eta_conc_all.segment(warm_start_params, sz).array() + b_vals[k];
				double lse_sum = 0.0;
				for (int j = 0; j < sz; ++j) lse_sum += log1pexp_safe(m_eta_k_g[j]);
				m_log_terms_mat(gi, k) = gh.log_norm_weights[k]
				                         + m_grp_y_eta0_conc[gi]
				                         + b_vals[k] * m_grp_y_sum_conc[gi]
				                         - lse_sum;
			}
		}
		double total_ll = 0.0;
		for (int gi = 0; gi < m_G_conc; ++gi) {
			const double ll_g = log_sum_exp_cpp(m_log_terms_mat.row(gi));
			if (!std::isfinite(ll_g)) return 1e100;
			total_ll += ll_g;
		}
		return -total_ll;
	}

	double value(const Eigen::Ref<const Eigen::VectorXd>& par) const {
		double out = 0.0;
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			out += neg_clogit(beta_no_intercept);
		}
		if (has_concordant) out += neg_glmm(par);
		if (!std::isfinite(out)) return 1e100;
		return out;
	}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& par, Eigen::VectorXd& grad) {
		const int p_full = par.size();
		grad.setZero(p_full);
		double total_nll = 0.0;

		// 1. Conditional Logistic component
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			
			const Eigen::VectorXd eta_d = X_disc * beta_no_intercept;
			const Eigen::VectorXd mu_d = plogis_array_safe(eta_d.array()).matrix();
			total_nll += (log1pexp_array_safe(eta_d.array()) - y_disc.array() * eta_d.array()).sum();

			Eigen::VectorXd grad_clogit = X_disc.transpose() * (mu_d - y_disc);
			if (has_concordant) {
				grad.segment(1, q) += grad_clogit;
			} else {
				grad.head(q) += grad_clogit;
			}
		}

		// 2. GLMM component
		if (has_concordant) {
			const double log_sigma = par[p_full - 1];
			const double sigma = std::exp(log_sigma);
			const Eigen::VectorXd beta = par.head(p_full - 1);
			const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
			const int n_nodes = (int)b_vals.size();

			m_eta_conc_all.noalias() = X_conc * beta;
			for (int gi = 0; gi < m_G_conc; ++gi) {
				const int warm_start_params = m_grp_start_conc[gi];
				const int sz = m_grp_size_conc[gi];
				m_grp_y_eta0_conc[gi] = y_conc.segment(warm_start_params, sz).dot(m_eta_conc_all.segment(warm_start_params, sz));
			}

			for (int k = 0; k < n_nodes; ++k) {
				for (int gi = 0; gi < m_G_conc; ++gi) {
					const int warm_start_params = m_grp_start_conc[gi];
					const int sz = m_grp_size_conc[gi];
					m_eta_k_g.head(sz) = m_eta_conc_all.segment(warm_start_params, sz).array() + b_vals[k];
					m_mu_conc_all_k_mat.col(k).segment(warm_start_params, sz) =
						(1.0 / (1.0 + (-m_eta_k_g.head(sz)).exp())).matrix();
					double lse_sum = 0.0;
					for (int j = 0; j < sz; ++j) lse_sum += log1pexp_safe(m_eta_k_g[j]);
					m_log_terms_mat(gi, k) = gh.log_norm_weights[k] +
					                       m_grp_y_eta0_conc[gi] +
					                       b_vals[k] * m_grp_y_sum_conc[gi] -
					                       lse_sum;
				}
			}

			for (int gi = 0; gi < m_G_conc; ++gi) {
				m_ll_g_vec[gi] = log_sum_exp_cpp(m_log_terms_mat.row(gi));
				total_nll -= m_ll_g_vec[gi];
			}

			m_grad_beta_conc.setZero();
			double grad_log_sigma_conc = 0.0;

			for (int gi = 0; gi < m_G_conc; ++gi) {
				const int warm_start_params = m_grp_start_conc[gi];
				const int sz = m_grp_size_conc[gi];
				m_weighted_res.head(sz).setZero();
				double grad_log_sigma_gi = 0.0;

				for (int k = 0; k < n_nodes; ++k) {
					const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
					if (pk < 1e-15) continue;
					m_weighted_res.head(sz).array() += pk * (
						y_conc.segment(warm_start_params, sz).array() -
						m_mu_conc_all_k_mat.block(warm_start_params, k, sz, 1).array()
					);
					const double res_sum = m_grp_y_sum_conc[gi] - m_mu_conc_all_k_mat.block(warm_start_params, k, sz, 1).sum();
					grad_log_sigma_gi += pk * res_sum * b_vals[k];
				}

				m_grad_beta_conc.noalias() -= X_conc.middleRows(warm_start_params, sz).transpose() * m_weighted_res.head(sz);
				grad_log_sigma_conc -= grad_log_sigma_gi;
			}
			grad.head(p_full - 1) += m_grad_beta_conc;
			grad[p_full - 1] += grad_log_sigma_conc;
		}

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& par) {
		const int p_full = par.size();
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(p_full, p_full);

		// 1. Conditional Logistic component
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			const Eigen::VectorXd eta_d = X_disc * beta_no_intercept;
			const Eigen::VectorXd mu_d = plogis_array_safe(eta_d.array()).matrix();
			Eigen::VectorXd w_d = (mu_d.array() * (1.0 - mu_d.array())).matrix();
			const Eigen::MatrixXd H_clogit = weighted_crossprod(X_disc, w_d);
			if (has_concordant) {
				H.block(1, 1, q, q).noalias() += H_clogit;
			} else {
				H.topLeftCorner(q, q).noalias() += H_clogit;
			}
		}

		// 2. GLMM component
		if (has_concordant) {
			const double log_sigma = par[p_full - 1];
			if (!std::isfinite(log_sigma) || std::abs(log_sigma) > max_abs_log_sigma) {
				H.setConstant(NA_REAL);
				return H;
			}

			const double sigma = std::exp(log_sigma);
			const Eigen::VectorXd beta = par.head(p_full - 1);
			const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
			const int n_nodes = (int)b_vals.size();
			const int n_beta = p_full - 1;

			// Reuse mutable members (all already sized in constructor)
			m_eta_conc_all.noalias() = X_conc * beta;
			for (int gi = 0; gi < m_G_conc; ++gi) {
				const int s = m_grp_start_conc[gi];
				const int sz = m_grp_size_conc[gi];
				m_grp_y_eta0_conc[gi] = y_conc.segment(s, sz).dot(m_eta_conc_all.segment(s, sz));
			}

			for (int k = 0; k < n_nodes; ++k) {
				// mu and w for all observations at node k — no per-k temporaries
				m_mu_conc_all_k_mat.col(k) =
					(1.0 / (1.0 + (-(m_eta_conc_all.array() + b_vals[k])).exp())).matrix();
				m_w_conc_all_k_mat.col(k) =
					(m_mu_conc_all_k_mat.col(k).array() * (1.0 - m_mu_conc_all_k_mat.col(k).array())).matrix();
				for (int gi = 0; gi < m_G_conc; ++gi) {
					const int s = m_grp_start_conc[gi];
					const int sz = m_grp_size_conc[gi];
					m_eta_k_g.head(sz) = m_eta_conc_all.segment(s, sz).array() + b_vals[k];
					double lse_sum = 0.0;
					for (int j = 0; j < sz; ++j) lse_sum += log1pexp_safe(m_eta_k_g[j]);
					const double mu_sum = m_mu_conc_all_k_mat.col(k).segment(s, sz).sum();
					m_mu_group_sum_mat(gi, k) = mu_sum;
					m_w_group_sum_mat(gi, k) = m_w_conc_all_k_mat.col(k).segment(s, sz).sum();
					m_res_group_sum_mat(gi, k) = m_grp_y_sum_conc[gi] - mu_sum;
					m_log_terms_mat(gi, k) = gh.log_norm_weights[k]
					                         + m_grp_y_eta0_conc[gi]
					                         + b_vals[k] * m_grp_y_sum_conc[gi]
					                         - lse_sum;
				}
			}

			for (int gi = 0; gi < m_G_conc; ++gi)
				m_ll_g_vec[gi] = log_sum_exp_cpp(m_log_terms_mat.row(gi));

			Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(p_full, p_full);
			Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(p_full, p_full);
			Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(p_full, p_full);

			for (int k = 0; k < n_nodes; k++) {
				for (int gi = 0; gi < m_G_conc; gi++)
					m_hess_pk_vec[gi] = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);

				for (int gi = 0; gi < m_G_conc; gi++) {
					const double pk = m_hess_pk_vec[gi];
					if (pk < 1e-15) continue;
					const int s = m_grp_start_conc[gi];
					const int sz = m_grp_size_conc[gi];
					// X_conc.middleRows() is MatrixBase<Derived> — no copy needed
					const auto w_seg = m_w_conc_all_k_mat.col(k).segment(s, sz);
					const double b = b_vals[k];
					E_Hik_sum.topLeftCorner(n_beta, n_beta).noalias() -=
						pk * weighted_crossprod(X_conc.middleRows(s, sz), w_seg);
					E_Hik_sum(p_full - 1, p_full - 1) +=
						pk * (b * m_res_group_sum_mat(gi, k) - b * b * m_w_group_sum_mat(gi, k));
					m_hess_d2L.noalias() = -b * (X_conc.middleRows(s, sz).transpose() * w_seg);
					E_Hik_sum.block(0, p_full - 1, n_beta, 1).noalias() += pk * m_hess_d2L;
				}
			}

			for (int gi = 0; gi < m_G_conc; gi++) {
				m_hess_beta_avg.setZero();
				double sigma_avg_gi = 0.0;
				m_hess_beta_beta.setZero();
				m_hess_beta_sigma.setZero();
				double sigma_sigma_gi = 0.0;
				const int s = m_grp_start_conc[gi];
				const int sz = m_grp_size_conc[gi];

				for (int k = 0; k < n_nodes; k++) {
					const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
					if (pk < 1e-15) continue;
					// g_beta: p-vector, computed directly into preallocated buffer
					m_hess_g_beta.noalias() = X_conc.middleRows(s, sz).transpose() *
						(y_conc.segment(s, sz) - m_mu_conc_all_k_mat.col(k).segment(s, sz));
					const double g_sigma = m_res_group_sum_mat(gi, k) * b_vals[k];

					m_hess_beta_avg.noalias() += pk * m_hess_g_beta;
					sigma_avg_gi += pk * g_sigma;
					m_hess_beta_beta.noalias() += pk * (m_hess_g_beta * m_hess_g_beta.transpose());
					m_hess_beta_sigma.noalias() += pk * (m_hess_g_beta * g_sigma);
					sigma_sigma_gi += pk * g_sigma * g_sigma;
				}
				E_GiGiT_sum.topLeftCorner(n_beta, n_beta).noalias() += m_hess_beta_beta;
				E_GiGiT_sum.block(0, p_full - 1, n_beta, 1).noalias() += m_hess_beta_sigma;
				E_GiGiT_sum.block(p_full - 1, 0, 1, n_beta).noalias() += m_hess_beta_sigma.transpose();
				E_GiGiT_sum(p_full - 1, p_full - 1) += sigma_sigma_gi;
				G_avg_outer_sum.topLeftCorner(n_beta, n_beta).noalias() +=
					m_hess_beta_avg * m_hess_beta_avg.transpose();
				G_avg_outer_sum.block(0, p_full - 1, n_beta, 1).noalias() +=
					m_hess_beta_avg * sigma_avg_gi;
				G_avg_outer_sum.block(p_full - 1, 0, 1, n_beta).noalias() +=
					(m_hess_beta_avg * sigma_avg_gi).transpose();
				G_avg_outer_sum(p_full - 1, p_full - 1) += sigma_avg_gi * sigma_avg_gi;
			}

			for (int r = 0; r < n_beta; r++) for (int c = 0; c < r; c++) E_Hik_sum(r, c) = E_Hik_sum(c, r);
			E_Hik_sum.block(p_full - 1, 0, 1, n_beta) = E_Hik_sum.block(0, p_full - 1, n_beta, 1).transpose();

			H.noalias() -= (E_Hik_sum + E_GiGiT_sum - G_avg_outer_sum);
		}

		return 0.5 * (H + H.transpose());
	}
};

} // namespace

// [[Rcpp::export]]
SEXP get_clogit_plus_glmm_score_cpp(
	const NumericMatrix& X_disc_r,
	const NumericVector& y_disc_r,
	const NumericMatrix& X_conc_r,
	const NumericVector& y_conc_r,
	const IntegerVector& group_conc_r,
	const NumericVector& params_r,
	bool has_discordant,
	bool has_concordant,
	double max_abs_log_sigma = 8.0
) {
	Eigen::Map<const Eigen::MatrixXd> X_disc(X_disc_r.begin(), X_disc_r.rows(), X_disc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_disc(y_disc_r.begin(), y_disc_r.size());
	Eigen::Map<const Eigen::MatrixXd> X_conc(X_conc_r.begin(), X_conc_r.rows(), X_conc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_conc(y_conc_r.begin(), y_conc_r.size());
	Eigen::Map<const Eigen::VectorXi> group_conc(group_conc_r.begin(), group_conc_r.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return wrap(-grad);
}

// [[Rcpp::export]]
SEXP get_clogit_plus_glmm_hessian_cpp(
	const NumericMatrix& X_disc_r,
	const NumericVector& y_disc_r,
	const NumericMatrix& X_conc_r,
	const NumericVector& y_conc_r,
	const IntegerVector& group_conc_r,
	const NumericVector& params_r,
	bool has_discordant,
	bool has_concordant,
	double max_abs_log_sigma = 8.0
) {
	Eigen::Map<const Eigen::MatrixXd> X_disc(X_disc_r.begin(), X_disc_r.rows(), X_disc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_disc(y_disc_r.begin(), y_disc_r.size());
	Eigen::Map<const Eigen::MatrixXd> X_conc(X_conc_r.begin(), X_conc_r.rows(), X_conc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_conc(y_conc_r.begin(), y_conc_r.size());
	Eigen::Map<const Eigen::VectorXi> group_conc(group_conc_r.begin(), group_conc_r.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());

	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);
	return wrap(-obj.hessian(params));
}

// [[Rcpp::export]]
SEXP fast_clogit_plus_glmm_cpp(
	const NumericMatrix& X_disc_r,
	const NumericVector& y_disc_r,
	const NumericMatrix& X_conc_r,
	const NumericVector& y_conc_r,
	const IntegerVector& group_conc_r,
	bool has_discordant,
	bool has_concordant,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
	bool estimate_only = false,
	double max_abs_log_sigma = 8.0,
	int maxit = 200,
	double eps_g = 1e-5,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	Eigen::Map<const Eigen::MatrixXd> X_disc(X_disc_r.begin(), X_disc_r.rows(), X_disc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_disc(y_disc_r.begin(), y_disc_r.size());
	Eigen::Map<const Eigen::MatrixXd> X_conc(X_conc_r.begin(), X_conc_r.rows(), X_conc_r.cols());
	Eigen::Map<const Eigen::VectorXd> y_conc(y_conc_r.begin(), y_conc_r.size());
	Eigen::Map<const Eigen::VectorXi> group_conc(group_conc_r.begin(), group_conc_r.size());

	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);

	int p_disc = has_discordant ? (int)X_disc.cols() : 0;
	int p_conc = has_concordant ? (int)X_conc.cols() : 0;
	int n_par = (has_discordant && has_concordant) ? p_conc + 1 : (has_concordant ? p_conc + 1 : p_disc);

	Eigen::VectorXd par = Eigen::VectorXd::Zero(n_par);
	
	if (warm_start_params.isNotNull()) {
		par = as<Eigen::VectorXd>(warm_start_params);
		if (par.size() != n_par) stop("warm_start_params size mismatch");
	} else if (warm_start_beta.isNotNull()) {
		VectorXd sb = as<VectorXd>(warm_start_beta);
		if (sb.size() == n_par) {
			par = sb;
		} else {
			// Try to fill beta part
			int n_beta = has_concordant ? p_conc : p_disc;
			if (sb.size() == n_beta) {
				par.head(n_beta) = sb;
			}
		}
	}

	FixedParamSpec fixed_spec = make_fixed_param_spec(par.size(), fixed_idx, fixed_values);

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
			Named("params") = par,
			Named("b") = par,
			Named("beta_T") = NA_REAL,
			Named("se_beta_T") = NA_REAL,
			Named("ssq_b_j") = NA_REAL,
			Named("converged") = false,
			Named("neg_loglik") = NA_REAL
		);
	}

		const int j_beta_T = has_concordant ? 1 : 0; // 0-based
		double ssq_b_j = NA_REAL;
		if (estimate_only) {
			return List::create(
				Named("params") = par,
				Named("b") = par,
				Named("beta_T") = par[j_beta_T],
				Named("se_beta_T") = NA_REAL,
				Named("ssq_b_j") = NA_REAL,
				Named("vcov") = R_NilValue,
				Named("score") = R_NilValue,
				Named("observed_information") = R_NilValue,
				Named("information") = R_NilValue,
				Named("information_type") = "observed",
				Named("hessian") = R_NilValue,
				Named("converged") = converged,
				Named("neg_loglik") = neg_ll,
				Named("neg_ll") = neg_ll,
				Named("loglik") = R_finite(neg_ll) ? -neg_ll : NA_REAL
			);
		}

		Eigen::MatrixXd info = obj.hessian(par);
		Eigen::VectorXd score(par.size());
		obj(par, score);
		score = -score;
		Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(par.size(), par.size(), NA_REAL);
		if (converged) {
		Eigen::MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(info_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(info_free.rows(), info_free.cols()));
			vcov = expand_free_covariance(par.size(), fixed_spec, inv_free, true);
			if (vcov.allFinite()) ssq_b_j = vcov(j_beta_T, j_beta_T);
		}
	}
	double se_beta_T = (std::isfinite(ssq_b_j) && ssq_b_j > 0.0) ? std::sqrt(ssq_b_j) : NA_REAL;

	return List::create(
		Named("params") = par,
		Named("b") = par,
		Named("beta_T") = par[j_beta_T],
		Named("se_beta_T") = se_beta_T,
		Named("ssq_b_j") = ssq_b_j,
		Named("vcov") = vcov,
		Named("score") = score,
		Named("observed_information") = info,
		Named("information") = info,
		Named("information_type") = "observed",
		Named("hessian") = -info,
		Named("converged") = converged,
		Named("neg_loglik") = neg_ll,
		Named("neg_ll") = neg_ll,
		Named("loglik") = R_finite(neg_ll) ? -neg_ll : NA_REAL,
		Named("fisher_information") = info
	);
}
