// Ordinal Cumulative-Logit GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  logit P(Y_ij <= k | u_i) = alpha_k - X_ij' beta - u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij \in {1, ..., K}

#include "_helper_functions.h"
#include "_glmm_engine.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

namespace {

// Maps ordinal y (1..K) to indicators I(y <= k) for k=1..K-1
struct OrdinalGLMMData {
	const Eigen::Ref<const Eigen::MatrixXd> X;
	const std::vector<int>& y;
	const std::vector<int>& group_id;
	int K;
	int n_gh;
	double max_abs_log_sigma;

	int n;
	int p;
	int m_n_groups;
	std::vector<int> m_group_start;
	std::vector<int> m_group_end;
	glmm::GHRule m_gh;

	OrdinalGLMMData(const Eigen::Ref<const Eigen::MatrixXd>& X_in, const std::vector<int>& y_in, 
	                const std::vector<int>& gid_in, int K_in, int n_gh_in, double max_als)
		: X(X_in), y(y_in), group_id(gid_in), K(K_in), n_gh(n_gh_in), max_abs_log_sigma(max_als), m_gh(glmm::gauss_hermite_rule(n_gh_in)) {
		n = (int)y.size();
		p = (int)X.cols();
		
		// Map group_id to contiguous blocks (assuming already sorted by group_id)
		m_n_groups = 0;
		if (n > 0) {
			m_n_groups = 1;
			m_group_start.push_back(0);
			for (int i = 1; i < n; ++i) {
				if (group_id[i] != group_id[i - 1]) {
					m_group_end.push_back(i);
					m_group_start.push_back(i);
					m_n_groups++;
				}
			}
			m_group_end.push_back(n);
		}
	}
};

class OrdinalGLMMObjective {
	const OrdinalGLMMData& m_dat;
	VectorXd m_alpha;
	VectorXd m_eta_fixed;
	VectorXd m_log_lik_k;
	MatrixXd m_grad_k;
	VectorXd m_p_plus;
	VectorXd m_p_minus;
	VectorXd m_g_plus;
	VectorXd m_g_minus;
public:
	OrdinalGLMMObjective(const OrdinalGLMMData& dat)
		: m_dat(dat),
		  m_alpha(dat.K - 1),
		  m_eta_fixed(dat.n),
		  m_log_lik_k(dat.n_gh),
		  m_grad_k(dat.K - 1 + dat.p + 1, dat.n_gh),
		  m_p_plus(dat.K - 1 + dat.p + 1),
		  m_p_minus(dat.K - 1 + dat.p + 1),
		  m_g_plus(dat.K - 1 + dat.p + 1),
		  m_g_minus(dat.K - 1 + dat.p + 1) {}

	// par = [alpha_1, log_diff_alpha_2, ..., log_diff_alpha_{K-1}, beta_1, ..., beta_p, log_sigma]
	double operator()(const Eigen::Ref<const VectorXd>& par, Eigen::Ref<VectorXd> grad) {
		const int n_alpha = m_dat.K - 1;
		const int p = m_dat.p;
		const int total = n_alpha + p + 1;

		// Recover actual cutpoints alpha
		m_alpha[0] = par[0];
		for (int k = 1; k < n_alpha; ++k) m_alpha[k] = m_alpha[k - 1] + std::exp(par[k]);

		const auto beta = par.segment(n_alpha, p);
		const double log_sigma = std::clamp(par[total - 1], -m_dat.max_abs_log_sigma, m_dat.max_abs_log_sigma);
		const double sigma = std::exp(log_sigma);

		m_eta_fixed.noalias() = m_dat.X * beta;
		const VectorXd& nodes = m_dat.m_gh.nodes;
		const VectorXd& log_weights = m_dat.m_gh.log_norm_weights;

		double total_neg_ll = 0.0;
		grad.setZero();

		for (int gi = 0; gi < m_dat.m_n_groups; ++gi) {
			const int start = m_dat.m_group_start[gi];
			const int end   = m_dat.m_group_end[gi];
			
			m_grad_k.setZero();

			for (int k = 0; k < m_dat.n_gh; ++k) {
				const double u = std::sqrt(2.0) * sigma * nodes[k];
				double ll_k = 0.0;
				for (int i = start; i < end; ++i) {
					const double eta_ij = m_eta_fixed[i] + u;
					const int y_ij = m_dat.y[i];
					
					double prob_ij;
					if (y_ij == 1) {
						prob_ij = plogis_safe(m_alpha[0] - eta_ij);
					} else if (y_ij == m_dat.K) {
						prob_ij = 1.0 - plogis_safe(m_alpha[n_alpha - 1] - eta_ij);
					} else {
						prob_ij = plogis_safe(m_alpha[y_ij - 1] - eta_ij) - plogis_safe(m_alpha[y_ij - 2] - eta_ij);
					}
					prob_ij = std::max(prob_ij, 1e-15);
					ll_k += std::log(prob_ij);

					// Gradient w.r.t. beta
					double dprob_deta = 0.0;
					if (y_ij == 1) {
						dprob_deta = -dplogis_safe(m_alpha[0] - eta_ij);
					} else if (y_ij == m_dat.K) {
						dprob_deta = dplogis_safe(m_alpha[n_alpha - 1] - eta_ij);
					} else {
						dprob_deta = dplogis_safe(m_alpha[y_ij - 2] - eta_ij) - dplogis_safe(m_alpha[y_ij - 1] - eta_ij);
					}
					m_grad_k.col(k).segment(n_alpha, p) -= (dprob_deta / prob_ij) * m_dat.X.row(i).transpose();
					
					// Gradient w.r.t. log_sigma (via u)
					m_grad_k(total - 1, k) -= (dprob_deta / prob_ij) * u;

					// Gradient w.r.t. alpha params
					if (y_ij == 1) {
						m_grad_k(0, k) -= dplogis_safe(m_alpha[0] - eta_ij) / prob_ij;
					} else if (y_ij == m_dat.K) {
						double d_alpha_Km1 = dplogis_safe(m_alpha[n_alpha - 1] - eta_ij);
						for (int j = 0; j < n_alpha; ++j) {
							double d_cut = (j == 0) ? 1.0 : std::exp(par[j]);
							m_grad_k(j, k) += (d_alpha_Km1 / prob_ij) * d_cut;
						}
					} else {
						// d(plogis(alpha_{y-1} - eta) - plogis(alpha_{y-2} - eta))
						double d1 = dplogis_safe(m_alpha[y_ij - 1] - eta_ij);
						double d2 = dplogis_safe(m_alpha[y_ij - 2] - eta_ij);
						for (int j = 0; j < n_alpha; ++j) {
							double d_cut = (j == 0) ? 1.0 : std::exp(par[j]);
							if (j <= y_ij - 1) m_grad_k(j, k) -= (d1 / prob_ij) * d_cut;
							if (j <= y_ij - 2) m_grad_k(j, k) += (d2 / prob_ij) * d_cut;
						}
					}
				}
				m_log_lik_k[k] = ll_k + log_weights[k];
			}
			
			double max_ll = m_log_lik_k.maxCoeff();
			double sum_exp = 0.0;
			for (int k = 0; k < m_dat.n_gh; ++k) sum_exp += std::exp(m_log_lik_k[k] - max_ll);
			double ll_gi = max_ll + std::log(sum_exp);
			total_neg_ll -= ll_gi;

			for (int k = 0; k < m_dat.n_gh; ++k) {
				double pk = std::exp(m_log_lik_k[k] - ll_gi);
				grad += pk * m_grad_k.col(k);
			}
		}
		return total_neg_ll;
	}

	MatrixXd hessian(const Eigen::Ref<const VectorXd>& par) {
		const int total = (int)par.size();
		MatrixXd H = MatrixXd::Zero(total, total);
		double h = 1e-4;
		for (int j = 0; j < total; ++j) {
			m_p_plus = par; m_p_plus[j] += h;
			(*this)(m_p_plus, m_g_plus);
			m_p_minus = par; m_p_minus[j] -= h;
			(*this)(m_p_minus, m_g_minus);
			H.col(j) = (m_g_plus - m_g_minus) / (2.0 * h);
		}
		return (H + H.transpose()) / 2.0;
	}
};

VectorXd ols_start_beta(const Eigen::Ref<const MatrixXd>& X, const Eigen::Ref<const VectorXd>& y) {
	return (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

} // namespace

//' @title Fast Ordinal GLMM (C++)
//' @description High-performance ordinal cumulative-logit GLMM fitting using Gauss-Hermite quadrature and L-BFGS.
//' @param X A numeric matrix of predictors (no intercept).
//' @param y A numeric vector of ordinal responses (1, 2, ...).
//' @param group_id A numeric vector of group identifiers.
//' @param K Number of ordinal levels.
//' @param j_T 0-based index of the treatment effect in the beta vector.
//' @param warm_start_params Optional starting values for all parameters [alpha, beta, log_sigma]. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_beta Optional starting values for coefficients. If provided (and \code{warm_start_params} is not), \code{smart_cold_start} is ignored. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param estimate_only If TRUE, skip variance component calculations.
//' @param n_gh Number of Gauss-Hermite nodes.
//' @param max_abs_log_sigma Maximum allowed value for log(sigma).
//' @param maxit Maximum number of iterations.
//' @param eps_g Convergence tolerance.
//' @param optimization_alg Optimization algorithm.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first iteration.
//' @return A list containing coefficients, thresholds, log_sigma, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_glmm_cpp(
	const Rcpp::NumericMatrix& X,     // n x p, NO intercept column; treatment at column j_T (0-based)
	const Rcpp::IntegerVector& y,     // 1-indexed ordinal outcomes, length n
	const Rcpp::IntegerVector& group_id, // group IDs, length n (sorted internally)
	int K,                        // number of ordinal levels
	int j_T,                      // 0-based treatment column index in X
	bool smart_cold_start = true,
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXi> map_y(y.begin(), y.size());
    Eigen::Map<const Eigen::VectorXi> map_group_id(group_id.begin(), group_id.size());

	const int n = map_X.rows();
	const int p = map_X.cols();
	const int n_alpha = K - 1;
	const int total = n_alpha + p + 1; // cutpoint params + betas + log_sigma
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	// Convert Eigen/R vectors to std::vector for OrdinalGLMMData
	std::vector<int> y_v(n), gid_v(n);
	for (int i = 0; i < n; ++i) { y_v[i] = map_y[i]; gid_v[i] = map_group_id[i]; }

	OrdinalGLMMData dat(map_X, y_v, gid_v, K, n_gh, max_abs_log_sigma);

	// Initialize parameters
	Eigen::VectorXd par(total);
	if (warm_start_params.isNotNull()) {
		Rcpp::NumericVector sv(warm_start_params);
		if (sv.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sv[i];
		} else {
			// Fallback: zero-initialize
			par.head(n_alpha).setZero();
			par.segment(n_alpha, p).setZero();
			par[total - 1] = -3.0;
		}
	} else if (warm_start_beta.isNotNull()) {
		VectorXd sb = as<VectorXd>(warm_start_beta);
		if (sb.size() == total) {
			par = sb;
		} else if (sb.size() == p) {
			par.head(n_alpha).setZero();
			par.segment(n_alpha, p) = sb;
			par[total - 1] = -3.0;
		} else {
            par.head(n_alpha).setZero();
			par.segment(n_alpha, p).setZero();
			par[total - 1] = -3.0;
        }
	} else if (smart_cold_start) {
		// Cutpoints: alpha_1 = 0, log_diffs = 0 (evenly spaced by 1)
		par.head(n_alpha).setZero();
		// Betas: OLS on y (rough)
		Eigen::VectorXd y_double = map_y.cast<double>();
		par.segment(n_alpha, p) = ols_start_beta(map_X, y_double);
		par[total - 1] = -3.0;
	} else {
		par.head(n_alpha).setZero();
		par.segment(n_alpha, p).setZero();
		par[total - 1] = -3.0;
	}

	OrdinalGLMMObjective obj(dat);

	Eigen::MatrixXd info_start;
	const Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par = fit.params;
		neg_ll = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("b")          = par.segment(n_alpha, p),
			Named("alpha")      = par.head(n_alpha),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const int j_T_full = n_alpha + j_T;
	Eigen::MatrixXd information = obj.hessian(par);
	double ssq_b_T = NA_REAL;
	if (!estimate_only && converged) {
		Eigen::MatrixXd H_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
			if (inv_free.allFinite()) {
				Eigen::MatrixXd inv = expand_free_covariance(total, fixed_spec, inv_free, true);
				if (j_T_full < total) ssq_b_T = inv(j_T_full, j_T_full);
			}
		}
	}

	return List::create(
		Named("b")          = par.segment(n_alpha, p),
		Named("alpha")      = par.head(n_alpha),
		Named("params")     = par,
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("converged")  = converged,
		Named("neg_loglik") = neg_ll,
		Named("fisher_information") = information
	);
}
