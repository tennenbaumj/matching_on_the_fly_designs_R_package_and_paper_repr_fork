#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <unordered_map>

using namespace Rcpp;
using namespace Eigen;

// Forward declaration from fast_logistic_regression.cpp
ModelResult fast_logistic_regression_internal(
	const Eigen::Ref<const Eigen::MatrixXd>& X,
	const Eigen::Ref<const Eigen::VectorXd>& y,
	const Eigen::Ref<const Eigen::VectorXd>& weights = Eigen::VectorXd(),
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
	bool smart_cold_start = true,
	int maxit = 100,
	double tol = 1e-8,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
	bool estimate_only = false);

namespace {

double clamp_exp_arg_hnb(double eta) {
	return std::min(eta, 700.0);
}

class TruncatedNegBinCount;

void validate_truncated_negbin_inputs(const Eigen::Ref<const MatrixXd>& X,
									  const Eigen::Ref<const VectorXd>& y,
									  const Nullable<NumericVector>& warm_start_params,
									  const Nullable<NumericMatrix>& warm_start_fisher_info,
									  const Nullable<IntegerVector>& fixed_idx,
									  const Nullable<NumericVector>& fixed_values) {
	if (X.rows() != y.size()) {
		Rcpp::stop("X and y must have compatible dimensions");
	}
	if (!X.allFinite()) {
		Rcpp::stop("X must contain only finite values");
	}
	if (!y.allFinite()) {
		Rcpp::stop("y must contain only finite values");
	}
	for (int i = 0; i < y.size(); ++i) {
		if (y[i] <= 0.0) {
			Rcpp::stop("y must contain only positive counts for truncated negative binomial regression");
		}
		const double y_round = std::round(y[i]);
		if (std::fabs(y[i] - y_round) > 1e-8) {
			Rcpp::stop("y must contain only integer-valued counts for truncated negative binomial regression");
		}
	}

	const int n_params = (int)X.cols() + 1;
	if (warm_start_params.isNotNull()) {
		NumericVector warm = NumericVector(warm_start_params);
		if (warm.size() != n_params) {
			Rcpp::stop("warm_start_params must have length ncol(X) + 1");
		}
		for (int i = 0; i < warm.size(); ++i) {
			if (!R_finite(warm[i])) {
				Rcpp::stop("warm_start_params must contain only finite values");
			}
		}
	}
	if (warm_start_fisher_info.isNotNull()) {
		NumericMatrix warm_info(warm_start_fisher_info);
		if (warm_info.nrow() != n_params || warm_info.ncol() != n_params) {
			Rcpp::stop("warm_start_fisher_info must be a square matrix with ncol(X) + 1 rows");
		}
		for (int j = 0; j < warm_info.ncol(); ++j) {
			for (int i = 0; i < warm_info.nrow(); ++i) {
				if (!R_finite(warm_info(i, j))) {
					Rcpp::stop("warm_start_fisher_info must contain only finite values");
				}
			}
		}
	}

	make_fixed_param_spec(n_params, fixed_idx, fixed_values);
}

bool validate_start_vector(const Eigen::Ref<const VectorXd>& params,
						   int expected_size,
						   std::string* reason = nullptr) {
	auto fail = [&](const std::string& msg) {
		if (reason != nullptr) *reason = msg;
		return false;
	};
	if (params.size() != expected_size) {
		return fail("start vector has incompatible length");
	}
	if (!params.allFinite()) {
		return fail("start vector must contain only finite values");
	}
	const double log_theta = params[expected_size - 1];
	if (!std::isfinite(log_theta) || log_theta < -20.0 || log_theta > 20.0) {
		return fail("start vector log-theta is outside the supported range");
	}
	return true;
}

double smart_truncated_negbin_theta_start_from_beta(const Eigen::Ref<const MatrixXd>& X,
													const Eigen::Ref<const VectorXi>& y,
													const Eigen::Ref<const VectorXd>& beta,
													double legacy_theta_start) {
	const int n = (int)X.rows();
	const int p = (int)X.cols();
	const int df = std::max(1, n - p);
	if (beta.size() != p || !beta.allFinite()) return legacy_theta_start;

	VectorXd eta = (X * beta).array().min(700.0).matrix();
	VectorXd mu = eta.array().exp().max(1e-8).matrix();
	double alpha_sum = 0.0;

	for (int i = 0; i < n; ++i) {
		const double yi = static_cast<double>(y[i]);
		const double mui = mu[i];
		alpha_sum += ((yi - mui) * (yi - mui) - mui) / (mui * mui);
	}

	const double alpha_hat = alpha_sum / static_cast<double>(df);
	if (!std::isfinite(alpha_hat) || alpha_hat <= 0.0) return legacy_theta_start;

	const double theta_hat = 1.0 / alpha_hat;
	if (!std::isfinite(theta_hat) || theta_hat <= 0.0) return legacy_theta_start;
	return std::max(0.1, theta_hat);
}

VectorXd make_truncated_negbin_start(const Eigen::Ref<const MatrixXd>& X_pos, const Eigen::Ref<const VectorXi>& y_pos) {
	const int p = (int)X_pos.cols();
	const VectorXd y_double = y_pos.cast<double>();
	const double mean_y = y_double.mean();
	const double var_y = (y_double.array() - mean_y).square().sum() /
		static_cast<double>(std::max(1, static_cast<int>(y_pos.size()) - 1));
	const double theta_start = (var_y > mean_y && mean_y > 0.0) ?
		std::max(0.1, mean_y * mean_y / (var_y - mean_y)) : 10.0;

	VectorXd legacy_beta = VectorXd::Zero(p);
	if (p > 0 && X_pos.col(0).array().isApprox(ArrayXd::Ones(X_pos.rows()), 1e-12)) {
		legacy_beta[0] = std::log(std::max(mean_y, 1e-8));
	}
	VectorXd beta_smart = ols_smart_cold_start_beta_on_log1p(X_pos, y_double);
	VectorXd beta_start = vector_is_usable_start(beta_smart, p) ? beta_smart : legacy_beta;
	VectorXd params = VectorXd::Zero(p + 1);
	params.head(p) = beta_start;
	params[p] = std::log(smart_truncated_negbin_theta_start_from_beta(X_pos, y_pos, beta_start, theta_start));
	return params;
}

std::vector<VectorXd> make_truncated_negbin_candidate_starts(const Eigen::Ref<const MatrixXd>& X_pos, const Eigen::Ref<const VectorXi>& y_pos) {
	const int p = (int)X_pos.cols();
	const VectorXd y_double = y_pos.cast<double>();
	const double mean_y = y_double.mean();
	const double var_y = (y_double.array() - mean_y).square().sum() /
		static_cast<double>(std::max(1, static_cast<int>(y_pos.size()) - 1));
	const double theta_moment = (var_y > mean_y && mean_y > 0.0) ?
		std::max(0.1, mean_y * mean_y / (var_y - mean_y)) : 10.0;

	VectorXd legacy_beta = VectorXd::Zero(p);
	if (p > 0 && X_pos.col(0).array().isApprox(ArrayXd::Ones(X_pos.rows()), 1e-12)) {
		legacy_beta[0] = std::log(std::max(mean_y, 1e-8));
	}

	std::vector<VectorXd> beta_candidates;
	beta_candidates.push_back(legacy_beta);

	VectorXd beta_log1p = ols_smart_cold_start_beta_on_log1p(X_pos, y_double);
	if (vector_is_usable_start(beta_log1p, p)) beta_candidates.push_back(beta_log1p);

	VectorXd beta_log = safe_ols_solve(X_pos, y_double.array().max(1.0).log().matrix());
	if (vector_is_usable_start(beta_log, p)) beta_candidates.push_back(beta_log);

	std::vector<VectorXd> starts;
	auto add_start = [&](const VectorXd& beta, double theta_hint) {
		if (!vector_is_usable_start(beta, p) || !std::isfinite(theta_hint) || theta_hint <= 0.0) return;
		VectorXd candidate(p + 1);
		candidate.head(p) = beta;
		candidate[p] = std::log(theta_hint);
		for (const auto& existing : starts) {
			if (existing.size() == candidate.size() && existing.isApprox(candidate, 1e-8)) return;
		}
		starts.push_back(candidate);
	};

	for (const VectorXd& beta : beta_candidates) {
		const double theta_smart = smart_truncated_negbin_theta_start_from_beta(X_pos, y_pos, beta, theta_moment);
		add_start(beta, theta_smart);
		add_start(beta, theta_moment);
		add_start(beta, 10.0);
		add_start(beta, 25.0);
	}

	if (starts.empty()) {
		starts.push_back(make_truncated_negbin_start(X_pos, y_pos));
	}
	return starts;
}

LikelihoodFitResult fit_truncated_negbin_with_fallback(TruncatedNegBinCount& fun,
														const Eigen::Ref<const VectorXd>& params,
														const FixedParamSpec& fixed_spec,
														int maxit,
														double tol,
														const std::string& optimization_alg,
														const Eigen::MatrixXd* info_start_ptr = nullptr,
														const VectorXd* fallback_start = nullptr,
														const std::vector<VectorXd>* extra_starts = nullptr,
														std::string* failure_message = nullptr) {
	LikelihoodFitResult best;
	best.params = params;
	best.value = std::numeric_limits<double>::infinity();
	std::string last_error;
	const int expected_size = (int)params.size();

	if (!validate_start_vector(params, expected_size, &last_error)) {
		Rcpp::stop(last_error);
	}

	auto maybe_keep = [&](const LikelihoodFitResult& candidate) {
		if (candidate.params.size() == 0 || !candidate.params.allFinite()) return;
		if (candidate.converged && !best.converged) {
			best = candidate;
			return;
		}
		if (candidate.converged == best.converged &&
			std::isfinite(candidate.value) &&
			(!std::isfinite(best.value) || candidate.value < best.value)) {
			best = candidate;
		}
	};

	auto try_alg = [&](const std::string& alg_try,
					   const VectorXd& start_try,
					   const Eigen::MatrixXd* warm_info_try = nullptr) -> bool {
		std::string start_reason;
		if (!validate_start_vector(start_try, expected_size, &start_reason)) {
			last_error = start_reason;
			return false;
		}
		try {
			LikelihoodFitResult fit = optimize_fixed_likelihood(
				fun,
				start_try,
				fixed_spec,
				maxit,
				tol,
				alg_try,
				"lbfgs",
				0,
				warm_info_try
			);
			maybe_keep(fit);
			return fit.converged;
		} catch (const std::exception& e) {
			last_error = e.what();
			return false;
		}
	};

	const bool has_distinct_fallback_start =
		fallback_start != nullptr &&
		fallback_start->size() == params.size() &&
		!fallback_start->isApprox(params, 1e-8);

	auto try_extra_starts = [&]() -> bool {
		if (extra_starts == nullptr) return false;
		for (const auto& start_try : *extra_starts) {
			if (start_try.size() != params.size()) continue;
			if (start_try.isApprox(params, 1e-8)) continue;
			if (has_distinct_fallback_start && start_try.isApprox(*fallback_start, 1e-8)) continue;
			if (try_alg("newton_raphson", start_try, nullptr)) return true;
			if (try_alg("lbfgs", start_try, nullptr)) return true;
		}
		return false;
	};

	if (optimization_alg == "lbfgs") {
		if (try_alg("lbfgs", params, info_start_ptr)) return best;
		VectorXd retry_start = (best.params.size() == params.size() && best.params.allFinite()) ? best.params : params.matrix();
		if (try_alg("newton_raphson", retry_start, info_start_ptr)) return best;
		if (has_distinct_fallback_start) {
			if (try_alg("newton_raphson", *fallback_start, nullptr)) return best;
			if (try_alg("lbfgs", *fallback_start, nullptr)) return best;
		}
		if (try_extra_starts()) return best;
	} else {
		if (try_alg("newton_raphson", params, info_start_ptr)) return best;
		if (has_distinct_fallback_start && try_alg("newton_raphson", *fallback_start, nullptr)) return best;
		if (try_extra_starts()) return best;
	}

	if (failure_message != nullptr) {
		*failure_message = last_error;
	}
	return best;
}

class TruncatedNegBinCount {
private:
	const MatrixXd m_X;
	const VectorXi m_y;
	const int m_n;
	const int m_p;

	// Distinct-y precomputation (built once at construction)
	std::vector<int>    m_y_slot;
	std::vector<double> m_distinct_y;
	std::vector<double> m_lgamma_y1;        // lgamma(y+1) per distinct y — constant
	std::vector<double> m_lgamma_yptheta;   // preallocated; filled per operator() call
	std::vector<double> m_digamma_yptheta;  // preallocated; filled per operator() and hessian()
	std::vector<double> m_trigamma_yptheta; // preallocated; filled per hessian() call

public:
	TruncatedNegBinCount(const MatrixXd& X, const VectorXi& y)
		: m_X(X), m_y(y), m_n((int)X.rows()), m_p((int)X.cols()), m_y_slot(X.rows(), -1)
	{
		std::unordered_map<int, int> seen;
		for (int i = 0; i < m_n; ++i) {
			const int yi = m_y[i];
			auto it = seen.find(yi);
			if (it == seen.end()) {
				const int slot = (int)m_distinct_y.size();
				seen[yi] = slot;
				m_distinct_y.push_back((double)yi);
				m_lgamma_y1.push_back(std::lgamma((double)yi + 1.0));
				m_y_slot[i] = slot;
			} else {
				m_y_slot[i] = it->second;
			}
		}
		const int nd = (int)m_distinct_y.size();
		m_lgamma_yptheta.resize(nd);
		m_digamma_yptheta.resize(nd);
		m_trigamma_yptheta.resize(nd);
	}

	double operator()(const Eigen::Ref<const VectorXd>& params, Eigen::Ref<VectorXd> grad) {
		const VectorXd beta = params.head(m_p);
		const double theta = std::exp(params[m_p]);
		const double log_r = std::log(theta);
		const double lgamma_r = std::lgamma(theta);
		const double digamma_r = fast_digamma(theta);

		VectorXd eta = (m_X * beta).array().min(700.0).matrix();
		VectorXd mu = eta.array().exp().max(1e-10).matrix();

		// Fill per-call distinct-y tables
		const int nd = (int)m_distinct_y.size();
		for (int k = 0; k < nd; ++k) {
			const double ypt = m_distinct_y[k] + theta;
			m_lgamma_yptheta[k]  = std::lgamma(ypt);
			m_digamma_yptheta[k] = fast_digamma(ypt);
		}

		double neg_ll = 0.0;
		VectorXd eta_score = VectorXd::Zero(m_n);
		double score_log_theta = 0.0;

		for (int i = 0; i < m_n; ++i) {
			const double mu_i = mu[i];
			const double yi = m_distinct_y[m_y_slot[i]];
			const int slot = m_y_slot[i];
			const double denom = theta + mu_i;
			const double log_denom = std::log(denom);

			const double log_p0 = theta * (log_r - log_denom);
			double p0 = std::exp(log_p0);
			p0 = std::min(std::max(p0, 1e-12), 1.0 - 1e-12);
			const double trunc_denom = 1.0 - p0;
			// log(1-p0) ≈ -p0 when p0 < 1e-7 (error < 5e-15); avoids std::log
			const double log_trunc_denom = (p0 < 1e-7) ? -p0 : std::log(trunc_denom);

			// Explicit NB log-PMF — replaces R::dnbinom_mu overhead
			const double log_mu_i = std::log(mu_i);
			neg_ll -= m_lgamma_yptheta[slot] - lgamma_r - m_lgamma_y1[slot]
			        + theta * (log_r - log_denom)
			        + yi * (log_mu_i - log_denom)
			        - log_trunc_denom;

			const double standard_eta_score = yi - mu_i * (yi + theta) / denom;
			const double trunc_eta_corr = -mu_i * p0 * theta / (trunc_denom * denom);
			eta_score[i] = standard_eta_score + trunc_eta_corr;

			const double dlogf_dr =
				m_digamma_yptheta[slot] - digamma_r +
				log_r - log_denom +
				1.0 - (yi + theta) / denom;
			const double dlogp0_dr = log_r - log_denom + 1.0 - theta / denom;
			score_log_theta += theta * (dlogf_dr + (p0 * dlogp0_dr) / trunc_denom);
		}

		grad.resize(m_p + 1);
		grad.head(m_p) = -(m_X.transpose() * eta_score);
		grad[m_p] = -score_log_theta;

		return neg_ll;
	}

	MatrixXd hessian(const Eigen::Ref<const VectorXd>& params) {
		const int total_p = m_p + 1;
		MatrixXd H = MatrixXd::Zero(total_p, total_p);
		const VectorXd beta = params.head(m_p);
		const double r = std::exp(params[m_p]);
		const double log_r = std::log(r);
		const double digamma_r = fast_digamma(r);
		const double trigamma_r = R::trigamma(r);
		VectorXd eta = (m_X * beta).array().min(700.0).matrix();
		VectorXd beta_weights = VectorXd::Zero(m_n);
		VectorXd cross_weights = VectorXd::Zero(m_n);

		// Fill per-call distinct-y tables for hessian
		const int nd = (int)m_distinct_y.size();
		for (int k = 0; k < nd; ++k) {
			const double ypt = m_distinct_y[k] + r;
			m_digamma_yptheta[k]  = fast_digamma(ypt);
			m_trigamma_yptheta[k] = R::trigamma(ypt);
		}

		for (int i = 0; i < m_n; ++i) {
			const double mu_i = std::max(std::exp(eta[i]), 1e-10);
			const double yi = m_distinct_y[m_y_slot[i]];
			const int slot = m_y_slot[i];
			const double denom = r + mu_i;
			const double denom_sq = denom * denom;
			const double log_denom = std::log(denom);

			double log_p0 = r * (log_r - log_denom);
			double p0 = std::exp(log_p0);
			p0 = std::min(std::max(p0, 1e-12), 1.0 - 1e-12);
			const double q0 = 1.0 - p0;

			const double c_eta = r * mu_i / denom;
			const double dlogp0_dr = log_r - log_denom + 1.0 - r / denom;

			const double d_score_eta_d_eta =
				- (yi + r) * r * mu_i / denom_sq
				- r * r * mu_i * p0 / (denom_sq * q0)
				+ c_eta * c_eta * p0 / (q0 * q0);
			beta_weights[i] = -d_score_eta_d_eta;

			const double d_score_eta_d_log_r =
				r * mu_i * (yi - mu_i) / denom_sq
				- r * mu_i * mu_i * p0 / (denom_sq * q0)
				- r * r * mu_i * p0 * dlogp0_dr / (denom * q0 * q0);
			cross_weights[i] = -d_score_eta_d_log_r;

			const double dlogf_dr =
				m_digamma_yptheta[slot] - digamma_r +
				log_r - log_denom +
				1.0 - (yi + r) / denom;
			const double d2logf_dr2 =
				m_trigamma_yptheta[slot] - trigamma_r +
				1.0 / r - 1.0 / denom + (yi - mu_i) / denom_sq;
			const double d2logp0_dr2 = 1.0 / r - 1.0 / denom - mu_i / denom_sq;
			const double p0_over_q0 = p0 / q0;
			const double d_score_log_r_d_log_r =
				r * (dlogf_dr + p0_over_q0 * dlogp0_dr) +
				r * r * (
					d2logf_dr2 +
					p0 * dlogp0_dr * dlogp0_dr / (q0 * q0) +
					p0_over_q0 * d2logp0_dr2
				);
			H(m_p, m_p) -= d_score_log_r_d_log_r;
		}

		H.topLeftCorner(m_p, m_p).noalias() = weighted_crossprod(m_X, beta_weights);
		H.topRightCorner(m_p, 1).noalias() = m_X.transpose() * cross_weights;
		H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
		return H;
	}

	MatrixXd expected_hessian(const Eigen::Ref<const VectorXd>& params) {
		const int total_p = m_p + 1;
		MatrixXd H = MatrixXd::Zero(total_p, total_p);
		const VectorXd beta = params.head(m_p);
		const double r = std::exp(params[m_p]);
		VectorXd eta = (m_X * beta).array().min(700.0).matrix();

		for (int i = 0; i < m_n; ++i) {
			const double mu_i = std::max(std::exp(eta[i]), 1e-10);
			const double p0 = std::pow(r / (r + mu_i), r);
			const double q_pos = std::max(1.0 - p0, 1e-12);
			const double sd = std::sqrt(mu_i + mu_i * mu_i / r);
			const int soft_min_y = std::max(20, static_cast<int>(std::ceil(mu_i + 12.0 * sd)));
			const int max_y = std::max(soft_min_y, 10000);
			double cum = 0.0;

			MatrixXd Xi(1, m_p);
			Xi.row(0) = m_X.row(i);
			for (int yv = 1; yv <= max_y; ++yv) {
				const double prob = R::dnbinom_mu(static_cast<double>(yv), r, mu_i, false) / q_pos;
				if (prob > 0.0 && std::isfinite(prob)) {
					VectorXi yi(1);
					yi[0] = yv;
					TruncatedNegBinCount one_row(Xi, yi);
					H.noalias() += prob * one_row.hessian(params);
					cum += prob;
				}
				if (yv >= soft_min_y && 1.0 - cum < 1e-10) break;
			}
		}
		return (H + H.transpose()) / 2.0;
	}
};

}

static List build_positive_hurdle_negbin_data(const Eigen::Ref<const MatrixXd>& X,
											  const Eigen::Ref<const VectorXd>& y) {
	std::vector<int> pos_rows;
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_rows.push_back(i);
	}

	const int p = (int)X.cols();
	MatrixXd X_pos((int)pos_rows.size(), p);
	VectorXi y_pos((int)pos_rows.size());
	for (size_t k = 0; k < pos_rows.size(); ++k) {
		const int i = pos_rows[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = static_cast<int>(y[i]);
	}
	return List::create(Named("X_pos") = X_pos, Named("y_pos") = y_pos);
}

// [[Rcpp::export]]
Eigen::VectorXd get_hurdle_negbin_count_score_cpp(SEXP X_r,
												  SEXP y_r,
												  SEXP params_sexp) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	List pos = build_positive_hurdle_negbin_data(X, y);
	MatrixXd X_pos = pos["X_pos"];
	VectorXi y_pos = pos["y_pos"];
	TruncatedNegBinCount fun(X_pos, y_pos);
	VectorXd grad(params.size());
	fun(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_hurdle_negbin_count_hessian_cpp(SEXP X_r,
													SEXP y_r,
													SEXP params_sexp) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	NumericVector params_r(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXd> params(params_r.begin(), params_r.size());
	List pos = build_positive_hurdle_negbin_data(X, y);
	MatrixXd X_pos = pos["X_pos"];
	VectorXi y_pos = pos["y_pos"];
	TruncatedNegBinCount fun(X_pos, y_pos);
	return -fun.hessian(params);
}

// [[Rcpp::export]]
List fast_hurdle_negbin_cpp(SEXP X_r,
						   SEXP y_r,
						   SEXP X_hurdle_r,
						   Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
						   bool smart_cold_start = true,
						   int maxit = 1000,
						   double tol = 1e-8,
						   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
						   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
						   std::string optimization_alg = "lbfgs",
						   Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
						   Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_hurdle_fisher_info = R_NilValue,
						   bool estimate_only = false) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	NumericMatrix X_hurdle_mat(X_hurdle_r);
	Eigen::Map<const Eigen::MatrixXd> X_hurdle(X_hurdle_mat.begin(), X_hurdle_mat.nrow(), X_hurdle_mat.ncol());

	const int n = (int)X.rows();
	const int p = (int)X.cols();
	const int p_hurdle = (int)X_hurdle.cols();
	std::string alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", false);
	FixedParamSpec count_fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

	VectorXd y_pos_ind = (y.array() > 0.0).cast<double>();
	Eigen::VectorXd hurdle_b = Eigen::VectorXd::Constant(p_hurdle, NA_REAL);
	Eigen::MatrixXd hurdle_XtWX;
	bool hurdle_converged = false;

	if (y_pos_ind.minCoeff() < y_pos_ind.maxCoeff()) {
		ModelResult hurdle_res = fast_logistic_regression_internal(X_hurdle, y_pos_ind, Eigen::VectorXd(), R_NilValue, true, 100, 1e-8, R_NilValue, R_NilValue, alg, R_NilValue, warm_start_hurdle_fisher_info, estimate_only);
		hurdle_b = hurdle_res.b;
		hurdle_XtWX = hurdle_res.XtWX;
		hurdle_converged = hurdle_res.converged;
	}

	std::vector<int> pos_rows;
	for (int i = 0; i < n; ++i) {
		if (y[i] > 0.0) pos_rows.push_back(i);
	}

	if (static_cast<int>(pos_rows.size()) <= p) {
		return List::create(
			Named("b") = NumericVector(p, NA_REAL),
			Named("theta_hat") = NA_REAL,
			Named("converged") = false,
			Named("hurdle_b") = hurdle_b,
			Named("hurdle_converged") = hurdle_converged
		);
	}

	MatrixXd X_pos((int)pos_rows.size(), p);
	VectorXi y_pos((int)pos_rows.size());
	for (size_t k = 0; k < pos_rows.size(); ++k) {
		const int i = pos_rows[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = static_cast<int>(y[i]);
	}

	std::vector<VectorXd> start_candidates = make_truncated_negbin_candidate_starts(X_pos, y_pos);
	VectorXd params = start_candidates.front();
	if (warm_start_params.isNotNull()) {
		params = as<VectorXd>(NumericVector(warm_start_params));
	} else if (!smart_cold_start) {
		params.setZero();
		params[p] = std::log(1.0);
	}

	TruncatedNegBinCount fun(X_pos, y_pos);
	double neg_ll = NA_REAL;
	bool converged = false;

	Eigen::MatrixXd info_start;
	const Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	std::string failure_message;
	LikelihoodFitResult fit = fit_truncated_negbin_with_fallback(
		fun,
		params,
		count_fixed_spec,
		maxit,
		tol,
		alg,
		info_start_ptr,
		nullptr,
		&start_candidates,
		&failure_message
	);
	params = fit.params;
	neg_ll = fit.value;
	converged = fit.converged;
	if (!converged && !failure_message.empty()) {
		Rcpp::warning(failure_message);
	}

	VectorXd beta = params.head(p);
	double theta_hat = std::exp(params[p]);

	if (estimate_only) {
		return List::create(
			Named("b") = beta,
			Named("theta_hat") = theta_hat,
			Named("converged") = converged,
			Named("hurdle_b") = hurdle_b,
			Named("hurdle_converged") = hurdle_converged,
			Named("neg_ll") = neg_ll
		);
	}

	MatrixXd observed_information = fun.hessian(params);

	return List::create(
		Named("b") = beta,
		Named("theta_hat") = theta_hat,
		Named("converged") = converged,
		Named("hurdle_b") = hurdle_b,
		Named("hurdle_converged") = hurdle_converged,
		Named("observed_information") = observed_information,
		Named("fisher_information") = observed_information,
		Named("information") = observed_information,
		Named("information_type") = "observed",
		Named("hessian") = -observed_information,
		Named("hurdle_fisher_information") = hurdle_XtWX
	);
}

//' @title Fast Hurdle Negative Binomial Regression with Variance (C++)
//' @description Hurdle NB regression with full variance-covariance matrix.
//' @param X Matrix of predictors for the count component.
//' @param y Vector of responses.
//' @param X_hurdle Matrix of predictors for the hurdle component.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param warm_start_params Optional starting values for count parameters. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first iteration.
//' @param warm_start_hurdle_fisher_info Optional initial Fisher Information matrix for the hurdle iteration.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_hurdle_negbin_with_var_cpp(SEXP X_r,
									 SEXP y_r,
									 SEXP X_hurdle_r,
									 int j = 2,
									 Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
									 bool smart_cold_start = true,
									 int maxit = 1000,
									 double tol = 1e-8,
									 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
									 std::string optimization_alg = "lbfgs",
									 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_hurdle_fisher_info = R_NilValue) {
	NumericMatrix X_mat(X_r);
	NumericVector y_vec(y_r);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	NumericMatrix X_hurdle_mat(X_hurdle_r);
	Eigen::Map<const Eigen::MatrixXd> X_hurdle(X_hurdle_mat.begin(), X_hurdle_mat.nrow(), X_hurdle_mat.ncol());

	List fit = fast_hurdle_negbin_cpp(X_r, y_r, X_hurdle_r, warm_start_params, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, warm_start_hurdle_fisher_info);
	SEXP b_sexp = fit["b"];
	NumericVector b_nv(b_sexp);
	const int p = (int)b_nv.size();

	double hurdle_ssq_b_j = NA_REAL;
	double hurdle_ssq_b_2 = NA_REAL;
	if (fit.containsElementNamed("hurdle_b")) {
		VectorXd y_pos_ind = (y.array() > 0.0).cast<double>();
		if (y_pos_ind.minCoeff() < y_pos_ind.maxCoeff()) {
			MatrixXd hurdle_XtWX;
			bool have_hurdle_info = false;
			if (fit.containsElementNamed("hurdle_fisher_information")) {
				hurdle_XtWX = as<MatrixXd>(fit["hurdle_fisher_information"]);
				have_hurdle_info = hurdle_XtWX.allFinite() && hurdle_XtWX.rows() == X_hurdle.cols();
			}
			if (!have_hurdle_info) {
				ModelResult hurdle_res = fast_logistic_regression_internal(X_hurdle, y_pos_ind, Eigen::VectorXd(), R_NilValue, true, 100, 1e-8, R_NilValue, R_NilValue, optimization_alg, R_NilValue, R_NilValue);
				hurdle_XtWX = hurdle_res.XtWX;
			}
			FixedParamSpec hurdle_spec = make_fixed_param_spec((int)X_hurdle.cols(), R_NilValue, R_NilValue);
			MatrixXd info_free = subset_matrix(hurdle_XtWX, hurdle_spec.free_idx, hurdle_spec.free_idx);
			auto hurdle_free_idx_of = [&](int k) -> int {
				for (int jj = 0; jj < (int)hurdle_spec.free_idx.size(); ++jj)
					if (hurdle_spec.free_idx[jj] == k) return jj + 1;
				return -1;
			};
			int hfree_j = (j > 0 && j <= X_hurdle.cols()) ? hurdle_free_idx_of(j - 1) : -1;
			if (hfree_j > 0) hurdle_ssq_b_j = compute_diagonal_inverse_entry(info_free, hfree_j);
			int hfree_2 = (X_hurdle.cols() >= 2) ? hurdle_free_idx_of(1) : -1;
			if (hfree_2 > 0) hurdle_ssq_b_2 = compute_diagonal_inverse_entry(info_free, hfree_2);
		}
	}

	double ssq_b_j = NA_REAL;
	double ssq_b_2 = NA_REAL;
	if (p > 0 && fit.containsElementNamed("theta_hat") && fit.containsElementNamed("observed_information")) {
		double theta_hat = as<double>(fit["theta_hat"]);
		if (R_finite(theta_hat) && p >= j) {
			MatrixXd H = as<MatrixXd>(fit["observed_information"]);
			if (H.allFinite() && H.rows() == p + 1) {
				FixedParamSpec count_fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
				MatrixXd H_free = subset_matrix(H, count_fixed_spec.free_idx, count_fixed_spec.free_idx);
				auto cnt_free_idx_of = [&](int k) -> int {
					for (int jj = 0; jj < (int)count_fixed_spec.free_idx.size(); ++jj)
						if (count_fixed_spec.free_idx[jj] == k) return jj + 1;
					return -1;
				};
				int cfree_j = (j > 0 && j <= p + 1) ? cnt_free_idx_of(j - 1) : -1;
				if (cfree_j > 0) ssq_b_j = compute_diagonal_inverse_entry(H_free, cfree_j);
				int cfree_2 = (p >= 2) ? cnt_free_idx_of(1) : -1;
				if (cfree_2 > 0) ssq_b_2 = compute_diagonal_inverse_entry(H_free, cfree_2);
			}
		}
	}

	fit["ssq_b_j"] = ssq_b_j;
	fit["ssq_b_2"] = ssq_b_2;
	fit["hurdle_ssq_b_j"] = hurdle_ssq_b_j;
	fit["hurdle_ssq_b_2"] = hurdle_ssq_b_2;
	return fit;
}

// [[Rcpp::export]]
List fast_truncated_negbin_count_cpp(SEXP X_r,
									 SEXP y_r,
									 Nullable<NumericVector> warm_start_params = R_NilValue,
									 bool smart_cold_start = true,
									 bool estimate_only = false,
									 int maxit = 1000,
									 double tol = 1e-8,
									 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
									 std::string optimization_alg = "lbfgs",
									 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());

    optimization_alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", false);
	validate_truncated_negbin_inputs(X, y, warm_start_params, warm_start_fisher_info, fixed_idx, fixed_values);

    List pos = build_positive_hurdle_negbin_data(X, y);
    MatrixXd X_pos = pos["X_pos"];
    VectorXi y_pos = pos["y_pos"];
    const int p = (int)X_pos.cols();
    if (X_pos.rows() <= p){
            return List::create(
                    Named("b") = NumericVector(p, NA_REAL),
                    Named("params") = NumericVector(p + 1, NA_REAL),
                    Named("converged") = false,
                    Named("neg_ll") = NA_REAL
            );
    }

    std::vector<VectorXd> start_candidates = make_truncated_negbin_candidate_starts(X_pos, y_pos);
    VectorXd heuristic_start = start_candidates.front();
    VectorXd params = heuristic_start;
    if (warm_start_params.isNotNull()) {
            params = as<VectorXd>(NumericVector(warm_start_params));
    } else if (!smart_cold_start) {
            params.setZero();
            params[p] = std::log(1.0);
    }
	std::string start_error;
	if (!validate_start_vector(params, p + 1, &start_error)) {
		Rcpp::stop(start_error);
	}

    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
    TruncatedNegBinCount fun(X_pos, y_pos);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    std::string failure_message;
    LikelihoodFitResult fit = fit_truncated_negbin_with_fallback(
            fun,
            params,
            fixed_spec,
            maxit,
            tol,
            optimization_alg,
            info_start_ptr,
            &heuristic_start,
            &start_candidates,
            &failure_message
    );
    if (!fit.converged && fit.params.size() == params.size()) {
            if (!fit.params.allFinite()) {
                    fit.params = params;
            }
    }
    if (!fit.converged && !failure_message.empty()) {
            Rcpp::warning(failure_message);
    }
    if (fit.params.size() != p + 1 || !fit.params.allFinite()) {
            return List::create(
                    Named("b") = NumericVector(p, NA_REAL),
                    Named("params") = NumericVector(p + 1, NA_REAL),
                    Named("converged") = false,
                    Named("neg_ll") = NA_REAL
            );
    }

    params = fit.params;
    VectorXd beta = params.head(p);
    if (estimate_only) {
            return List::create(
                    Named("b") = beta,
                    Named("params") = params,
                    Named("converged") = fit.converged,
                    Named("neg_ll") = fit.value
            );
    }

    MatrixXd observed_information = fun.hessian(params);
    List out = List::create(
            Named("b") = beta,
            Named("params") = params,
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.value,
            Named("observed_information") = observed_information,
            Named("fisher_information") = observed_information,
            Named("information") = observed_information,
            Named("information_type") = "observed",
            Named("hessian") = -observed_information
    );
    if (!fit.converged){
            return out;
    }

    VectorXd score = get_hurdle_negbin_count_score_cpp(X_r, y_r, Rcpp::wrap(params));
    out["score"] = score;
    return out;
}
