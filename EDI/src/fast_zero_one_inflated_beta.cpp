#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

struct DigammaFunctor {
	double operator()(double x) const {
		return fast_digamma(x);
	}
};

struct LgammaFunctor {
	double operator()(double x) const {
		return std::lgamma(x);
	}
};

inline Eigen::VectorXd logistic(const Eigen::VectorXd& x){
	return (1.0 / (1.0 + (-x).array().exp())).matrix();
}

inline void clamp_probs(Eigen::VectorXd& v, double lower = 1e-8, double upper = 1.0 - 1e-8){
	for (int i = 0; i < v.size(); ++i){
		if (v[i] < lower){
			v[i] = lower;
		} else if (v[i] > upper){
			v[i] = upper;
		}
	}
}

class ZeroOneInflatedBeta {
public:
	ZeroOneInflatedBeta(const Eigen::Ref<const Eigen::VectorXd>& y, 
	                   const Eigen::Ref<const Eigen::MatrixXd>& X, 
	                   const Eigen::Ref<const Eigen::MatrixXd>& X_zero_one):
		m_y(y),
		m_X(X),
		m_X_zero_one(X_zero_one),
		m_n(X.rows()),
		m_p(X.cols()),
		m_p_zero_one(X_zero_one.cols()),
		m_pi0(m_n),
		m_pi1(m_n),
		m_pib(m_n),
		m_grad_gamma0(m_p_zero_one),
		m_grad_gamma1(m_p_zero_one)
	{
		m_n_zero = 0;
		m_n_one = 0;
		std::vector<int> beta_idx;
		for (int i = 0; i < m_n; ++i){
			if (m_y[i] <= 0){
				++m_n_zero;
			} else if (m_y[i] >= 1){
				++m_n_one;
			} else {
				beta_idx.push_back(i);
			}
		}

		m_n_beta = beta_idx.size();
		m_X_beta.resize(m_n_beta, m_p);
		m_y_beta.resize(m_n_beta);
		for (int j = 0; j < m_n_beta; ++j){
			int row = beta_idx[j];
			m_X_beta.row(j) = m_X.row(row);
			m_y_beta[j] = m_y[row];
		}
		m_log_y_beta = m_y_beta.array().log().matrix();
		m_log1_y_beta = (1.0 - m_y_beta.array()).log().matrix();
	}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad){
		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		Eigen::VectorXd gamma0 = params.segment(m_p + 1, m_p_zero_one);
		Eigen::VectorXd gamma1 = params.tail(m_p_zero_one);

		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = logistic(eta);
		clamp_probs(mu);

		Eigen::VectorXd eta_beta = m_X_beta * beta;
		Eigen::VectorXd mu_beta = logistic(eta_beta);
		clamp_probs(mu_beta);

		Eigen::VectorXd eta0 = m_X_zero_one * gamma0;
		Eigen::VectorXd eta1 = m_X_zero_one * gamma1;
		double mixture_loglik = 0.0;
		m_grad_gamma0.setZero();
		m_grad_gamma1.setZero();
		double* g0 = m_grad_gamma0.data();
		double* g1 = m_grad_gamma1.data();
		for (int i = 0; i < m_n; ++i){
			double max_logit = std::max(std::max(eta0[i], eta1[i]), 0.0);
			double e0 = std::exp(eta0[i] - max_logit);
			double e1 = std::exp(eta1[i] - max_logit);
			double eb = std::exp(-max_logit);
			double denom = e0 + e1 + eb;
			m_pi0[i] = e0 / denom;
			m_pi1[i] = e1 / denom;
			m_pib[i] = eb / denom;
			const double* xzi_ptr = m_X_zero_one.data() + i;  // xzi_ptr[j * m_n] == X_zero_one(i,j)
			if (m_y[i] <= 0){
				mixture_loglik += std::log(m_pi0[i]);
				for (int j = 0; j < m_p_zero_one; ++j) {
					const double xj = xzi_ptr[j * m_n];
					g0[j] += xj * (1.0 - m_pi0[i]);
					g1[j] -= xj * m_pi1[i];
				}
			} else if (m_y[i] >= 1){
				mixture_loglik += std::log(m_pi1[i]);
				for (int j = 0; j < m_p_zero_one; ++j) {
					const double xj = xzi_ptr[j * m_n];
					g0[j] -= xj * m_pi0[i];
					g1[j] += xj * (1.0 - m_pi1[i]);
				}
			} else {
				mixture_loglik += std::log(m_pib[i]);
				for (int j = 0; j < m_p_zero_one; ++j) {
					const double xj = xzi_ptr[j * m_n];
					g0[j] -= xj * m_pi0[i];
					g1[j] -= xj * m_pi1[i];
				}
			}
		}

		Eigen::VectorXd mu_beta_phi = mu_beta.array() * phi;
		Eigen::VectorXd one_minus_mu_beta_phi = (1.0 - mu_beta.array()) * phi;
		double sum_loggamma_mu_phi = mu_beta_phi.unaryExpr(LgammaFunctor()).sum();
		double sum_loggamma_one_minus_mu_phi = one_minus_mu_beta_phi.unaryExpr(LgammaFunctor()).sum();
		double sum_term3 = ((mu_beta_phi.array() - 1.0) * m_log_y_beta.array()).sum();
		double sum_term4 = ((one_minus_mu_beta_phi.array() - 1.0) * m_log1_y_beta.array()).sum();
		double neg_ll_beta = -(
			m_n_beta * std::lgamma(phi) -
			sum_loggamma_mu_phi -
			sum_loggamma_one_minus_mu_phi +
			sum_term3 +
			sum_term4
		);

		double neg_ll = -mixture_loglik + neg_ll_beta;

		grad.resize(m_p + 1 + 2 * m_p_zero_one);

		Eigen::VectorXd d_mu_d_eta_beta = mu_beta.array() * (1.0 - mu_beta.array());
		Eigen::VectorXd digamma_mu_phi = mu_beta_phi.unaryExpr(DigammaFunctor());
		Eigen::VectorXd digamma_one_minus = one_minus_mu_beta_phi.unaryExpr(DigammaFunctor());

		Eigen::VectorXd d_neg_ll_d_mu_beta = (
			(digamma_mu_phi - digamma_one_minus - m_log_y_beta + m_log1_y_beta).array()
		) * phi;

		grad.head(m_p) = m_X_beta.transpose() *
			(d_neg_ll_d_mu_beta.array() * d_mu_d_eta_beta.array()).matrix();

		double sum_mu_digamma = (mu_beta.array() * digamma_mu_phi.array()).sum();
		double sum_one_minus_digamma = ((1.0 - mu_beta.array()) * digamma_one_minus.array()).sum();
		double sum_mu_logy = (mu_beta.array() * m_log_y_beta.array()).sum();
		double sum_one_minus_log1y = ((1.0 - mu_beta.array()) * m_log1_y_beta.array()).sum();

		double d_neg_ll_d_phi = -m_n_beta * fast_digamma(phi) +
			sum_mu_digamma +
			sum_one_minus_digamma -
			sum_mu_logy -
			sum_one_minus_log1y;
		grad[m_p] = d_neg_ll_d_phi * phi;
		grad.segment(m_p + 1, m_p_zero_one) = -m_grad_gamma0;
		grad.tail(m_p_zero_one) = -m_grad_gamma1;

		return neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params){
		Eigen::MatrixXd H(m_p + 1 + 2 * m_p_zero_one, m_p + 1 + 2 * m_p_zero_one);
		H.setZero();

		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		if (m_n_beta > 0){
			Eigen::VectorXd eta_beta = m_X_beta * beta;
			Eigen::VectorXd mu_beta = logistic(eta_beta);
			clamp_probs(mu_beta);

			double d_neg_ll_d_phi = -m_n_beta * fast_digamma(phi);
			double d2_neg_ll_d_phi2 = -m_n_beta * R::trigamma(phi);
			// total_H is the full H matrix column stride
			const int total_H = m_p + 1 + 2 * m_p_zero_one;
			double* H_data = H.data();

			for (int i = 0; i < m_n_beta; ++i){
				const double mu = mu_beta[i];
				const double one_minus_mu = 1.0 - mu;
				const double dmu_deta = mu * one_minus_mu;
				const double d2mu_deta2 = dmu_deta * (1.0 - 2.0 * mu);
				const double a = mu * phi;
				const double b = one_minus_mu * phi;
				const double digamma_a = fast_digamma(a);
				const double digamma_b = fast_digamma(b);
				const double trigamma_a = R::trigamma(a);
				const double trigamma_b = R::trigamma(b);
				const double c = digamma_a - digamma_b - m_log_y_beta[i] + m_log1_y_beta[i];
				const double dc_dmu = phi * (trigamma_a + trigamma_b);
				const double dc_dphi = mu * trigamma_a - one_minus_mu * trigamma_b;
				const double dscore_eta_deta =
					phi * (dc_dmu * dmu_deta * dmu_deta + c * d2mu_deta2);
				const double dscore_eta_dlogphi =
					phi * dmu_deta * (c + phi * dc_dphi);
				const double* xi = m_X_beta.data() + i;  // xi[j * m_n_beta] == X_beta(i,j)

				// Rank-1 update into top-left m_p x m_p block (upper triangle)
				for (int col = 0; col < m_p; ++col) {
					const double wxi_c = dscore_eta_deta * xi[col * m_n_beta];
					for (int row = 0; row <= col; ++row)
						H_data[row + col * total_H] += wxi_c * xi[row * m_n_beta];
				}
				// Vector update into top-right column m_p
				for (int row = 0; row < m_p; ++row)
					H_data[row + m_p * total_H] += dscore_eta_dlogphi * xi[row * m_n_beta];

				d_neg_ll_d_phi +=
					mu * digamma_a +
					one_minus_mu * digamma_b -
					mu * m_log_y_beta[i] -
					one_minus_mu * m_log1_y_beta[i];
				d2_neg_ll_d_phi2 +=
					mu * mu * trigamma_a +
					one_minus_mu * one_minus_mu * trigamma_b;
			}

			H(m_p, m_p) = phi * d_neg_ll_d_phi + phi * phi * d2_neg_ll_d_phi2;
			// Reflect upper triangle to lower for the beta block (top-left (m_p+1) x (m_p+1))
			for (int col = 0; col < m_p + 1; ++col)
				for (int row = 0; row < col; ++row)
					H_data[col + row * total_H] = H_data[row + col * total_H];
		}

		Eigen::VectorXd gamma0 = params.segment(m_p + 1, m_p_zero_one);
		Eigen::VectorXd gamma1 = params.tail(m_p_zero_one);
		Eigen::MatrixXd H00 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H11 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H01 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::VectorXd eta0 = m_X_zero_one * gamma0;
		Eigen::VectorXd eta1 = m_X_zero_one * gamma1;
		double* h00 = H00.data();
		double* h11 = H11.data();
		double* h01 = H01.data();
		for (int i = 0; i < m_n; ++i){
			double max_logit = std::max(std::max(eta0[i], eta1[i]), 0.0);
			double e0 = std::exp(eta0[i] - max_logit);
			double e1 = std::exp(eta1[i] - max_logit);
			double eb = std::exp(-max_logit);
			double denom = e0 + e1 + eb;
			double pi0 = e0 / denom;
			double pi1 = e1 / denom;
			const double* xzi_ptr = m_X_zero_one.data() + i;  // xzi_ptr[j * m_n] == X_zero_one(i,j)
			const double w00 = pi0 * (1.0 - pi0);
			const double w11 = pi1 * (1.0 - pi1);
			const double w01 = -pi0 * pi1;
			for (int c = 0; c < m_p_zero_one; ++c) {
				const double xj_c = xzi_ptr[c * m_n];
				for (int r = 0; r <= c; ++r) {
					const double xj_r = xzi_ptr[r * m_n];
					h00[r + c * m_p_zero_one] += w00 * xj_r * xj_c;
					h11[r + c * m_p_zero_one] += w11 * xj_r * xj_c;
					h01[r + c * m_p_zero_one] += w01 * xj_r * xj_c;
				}
			}
		}
		// Reflect upper triangle to lower for H00, H11, H01
		for (int c = 0; c < m_p_zero_one; ++c)
			for (int r = 0; r < c; ++r) {
				h00[c + r * m_p_zero_one] = h00[r + c * m_p_zero_one];
				h11[c + r * m_p_zero_one] = h11[r + c * m_p_zero_one];
				h01[c + r * m_p_zero_one] = h01[r + c * m_p_zero_one];
			}
		const int g0_start = m_p + 1;
		const int g1_start = m_p + 1 + m_p_zero_one;
		H.block(g0_start, g0_start, m_p_zero_one, m_p_zero_one) = H00;
		H.block(g1_start, g1_start, m_p_zero_one, m_p_zero_one) = H11;
		H.block(g0_start, g1_start, m_p_zero_one, m_p_zero_one) = H01;
		H.block(g1_start, g0_start, m_p_zero_one, m_p_zero_one) = H01.transpose();
		return H;
	}

	Eigen::MatrixXd expected_hessian(const Eigen::Ref<const Eigen::VectorXd>& params){
		Eigen::MatrixXd H(m_p + 1 + 2 * m_p_zero_one, m_p + 1 + 2 * m_p_zero_one);
		H.setZero();
		const int total_H = m_p + 1 + 2 * m_p_zero_one;
		double* H_data = H.data();

		Eigen::VectorXd beta = params.head(m_p);
		double phi = std::exp(params[m_p]);
		Eigen::VectorXd gamma0 = params.segment(m_p + 1, m_p_zero_one);
		Eigen::VectorXd gamma1 = params.tail(m_p_zero_one);
		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu_vec = logistic(eta);
		clamp_probs(mu_vec);
		Eigen::VectorXd eta0 = m_X_zero_one * gamma0;
		Eigen::VectorXd eta1 = m_X_zero_one * gamma1;

		Eigen::MatrixXd H00 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H11 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H01 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		double* h00 = H00.data();
		double* h11 = H11.data();
		double* h01 = H01.data();

		for (int i = 0; i < m_n; ++i){
			double max_logit = std::max(std::max(eta0[i], eta1[i]), 0.0);
			double e0 = std::exp(eta0[i] - max_logit);
			double e1 = std::exp(eta1[i] - max_logit);
			double eb = std::exp(-max_logit);
			double denom = e0 + e1 + eb;
			double pi0 = e0 / denom;
			double pi1 = e1 / denom;
			double pib = eb / denom;

			const double mu = mu_vec[i];
			const double one_minus_mu = 1.0 - mu;
			const double dmu_deta = mu * one_minus_mu;
			const double a = mu * phi;
			const double b = one_minus_mu * phi;
			const double trigamma_a = R::trigamma(a);
			const double trigamma_b = R::trigamma(b);
			const double h_beta = pib * phi * phi * (trigamma_a + trigamma_b) * dmu_deta * dmu_deta;
			const double h_cross = pib * phi * (a * trigamma_a - b * trigamma_b) * dmu_deta;
			const double h_logphi = pib * phi * phi * (
				-R::trigamma(phi) + mu * mu * trigamma_a + one_minus_mu * one_minus_mu * trigamma_b
			);
			const double* xi = m_X.data() + i;
			for (int col = 0; col < m_p; ++col) {
				const double wxi_c = h_beta * xi[col * m_n];
				for (int row = 0; row <= col; ++row)
					H_data[row + col * total_H] += wxi_c * xi[row * m_n];
				H_data[col + m_p * total_H] += h_cross * xi[col * m_n];
			}
			H(m_p, m_p) += h_logphi;

			const double* xzi_ptr = m_X_zero_one.data() + i;
			const double w00 = pi0 * (1.0 - pi0);
			const double w11 = pi1 * (1.0 - pi1);
			const double w01 = -pi0 * pi1;
			for (int c = 0; c < m_p_zero_one; ++c) {
				const double xj_c = xzi_ptr[c * m_n];
				for (int r = 0; r <= c; ++r) {
					const double xj_r = xzi_ptr[r * m_n];
					h00[r + c * m_p_zero_one] += w00 * xj_r * xj_c;
					h11[r + c * m_p_zero_one] += w11 * xj_r * xj_c;
					h01[r + c * m_p_zero_one] += w01 * xj_r * xj_c;
				}
			}
		}

		for (int col = 0; col < m_p + 1; ++col)
			for (int row = 0; row < col; ++row)
				H_data[col + row * total_H] = H_data[row + col * total_H];
		for (int c = 0; c < m_p_zero_one; ++c)
			for (int r = 0; r < c; ++r) {
				h00[c + r * m_p_zero_one] = h00[r + c * m_p_zero_one];
				h11[c + r * m_p_zero_one] = h11[r + c * m_p_zero_one];
				h01[c + r * m_p_zero_one] = h01[r + c * m_p_zero_one];
			}
		const int g0_start = m_p + 1;
		const int g1_start = m_p + 1 + m_p_zero_one;
		H.block(g0_start, g0_start, m_p_zero_one, m_p_zero_one) = H00;
		H.block(g1_start, g1_start, m_p_zero_one, m_p_zero_one) = H11;
		H.block(g0_start, g1_start, m_p_zero_one, m_p_zero_one) = H01;
		H.block(g1_start, g0_start, m_p_zero_one, m_p_zero_one) = H01.transpose();
		return H;
	}

private:
	Eigen::Ref<const Eigen::VectorXd> m_y;
	Eigen::Ref<const Eigen::MatrixXd> m_X;
	Eigen::Ref<const Eigen::MatrixXd> m_X_zero_one;
	Eigen::MatrixXd m_X_beta;
	Eigen::VectorXd m_y_beta;
	Eigen::VectorXd m_log_y_beta;
	Eigen::VectorXd m_log1_y_beta;
	int m_n;
	int m_p;
	int m_p_zero_one;
	int m_n_beta;
	int m_n_zero;
	int m_n_one;
	Eigen::VectorXd m_pi0;
	Eigen::VectorXd m_pi1;
	Eigen::VectorXd m_pib;
	Eigen::VectorXd m_grad_gamma0;
	Eigen::VectorXd m_grad_gamma1;
};

// [[Rcpp::export]]
SEXP get_zero_one_inflated_beta_score_cpp(SEXP X_sexp,
                                         SEXP X_zero_one_sexp,
                                         SEXP y_sexp,
                                         SEXP params_sexp){
	NumericMatrix X_mat(X_sexp);
	NumericMatrix X_zero_one_mat(X_zero_one_sexp);
	NumericVector y_vec(y_sexp);
	NumericVector params_vec(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::MatrixXd> X_zero_one(X_zero_one_mat.begin(), X_zero_one_mat.nrow(), X_zero_one_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());
	ZeroOneInflatedBeta fun(y, X, X_zero_one);
	Eigen::VectorXd grad(params.size());
	fun(params, grad);
	return wrap(-grad);
}

// [[Rcpp::export]]
SEXP get_zero_one_inflated_beta_hessian_cpp(SEXP X_sexp,
                                           SEXP X_zero_one_sexp,
                                           SEXP y_sexp,
                                           SEXP params_sexp){
	NumericMatrix X_mat(X_sexp);
	NumericMatrix X_zero_one_mat(X_zero_one_sexp);
	NumericVector y_vec(y_sexp);
	NumericVector params_vec(params_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
	Eigen::Map<const Eigen::MatrixXd> X_zero_one(X_zero_one_mat.begin(), X_zero_one_mat.nrow(), X_zero_one_mat.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
	Eigen::Map<const Eigen::VectorXd> params(params_vec.begin(), params_vec.size());
	ZeroOneInflatedBeta fun(y, X, X_zero_one);
	return wrap(-fun.hessian(params));
}

//' @title Fast Zero/One-Inflated Beta Regression (C++)
//' @description High-performance zero/one-inflated beta regression fitting using Newton-Raphson or L-BFGS.
//' @param X Matrix of predictors for the beta component.
//' @param X_zero_one Matrix of predictors for the zero and one inflation components.
//' @param y Vector of responses in [0, 1].
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first iteration.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_zero_one_inflated_beta_cpp(SEXP X_sexp,
									 SEXP X_zero_one_sexp,
									 SEXP y_sexp,
									 Nullable<NumericVector> warm_start_params = R_NilValue,
									 bool smart_cold_start = true,
									 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
									 std::string optimization_alg = "lbfgs",
									 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
									 bool estimate_only = false){

    NumericMatrix X_r(X_sexp);
    NumericMatrix X_zero_one_r(X_zero_one_sexp);
    NumericVector y_r(y_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::MatrixXd> X_zero_one(X_zero_one_r.begin(), X_zero_one_r.nrow(), X_zero_one_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y_eigen(y_r.begin(), y_r.size());

	int p = X.cols();
	int p_zero_one = X_zero_one.cols();
	int total = p + 1 + 2 * p_zero_one;
	Eigen::VectorXd params(total);
	
	if (warm_start_params.isNotNull()) {
		params = Rcpp::as<Eigen::VectorXd>(warm_start_params);
	} else if (smart_cold_start) {
		// Beta component: OLS on logit(y) for entries in (0, 1)
		std::vector<int> idx_mid;
		for(int i=0; i<y_eigen.size(); ++i) if(y_eigen[i] > 0 && y_eigen[i] < 1) idx_mid.push_back(i);
		if (idx_mid.size() > (size_t)p) {
			Eigen::MatrixXd X_mid(idx_mid.size(), p);
			Eigen::VectorXd y_logit(idx_mid.size());
			for(size_t i=0; i<idx_mid.size(); ++i) {
				X_mid.row(i) = X.row(idx_mid[i]);
				double yi = y_eigen[idx_mid[i]];
				y_logit[i] = std::log(yi / (1.0 - yi));
			}
			params.head(p) = safe_ols_solve(X_mid, y_logit);
		} else {
			params.head(p).setZero();
		}
		params[p] = 2.0; // log_phi warm_start_params
		
		// Zero/One components: OLS on indicators
		Eigen::VectorXd y_is_zero = (y_eigen.array() == 0.0).cast<double>();
		Eigen::VectorXd y_is_one  = (y_eigen.array() == 1.0).cast<double>();
		params.segment(p + 1, p_zero_one) = ols_smart_cold_start_beta(X_zero_one, y_is_zero);
		params.tail(p_zero_one)           = ols_smart_cold_start_beta(X_zero_one, y_is_one);
	} else {
		params.setZero();
		params[p] = 2.0;
	}
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	ZeroOneInflatedBeta fun(y_eigen, X, X_zero_one);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }
    
	LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, 1500, 1e-6, optimization_alg, "lbfgs", 0, info_start_ptr);
	params = fit.params;

	if (estimate_only) {
		return List::create(
			Named("b") = params.head(p),
			Named("log_phi") = params[p],
			Named("zero_one_b0") = params.segment(p + 1, p_zero_one),
			Named("zero_one_b1") = params.tail(p_zero_one),
			Named("params") = params,
			Named("neg_loglik") = fit.value,
			Named("converged") = fit.converged
		);
	}

	Eigen::MatrixXd observed_information = fun.hessian(params);

	int dim = params.size();
	NumericMatrix vcov_mat(dim, dim);
	bool has_vcov = false;
	if (observed_information.rows() == dim && observed_information.cols() == dim && observed_information.allFinite()){
		Eigen::MatrixXd H_free = subset_matrix(observed_information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
		if (lu.isInvertible()){
			Eigen::MatrixXd inv_free = lu.inverse();
			Eigen::MatrixXd inv = expand_free_covariance(dim, fixed_spec, inv_free, true);
			for (int i = 0; i < dim; ++i){
				for (int j = 0; j < dim; ++j){
					vcov_mat(i, j) = inv(i, j);
				}
			}
			has_vcov = true;
		}
	}
	if (!has_vcov){
		for (int i = 0; i < dim; ++i){
			for (int j = 0; j < dim; ++j){
				vcov_mat(i, j) = NA_REAL;
			}
		}
	}

	return List::create(
		Named("b") = params.head(p),
		Named("log_phi") = params[p],
		Named("zero_one_b0") = params.segment(p + 1, p_zero_one),
		Named("zero_one_b1") = params.tail(p_zero_one),
		Named("params") = params,
		Named("vcov") = vcov_mat,
		Named("neg_loglik") = fit.value,
		Named("observed_information") = observed_information,
		Named("fisher_information") = observed_information,
		Named("information") = observed_information,
		Named("information_type") = "observed",
		Named("hessian") = -observed_information
	);
}
