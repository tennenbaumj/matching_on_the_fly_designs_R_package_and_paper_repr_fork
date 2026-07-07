#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

double clamp_eta_for_exp(double eta) {
	return std::min(eta, 700.0);
}

class PoissonNegLogLik {
private:
    const Eigen::Ref<const MatrixXd> m_X;
    const Eigen::Ref<const VectorXd> m_y;
    const Eigen::Ref<const VectorXd> m_weights;
    const bool m_use_weights;
    const int m_n;

public:
    PoissonNegLogLik(const Eigen::Ref<const MatrixXd>& X,
                     const Eigen::Ref<const VectorXd>& y,
                     const Eigen::Ref<const VectorXd>& weights) :
        m_X(X), m_y(y), m_weights(weights), m_use_weights(weights.size() == X.rows()), m_n(X.rows()) {}

    double operator()(const VectorXd& beta, VectorXd& grad) {
        VectorXd eta = m_X * beta;
        double neg_ll = 0.0;
        VectorXd diff(m_n);
        
        if (m_use_weights) {
            for (int i = 0; i < m_n; ++i) {
                double ei = clamp_eta_for_exp(eta[i]);
                double mui = std::exp(ei);
                double wi = m_weights[i];
                neg_ll += wi * (mui - m_y[i] * ei);
                diff[i] = wi * (mui - m_y[i]);
            }
        } else {
            for (int i = 0; i < m_n; ++i) {
                double ei = clamp_eta_for_exp(eta[i]);
                double mui = std::exp(ei);
                neg_ll += (mui - m_y[i] * ei);
                diff[i] = (mui - m_y[i]);
            }
        }
        grad.noalias() = m_X.transpose() * diff;
        return neg_ll;
    }

    MatrixXd hessian(const VectorXd& beta) {
        VectorXd eta = m_X * beta;
        VectorXd w(m_n);
        if (m_use_weights) {
            for (int i = 0; i < m_n; ++i) {
                double ei = clamp_eta_for_exp(eta[i]);
                w[i] = std::exp(ei) * m_weights[i];
            }
        } else {
            for (int i = 0; i < m_n; ++i) {
                double ei = clamp_eta_for_exp(eta[i]);
                w[i] = std::exp(ei);
            }
        }
        return weighted_crossprod(m_X, w);
    }
};

ModelResult fast_poisson_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
							 const Eigen::Ref<const Eigen::VectorXd>& y,
                             const Eigen::Ref<const Eigen::VectorXd>& weights,
                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                             bool smart_cold_start = false,
							 int maxit = 100,
							 double tol = 1e-8,
                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                             std::string optimization_alg = "irls",
                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                             bool estimate_only = false) {
	const int n = X.rows();
	const int p = X.cols();
    bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::string alg = normalize_optimizer_algorithm(optimization_alg, "irls", true);
    VectorXd beta_start = VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        Rcpp::NumericVector ws_nv(warm_start_beta);
        if (ws_nv.size() > 0) {
            beta_start = as<Eigen::VectorXd>(ws_nv);
        }
    } else if (smart_cold_start) {
        beta_start = edi_opt::poisson_smart_cold_start(X, y);
    }
    beta_start = apply_fixed_values(beta_start, fixed_spec);

    if (alg != "irls") {
        VectorXd beta = beta_start;
        PoissonNegLogLik fun(X, y, weights);
        
        Eigen::MatrixXd H_start_val;
        const Eigen::MatrixXd* h_ptr = nullptr;
        if (warm_start_fisher_info.isNotNull()) {
            H_start_val = as<Eigen::MatrixXd>(warm_start_fisher_info);
            h_ptr = &H_start_val;
        } else if (smart_cold_start) {
            H_start_val = edi_opt::poisson_smart_hessian(X, beta_start);
            h_ptr = &H_start_val;
        }

        LikelihoodFitResult fit = (alg == "lbfgs") ?
            optimize_fixed_likelihood_lbfgs(fun, beta, fixed_spec, maxit, tol) :
            optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, alg, "newton_raphson", 0, h_ptr);

        ModelResult res;
        res.b = fit.params;
        res.neg_ll = fit.value;
        if (!estimate_only) {
            VectorXd eta = X * res.b;
            eta = eta.cwiseMin(700.0);
            res.mu = eta.array().exp().matrix();
            res.XtWX = fun.hessian(res.b);
            VectorXd diff = res.mu - y;
            if (use_weights) diff.array() *= weights.array();
            res.score = X.transpose() * diff;
        }
        res.iterations = fit.niter;
        res.converged = fit.converged;
        return res;
    }


    const int p_free = (int)fixed_spec.free_idx.size();
    const bool all_free = (p_free == p);
    
    MatrixXd X_free_storage;
    if (!all_free) {
        X_free_storage.resize(n, p_free);
        for (int j = 0; j < p_free; ++j) {
            X_free_storage.col(j) = X.col(fixed_spec.free_idx[j]);
        }
    }
    auto get_X_f = [&]() -> Eigen::Ref<const MatrixXd> {
        if (all_free) return X;
        return X_free_storage;
    };
    Eigen::Ref<const MatrixXd> X_f = get_X_f();

    VectorXd beta_free(p_free);
    for (int j = 0; j < p_free; ++j) beta_free[j] = beta_start[fixed_spec.free_idx[j]];
    
    VectorXd eta_fixed = VectorXd::Zero(n);
    if (fixed_spec.has_fixed) {
        for (int j = 0; j < (int)fixed_spec.fixed_idx.size(); ++j) {
            eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
        }
    }

    ModelResult res;
	res.b = VectorXd::Zero(p);
    VectorXd mu(n);
    VectorXd eta(n);
    VectorXd eta_try(n);
    VectorXd delta_eta(n);
    MatrixXd XtWX_free = MatrixXd::Zero(p_free, p_free);
    VectorXd score_free(p_free);

    auto compute_neg_loglik = [&](const VectorXd& b_free) {
        VectorXd e = eta_fixed + X_f * b_free;
        double nll = 0.0;
        for (int i = 0; i < n; ++i) {
            double ei = clamp_eta_for_exp(e[i]);
            double mui = std::exp(ei);
            double wi = use_weights ? weights[i] : 1.0;
            nll += wi * (mui - y[i] * ei);
        }
        return nll;
    };

    auto compute_neg_loglik_from_eta = [&](const VectorXd& e_full) {
        double nll = 0.0;
        for (int i = 0; i < n; ++i) {
            double ei = clamp_eta_for_exp(e_full[i]);
            double mui = std::exp(ei);
            double wi = use_weights ? weights[i] : 1.0;
            nll += wi * (mui - y[i] * ei);
        }
        return nll;
    };

    double current_nll = compute_neg_loglik(beta_free);

	for (int iter = 0; iter < maxit; ++iter) {
        res.iterations = iter + 1;
        
        eta.noalias() = eta_fixed;
        eta.noalias() += X_f * beta_free;
        mu = eta.array().cwiseMin(700.0).exp();
        
        VectorXd diff = y - mu;
        if (use_weights) diff.array() *= weights.array();
        score_free.noalias() = X_f.transpose() * diff;

        if (score_free.norm() < tol) {
            res.converged = true;
            break;
        }

        VectorXd w_tmp = use_weights ? mu.cwiseProduct(weights) : mu;
        w_tmp = w_tmp.array().cwiseMax(1e-10);

        const bool use_warm_start_xtwx = (iter == 0) &&
            (warm_start_fisher_info.isNotNull() ||
             (smart_cold_start && warm_start_beta.isNull()));

        if (use_warm_start_xtwx) {
            if (warm_start_fisher_info.isNotNull()) {
                Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
                if (info_full.rows() == p && info_full.cols() == p) {
                    XtWX_free = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
                } else {
                    XtWX_free = weighted_crossprod(X_f, w_tmp);
                }
            } else {
                Eigen::MatrixXd H_full = edi_opt::poisson_smart_hessian(X, beta_start);
                XtWX_free = subset_matrix(H_full, fixed_spec.free_idx, fixed_spec.free_idx);
            }
        } else {
            XtWX_free = weighted_crossprod(X_f, w_tmp);
        }
        
		Eigen::LDLT<MatrixXd> ldlt(XtWX_free);
        if (ldlt.info() != Eigen::Success) break;
        VectorXd step = ldlt.solve(score_free);
		if (!step.allFinite()) break;

        // Step-halving: delta_eta cached so each probe is O(n) not O(n·p)
        delta_eta.noalias() = X_f * step;
        double step_size = 1.0;
        bool step_ok = false;
        for (int s = 0; s < 10; ++s) {
            eta_try.noalias() = eta + step_size * delta_eta;
            double try_nll = compute_neg_loglik_from_eta(eta_try);
            if (try_nll < current_nll + 1e-10) {
                beta_free.noalias() += step_size * step;
                current_nll = try_nll;
                step_ok = true;
                break;
            }
            step_size *= 0.5;
        }
        
        if (!step_ok) break;
		if (step.norm() < tol) {
			res.converged = true;
			break;
		}
	}

    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
    for (int j = 0; j < (int)fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    if (!estimate_only) {
        eta.noalias() = X * res.b;
        res.mu = eta.array().cwiseMin(700.0).exp().matrix();
        res.mu = res.mu.array().cwiseMax(1e-10);
        
        VectorXd w_final = use_weights ? res.mu.cwiseProduct(weights) : res.mu;
        w_final = w_final.array().cwiseMax(1e-10);
        
        res.neg_ll = current_nll;
        MatrixXd info_free(p_free, p_free);
        info_free = weighted_crossprod(X_f, w_final);
        res.XtWX = expand_free_covariance(p, fixed_spec, info_free, false);
        
        VectorXd final_diff = y - res.mu;
        if (use_weights) final_diff.array() *= weights.array();
        res.score = X.transpose() * final_diff;
    }

	return res;
}

} // namespace

ModelResult fast_poisson_regression_internal(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                             const Eigen::Ref<const Eigen::VectorXd>& weights,
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                             bool smart_cold_start = false,
                                             int maxit = 100,
                                             double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "irls",
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                             bool estimate_only = false) {
    return fast_poisson_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
}

// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_score_cpp(SEXP X_sexp, SEXP y_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp); NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * (y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_hessian_cpp(SEXP X_sexp, SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp); NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return -weighted_crossprod(X, mu);
}

// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_weighted_score_cpp(SEXP X_sexp,
                                                          SEXP y_sexp,
                                                          SEXP weights_sexp,
                                                          SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp); NumericVector weights_r(weights_sexp); NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
    Eigen::Map<const Eigen::VectorXd> weights(weights_r.begin(), weights_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * weights.cwiseProduct(y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_weighted_hessian_cpp(SEXP X_sexp,
                                                            SEXP weights_sexp,
                                                            SEXP beta_sexp) {
    NumericMatrix X_r(X_sexp); NumericVector weights_r(weights_sexp); NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
    Eigen::Map<const Eigen::VectorXd> weights(weights_r.begin(), weights_r.size());
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    Eigen::VectorXd w = weights.cwiseProduct(mu);
    return -weighted_crossprod(X, w);
}

// [[Rcpp::export]]
List fast_poisson_regression_cpp(SEXP X_sexp, SEXP y_sexp,
                                 Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                 bool smart_cold_start = false,
                                 int maxit = 100,
                                 double tol = 1e-8,
                                 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "irls",
                                     Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                     bool estimate_only = false) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
	if (estimate_only) {
		return List::create(
			Named("b") = res.b,
			Named("converged") = res.converged,
			Named("iterations") = res.iterations
		);
	}
	return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}

// [[Rcpp::export]]
List fast_poisson_regression_weighted_cpp(SEXP X_sexp,
                                          SEXP y_sexp,
                                          SEXP weights_sexp,
                                          Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                          bool smart_cold_start = false,
                                          int maxit = 100,
                                          double tol = 1e-8,
                                          Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                          std::string optimization_alg = "irls",
                                          Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp); NumericVector weights_r(weights_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> weights(weights_r.begin(), weights_r.size());
	ModelResult res = fast_poisson_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
	return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("score") = res.score,
        Named("neg_ll") = res.neg_ll,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}

// [[Rcpp::export]]
List fast_poisson_regression_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int j = 2,
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                              bool smart_cold_start = false,
											  int maxit = 100,
											  double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls",
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj + 1;
        return -1;
    };
    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;
    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;

	return List::create(
		Named("b") = res.b,
		Named("params") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations,
		Named("score") = res.score,
		Named("observed_information") = res.XtWX,
		Named("fisher_information") = res.XtWX,
		Named("information") = res.XtWX,
		Named("information_type") = "fisher",
		Named("hessian") = -res.XtWX,
		Named("neg_loglik") = res.neg_ll,
		Named("neg_ll") = res.neg_ll,
		Named("loglik") = R_finite(res.neg_ll) ? -res.neg_ll : NA_REAL
	);
}

// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(SEXP X_sexp,
                                               SEXP y_sexp,
                                               int j = 2,
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                               bool smart_cold_start = false,
                                               int maxit = 100,
                                               double tol = 1e-8,
                                               Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                               std::string optimization_alg = "irls",
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj + 1;
        return -1;
    };

	const int df_resid = X.rows() - X.cols();
	if (df_resid > 0) {
        double pearson_sum = 0.0;
        int n = X.rows();
        for (int i = 0; i < n; ++i) {
            double diff = y[i] - res.mu[i];
            pearson_sum += (diff * diff) / res.mu[i];
        }
		res.dispersion = pearson_sum / static_cast<double>(df_resid);

		if (std::isfinite(res.dispersion) && res.dispersion > 0) {
            int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
            int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;

            if (free_j > 0) {
                res.ssq_b_j = res.dispersion * compute_diagonal_inverse_entry(info_free, free_j);
            } else {
                res.ssq_b_j = NA_REAL;
            }

            if (free_2 > 0) {
                if (free_2 == free_j) {
                    res.ssq_b_2 = res.ssq_b_j;
                } else {
                    res.ssq_b_2 = res.dispersion * compute_diagonal_inverse_entry(info_free, free_2);
                }
            }
		}
	}

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("dispersion") = res.dispersion,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}
