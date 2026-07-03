// [[Rcpp::depends(RcppEigen)]]
#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

constexpr double kEps = 1e-8;
constexpr double kMinMu = kEps;
constexpr double kMaxMu = 1.0 - kEps;
constexpr double kMaxEtaLog = -kEps;

enum class BinomialConstrainedLink {
  kLog = 0,
  kIdentity = 1
};

inline double clamp_prob(double mu) {
  if (mu <= kMinMu) return kMinMu;
  if (mu >= kMaxMu) return kMaxMu;
  return mu;
}

inline double safe_mu_from_eta(double eta, BinomialConstrainedLink link_type) {
  if (link_type == BinomialConstrainedLink::kLog) {
    if (eta >= kMaxEtaLog) return std::exp(kMaxEtaLog);
    return std::exp(eta);
  }
  return clamp_prob(eta);
}

double loglik_constrained_binomial(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                   const Eigen::Ref<const Eigen::VectorXd>& y,
                                   const Eigen::Ref<const Eigen::VectorXd>& beta,
                                   BinomialConstrainedLink link_type) {
  const int n = (int)X.rows();
  Eigen::VectorXd eta = X * beta;
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    double ei = eta[i];
    if (link_type == BinomialConstrainedLink::kLog) {
      if (ei >= kMaxEtaLog) return R_NegInf;
      const double mu = std::exp(ei);
      ll += y[i] * ei + (1.0 - y[i]) * std::log1p(-mu);
    } else {
      if (ei <= kMinMu || ei >= kMaxMu) return R_NegInf;
      ll += y[i] * std::log(ei) + (1.0 - y[i]) * std::log1p(-ei);
    }
  }
  return ll;
}

double weighted_loglik_constrained_binomial(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                            const Eigen::Ref<const Eigen::VectorXd>& y,
                                            const Eigen::Ref<const Eigen::VectorXd>& obs_weights,
                                            const Eigen::Ref<const Eigen::VectorXd>& beta,
                                            BinomialConstrainedLink link_type) {
  const int n = (int)X.rows();
  Eigen::VectorXd eta = X * beta;
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    const double wi = obs_weights[i];
    if (!R_finite(wi) || wi < 0.0) return R_NegInf;
    double ei = eta[i];
    if (link_type == BinomialConstrainedLink::kLog) {
      if (ei >= kMaxEtaLog) return R_NegInf;
      const double mu = std::exp(ei);
      ll += wi * (y[i] * ei + (1.0 - y[i]) * std::log1p(-mu));
    } else {
      if (ei <= kMinMu || ei >= kMaxMu) return R_NegInf;
      ll += wi * (y[i] * std::log(ei) + (1.0 - y[i]) * std::log1p(-ei));
    }
  }
  return ll;
}

bool all_finite_vec(const Eigen::Ref<const Eigen::VectorXd>& x) {
  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i])) return false;
  }
  return true;
}

bool all_finite_mat(const Eigen::Ref<const Eigen::MatrixXd>& X) {
  for (int j = 0; j < X.cols(); ++j) {
    for (int i = 0; i < X.rows(); ++i) {
      if (!R_finite(X(i, j))) return false;
    }
  }
  return true;
}

// Evaluate log-likelihood directly from precomputed eta (avoids X*beta GEMV).
inline double loglik_from_eta(const Eigen::VectorXd& eta,
                               const Eigen::Ref<const Eigen::VectorXd>& y,
                               BinomialConstrainedLink link_type) {
  const int n = (int)eta.size();
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ei = eta[i];
    if (link_type == BinomialConstrainedLink::kLog) {
      if (ei >= kMaxEtaLog) return R_NegInf;
      ll += y[i] * ei + (1.0 - y[i]) * std::log1p(-std::exp(ei));
    } else {
      if (ei <= kMinMu || ei >= kMaxMu) return R_NegInf;
      ll += y[i] * std::log(ei) + (1.0 - y[i]) * std::log1p(-ei);
    }
  }
  return ll;
}

// Weighted version of loglik_from_eta for bootstrap/IPW paths.
inline double weighted_loglik_from_eta(const Eigen::VectorXd& eta,
                                       const Eigen::Ref<const Eigen::VectorXd>& y,
                                       const Eigen::Ref<const Eigen::VectorXd>& obs_weights,
                                       BinomialConstrainedLink link_type) {
  const int n = (int)eta.size();
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    const double wi = obs_weights[i];
    if (!R_finite(wi) || wi < 0.0) return R_NegInf;

    const double ei = eta[i];
    if (link_type == BinomialConstrainedLink::kLog) {
      if (ei >= kMaxEtaLog) return R_NegInf;
      ll += wi * (y[i] * ei + (1.0 - y[i]) * std::log1p(-std::exp(ei)));
    } else {
      if (ei <= kMinMu || ei >= kMaxMu) return R_NegInf;
      ll += wi * (y[i] * std::log(ei) + (1.0 - y[i]) * std::log1p(-ei));
    }
  }
  return ll;
}

List fit_constrained_binomial_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                       const Eigen::Ref<const Eigen::VectorXd>& y,
                                       BinomialConstrainedLink link_type,
                                       int maxit,
                                       double tol,
                                       Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                       bool smart_cold_start = true,
                                       Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                       Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                       bool estimate_only = false) {
  const int n = (int)X.rows();
  const int p = (int)X.cols();
  if (y.size() != n) stop("dimension mismatch in constrained binomial regression");
  FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
  const int p_free = (int)fixed_spec.free_idx.size();
  Eigen::MatrixXd X_free(n, p_free);
  for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
  Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
  for (int j = 0; j < (int)fixed_spec.fixed_idx.size(); ++j) {
    eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  if (warm_start_beta.isNotNull()) {
    beta = as<Eigen::VectorXd>(warm_start_beta);
    if (beta.size() != p) stop("warm_start_beta must have length equal to ncol(X)");
  } else if (smart_cold_start) {
    if (link_type == BinomialConstrainedLink::kLog) {
      beta = ols_smart_cold_start_beta_on_log1p_or_legacy(X, y, Eigen::VectorXd::Zero(p), fixed_spec);
    } else {
      beta = ols_smart_cold_start_beta_or_legacy(X, y, Eigen::VectorXd::Zero(p), fixed_spec);
    }
  } else {
    const double y_mean = std::min(std::max(y.mean(), kMinMu), kMaxMu);
    beta[0] = (link_type == BinomialConstrainedLink::kLog) ? std::log(y_mean) : y_mean;
  }
  beta = apply_fixed_values(beta, fixed_spec);
  Eigen::VectorXd beta_free = subset_vector(beta, fixed_spec.free_idx);

  bool converged = false;
  Eigen::VectorXd mu = Eigen::VectorXd::Constant(n, y.mean());
  Eigen::VectorXd w = Eigen::VectorXd::Constant(n, 1.0);
  Eigen::VectorXd z_adj(n);

  // Cache ll at current beta — carried forward across iterations to avoid one
  // loglik_constrained_binomial call per accepted step.
  double ll_curr = loglik_constrained_binomial(X, y, beta, link_type);

  int iterations = 0;
  for (int iter = 0; iter < maxit; ++iter) {
    iterations = iter + 1;
    const Eigen::VectorXd eta = eta_fixed + X_free * beta_free;

    // In-place update of mu and w to avoid temporary arrays
    for (int i = 0; i < n; ++i) {
        double ei = eta[i];
        if (link_type == BinomialConstrainedLink::kLog) {
            if (ei > kMaxEtaLog) ei = kMaxEtaLog;
            double mui = std::exp(ei);
            if (mui < kEps) mui = kEps;
            mu[i] = mui;
            w[i] = std::max(mui / std::max(1.0 - mui, kEps), kEps);
            z_adj[i] = ei + (y[i] - mui) / mui - eta_fixed[i];
        } else {
            if (ei < kMinMu) ei = kMinMu;
            if (ei > kMaxMu) ei = kMaxMu;
            mu[i] = ei;
            w[i] = 1.0 / std::max(ei * (1.0 - ei), kEps);
            z_adj[i] = y[i] - eta_fixed[i];
        }
    }

    Eigen::MatrixXd XtWX;
    bool used_warm_fisher = false;
    if (iter == 0 && warm_start_fisher_info.isNotNull()) {
      Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
      if (info_full.rows() == p && info_full.cols() == p) {
        XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        used_warm_fisher = true;
      } else {
        XtWX = weighted_crossprod(X_free, w);
      }
    } else {
      XtWX = weighted_crossprod(X_free, w);
    }
    Eigen::VectorXd XtWz = weighted_crossprod_rhs(X_free, w, z_adj);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) {
      if (used_warm_fisher) {
        XtWX = weighted_crossprod(X_free, w);
        ldlt.compute(XtWX);
        if (ldlt.info() != Eigen::Success) {
          return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
        }
      } else {
        return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
      }
    }

    Eigen::VectorXd beta_free_target = ldlt.solve(XtWz);
    if (ldlt.info() != Eigen::Success || !all_finite_vec(beta_free_target)) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    // Precompute the step direction in eta-space once per IRLS iteration so
    // each backtracking probe is a cheap O(n) vector-add + scalar scan rather
    // than a full O(n*p) GEMV inside loglik_constrained_binomial.
    const Eigen::VectorXd delta_eta = X_free * (beta_free_target - beta_free);

    double step = 1.0;

    Eigen::VectorXd beta_new = beta;
    Eigen::VectorXd beta_free_new = beta_free;
    Eigen::VectorXd eta_try(n);
    bool accepted = false;
    while (step >= 1e-8) {
      eta_try.noalias() = eta + step * delta_eta;
      const double ll_new = loglik_from_eta(eta_try, y, link_type);
      if (R_finite(ll_new) && ll_new >= ll_curr - 1e-10) {
        beta_free_new = beta_free + step * (beta_free_target - beta_free);
        beta_new = expand_free_params(beta_free_new, beta, fixed_spec);
        ll_curr = ll_new;  // carry forward — avoids recomputing at next iter start
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;
    if ((beta_free_new - beta_free).norm() < tol) {
      beta = beta_new;
      beta_free = beta_free_new;
      converged = true;
      break;
    }
    beta = beta_new;
    beta_free = beta_free_new;
  }

  if (estimate_only) {
    return List::create(
      _["b"] = beta,
      _["converged"] = converged && all_finite_vec(beta),
      _["iterations"] = iterations
    );
  }

  Eigen::VectorXd eta = X * beta;
  if (link_type == BinomialConstrainedLink::kLog) {
    eta = eta.array().min(kMaxEtaLog).matrix();
    mu = eta.array().exp().matrix();
    w = (mu.array() / (1.0 - mu.array()).max(kEps)).max(kEps).matrix();
  } else {
    eta = eta.array().max(kMinMu).min(kMaxMu).matrix();
    mu = eta;
    w = (1.0 / (mu.array() * (1.0 - mu.array())).max(kEps)).max(kEps).matrix();
  }

  return List::create(
    _["b"] = beta,
    _["mu_hat"] = mu,
    _["working_weights"] = w,
    _["iterations"] = iterations,
    _["converged"] = converged && all_finite_vec(beta) && all_finite_vec(mu) && all_finite_vec(w),
    _["fisher_information"] = weighted_crossprod(X, w)
  );
}

List fit_constrained_binomial_weighted_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                const Eigen::Ref<const Eigen::VectorXd>& y,
                                                const Eigen::VectorXd& obs_weights,
                                                BinomialConstrainedLink link_type,
                                                int maxit,
                                                double tol,
                                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                bool smart_cold_start = true,
                                                Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                                bool estimate_only = false) {
  const int n = (int)X.rows();
  const int p = (int)X.cols();
  if (y.size() != n) stop("dimension mismatch in constrained binomial regression");
  if (obs_weights.size() != n) stop("weights length mismatch in constrained binomial regression");
  FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
  const int p_free = (int)fixed_spec.free_idx.size();
  Eigen::MatrixXd X_free(n, p_free);
  for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
  Eigen::VectorXd eta_fixed = Eigen::VectorXd::Zero(n);
  for (int j = 0; j < (int)fixed_spec.fixed_idx.size(); ++j) {
    eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  if (warm_start_beta.isNotNull()) {
    beta = as<Eigen::VectorXd>(warm_start_beta);
    if (beta.size() != p) stop("warm_start_beta must have length equal to ncol(X)");
  } else if (smart_cold_start) {
    if (link_type == BinomialConstrainedLink::kLog) {
      beta = ols_smart_cold_start_beta_on_log1p_or_legacy(X, y, Eigen::VectorXd::Zero(p), fixed_spec);
    } else {
      beta = ols_smart_cold_start_beta_or_legacy(X, y, Eigen::VectorXd::Zero(p), fixed_spec);
    }
  } else {
    const double y_mean = std::min(std::max(y.mean(), kMinMu), kMaxMu);
    beta[0] = (link_type == BinomialConstrainedLink::kLog) ? std::log(y_mean) : y_mean;
  }
  beta = apply_fixed_values(beta, fixed_spec);
  Eigen::VectorXd beta_free = subset_vector(beta, fixed_spec.free_idx);

  bool converged = false;
  Eigen::VectorXd mu = Eigen::VectorXd::Constant(n, y.mean());
  Eigen::VectorXd w = Eigen::VectorXd::Constant(n, 1.0);
  Eigen::VectorXd z_adj(n);
  Eigen::VectorXd w_eff(n);
  Eigen::VectorXd eta_try(n);
  Eigen::VectorXd warm_weights_vec;
  const bool has_warm_start_weights = warm_start_weights.isNotNull();
  if (has_warm_start_weights) {
    warm_weights_vec = as<Eigen::VectorXd>(warm_start_weights);
    if (warm_weights_vec.size() != n) stop("warm_start_weights must have length equal to nrow(X)");
  }

  double ll_curr = weighted_loglik_constrained_binomial(X, y, obs_weights, beta, link_type);

  int iterations = 0;
  for (int iter = 0; iter < maxit; ++iter) {
    iterations = iter + 1;
    const Eigen::VectorXd eta = eta_fixed + X_free * beta_free;
    const bool use_warm_weights = (iter == 0 && has_warm_start_weights);
    for (int i = 0; i < n; ++i) {
      double ei = eta[i];
      if (link_type == BinomialConstrainedLink::kLog) {
        if (ei > kMaxEtaLog) ei = kMaxEtaLog;
        double mui = std::exp(ei);
        if (mui < kEps) mui = kEps;
        mu[i] = mui;
        w[i] = use_warm_weights ? warm_weights_vec[i] : std::max(mui / std::max(1.0 - mui, kEps), kEps);
        z_adj[i] = ei + (y[i] - mui) / mui - eta_fixed[i];
      } else {
        if (ei < kMinMu) ei = kMinMu;
        if (ei > kMaxMu) ei = kMaxMu;
        mu[i] = ei;
        w[i] = use_warm_weights ? warm_weights_vec[i] : 1.0 / std::max(ei * (1.0 - ei), kEps);
        z_adj[i] = y[i] - eta_fixed[i];
      }
      w_eff[i] = obs_weights[i] * w[i];
    }

    Eigen::MatrixXd XtWX;
    bool used_warm_fisher_w = false;
    if (iter == 0 && warm_start_fisher_info.isNotNull()) {
      Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
      if (info_full.rows() == p && info_full.cols() == p) {
        XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        used_warm_fisher_w = true;
      } else {
        XtWX = weighted_crossprod(X_free, w_eff);
      }
    } else {
      XtWX = weighted_crossprod(X_free, w_eff);
    }
    Eigen::VectorXd XtWz = weighted_crossprod_rhs(X_free, w_eff, z_adj);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) {
      if (used_warm_fisher_w) {
        XtWX = weighted_crossprod(X_free, w_eff);
        ldlt.compute(XtWX);
        if (ldlt.info() != Eigen::Success) {
          return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
        }
      } else {
        return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
      }
    }

    Eigen::VectorXd beta_free_target = ldlt.solve(XtWz);
    if (ldlt.info() != Eigen::Success || !all_finite_vec(beta_free_target)) {
      return List::create(_["b"] = beta, _["mu_hat"] = mu, _["working_weights"] = w, _["converged"] = false);
    }

    const Eigen::VectorXd delta_beta = beta_free_target - beta_free;
    const Eigen::VectorXd delta_eta = X_free * delta_beta;
    double step = 1.0;
    Eigen::VectorXd beta_new = beta;
    Eigen::VectorXd beta_free_new = beta_free;
    bool accepted = false;
    while (step >= 1e-8) {
      eta_try.noalias() = eta + step * delta_eta;
      const double ll_new = weighted_loglik_from_eta(eta_try, y, obs_weights, link_type);
      if (R_finite(ll_new) && ll_new >= ll_curr - 1e-10) {
        beta_free_new = beta_free + step * delta_beta;
        beta_new = expand_free_params(beta_free_new, beta, fixed_spec);
        ll_curr = ll_new;
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) break;
    if ((beta_free_new - beta_free).norm() < tol) {
      beta = beta_new;
      beta_free = beta_free_new;
      converged = true;
      break;
    }
    beta = beta_new;
    beta_free = beta_free_new;
  }

  if (estimate_only) {
    return List::create(
      _["b"] = beta,
      _["converged"] = converged && all_finite_vec(beta),
      _["iterations"] = iterations
    );
  }

  Eigen::VectorXd eta = X * beta;
  if (link_type == BinomialConstrainedLink::kLog) {
    eta = eta.array().min(kMaxEtaLog).matrix();
    mu = eta.array().exp().matrix();
    w = (mu.array() / (1.0 - mu.array()).max(kEps)).max(kEps).matrix();
  } else {
    eta = eta.array().max(kMinMu).min(kMaxMu).matrix();
    mu = eta;
    w = (1.0 / (mu.array() * (1.0 - mu.array())).max(kEps)).max(kEps).matrix();
  }
  w_eff = obs_weights.cwiseProduct(w);

  return List::create(
    _["b"] = beta,
    _["mu_hat"] = mu,
    _["working_weights"] = w,
    _["iterations"] = iterations,
    _["converged"] = converged && all_finite_vec(beta) && all_finite_vec(mu) && all_finite_vec(w),
    _["fisher_information"] = weighted_crossprod(X, w_eff)
  );
}

List fit_constrained_binomial_with_var_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                 const Eigen::Ref<const Eigen::VectorXd>& y,
                                                 BinomialConstrainedLink link_type,
                                                 int j,
                                                 int maxit,
                                                 double tol,
                                                 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                 Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                 bool smart_cold_start = true,
                                                 Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                                 Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
  List fit = fit_constrained_binomial_cpp_impl(X, y, link_type, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info);
  const bool converged = as<bool>(fit["converged"]);
  Eigen::VectorXd beta = fit["b"];
  Eigen::VectorXd w = fit["working_weights"];

  if (!converged || !all_finite_vec(beta) || !all_finite_vec(w)) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  FixedParamSpec fixed_spec = make_fixed_param_spec((int)X.cols(), fixed_idx, fixed_values);
  Eigen::MatrixXd X_free(X.rows(), (int)fixed_spec.free_idx.size());
  for (int col = 0; col < (int)fixed_spec.free_idx.size(); ++col) X_free.col(col) = X.col(fixed_spec.free_idx[col]);
  Eigen::MatrixXd XtWX_free = weighted_crossprod(X_free, w);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX_free);
  if (ldlt.info() != Eigen::Success) {
    return List::create(
      _["b"] = beta,
      _["vcov"] = NumericMatrix(0, 0),
      _["std_err"] = NumericVector(0),
      _["z_vals"] = NumericVector(0),
      _["ssq_b_j"] = NA_REAL,
      _["converged"] = false
    );
  }

  int free_j = -1;
  for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
    if (fixed_spec.free_idx[jj] == j - 1) { free_j = jj + 1; break; }
  const double ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(XtWX_free, free_j) : NA_REAL;

  return List::create(
    _["b"] = beta,
    _["ssq_b_j"] = ssq_b_j,
    _["converged"] = true,
    _["fisher_information"] = weighted_crossprod(X, w),
    _["neg_ll"] = -loglik_constrained_binomial(X, y, beta, link_type),
    _["logLik"] = loglik_constrained_binomial(X, y, beta, link_type)
  );
}

Eigen::VectorXd constrained_binomial_score_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
													const Eigen::Ref<const Eigen::VectorXd>& y,
													const Eigen::Ref<const Eigen::VectorXd>& beta,
													BinomialConstrainedLink link_type) {
	const int p = (int)beta.size();
	Eigen::VectorXd score(p);
	const double h = 1e-6;
	for (int j = 0; j < p; ++j) {
		Eigen::VectorXd bp = beta;
		Eigen::VectorXd bm = beta;
		bp[j] += h;
		bm[j] -= h;
		score[j] = (loglik_constrained_binomial(X, y, bp, link_type) -
					loglik_constrained_binomial(X, y, bm, link_type)) / (2.0 * h);
	}
	return score;
}

Eigen::MatrixXd constrained_binomial_hessian_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
													  const Eigen::Ref<const Eigen::VectorXd>& y,
													  const Eigen::Ref<const Eigen::VectorXd>& beta,
													  BinomialConstrainedLink link_type) {
	const int p = (int)beta.size();
	Eigen::MatrixXd H(p, p);
	const double h = 1e-4;
	for (int i = 0; i < p; ++i) {
		for (int j = i; j < p; ++j) {
			Eigen::VectorXd bpp = beta; bpp[i] += h; bpp[j] += h;
			Eigen::VectorXd bpm = beta; bpm[i] += h; bpm[j] -= h;
			Eigen::VectorXd bmp = beta; bmp[i] -= h; bmp[j] += h;
			Eigen::VectorXd bmm = beta; bmm[i] -= h; bmm[j] -= h;
			H(i, j) = (loglik_constrained_binomial(X, y, bpp, link_type) -
					   loglik_constrained_binomial(X, y, bpm, link_type) -
					   loglik_constrained_binomial(X, y, bmp, link_type) +
					   loglik_constrained_binomial(X, y, bmm, link_type)) / (4.0 * h * h);
			H(j, i) = H(i, j);
		}
	}
	return H;
}

Eigen::VectorXd constrained_binomial_weighted_score_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                                             const Eigen::Ref<const Eigen::VectorXd>& weights,
                                                             const Eigen::Ref<const Eigen::VectorXd>& beta,
                                                             BinomialConstrainedLink link_type) {
  const int p = (int)beta.size();
  Eigen::VectorXd score(p);
  const double h = 1e-6;
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd bp = beta;
    Eigen::VectorXd bm = beta;
    bp[j] += h;
    bm[j] -= h;
    score[j] = (weighted_loglik_constrained_binomial(X, y, weights, bp, link_type) -
                weighted_loglik_constrained_binomial(X, y, weights, bm, link_type)) / (2.0 * h);
  }
  return score;
}

Eigen::MatrixXd constrained_binomial_weighted_hessian_cpp_impl(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                               const Eigen::Ref<const Eigen::VectorXd>& y,
                                                               const Eigen::Ref<const Eigen::VectorXd>& weights,
                                                               const Eigen::Ref<const Eigen::VectorXd>& beta,
                                                               BinomialConstrainedLink link_type) {
  const int p = (int)beta.size();
  Eigen::MatrixXd H(p, p);
  const double h = 1e-4;
  for (int i = 0; i < p; ++i) {
    for (int j = i; j < p; ++j) {
      Eigen::VectorXd bpp = beta; bpp[i] += h; bpp[j] += h;
      Eigen::VectorXd bpm = beta; bpm[i] += h; bpm[j] -= h;
      Eigen::VectorXd bmp = beta; bmp[i] -= h; bmp[j] += h;
      Eigen::VectorXd bmm = beta; bmm[i] -= h; bmm[j] -= h;
      H(i, j) = (weighted_loglik_constrained_binomial(X, y, weights, bpp, link_type) -
                 weighted_loglik_constrained_binomial(X, y, weights, bpm, link_type) -
                 weighted_loglik_constrained_binomial(X, y, weights, bmp, link_type) +
                 weighted_loglik_constrained_binomial(X, y, weights, bmm, link_type)) / (4.0 * h * h);
      H(j, i) = H(i, j);
    }
  }
  return H;
}

}  // namespace

//' @title Compute Log-Binomial Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_log_binomial_regression_score_cpp(SEXP X_r,
													  SEXP y_r,
													  SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
	return constrained_binomial_score_cpp_impl(X, y, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Log-Binomial Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_log_binomial_regression_hessian_cpp(SEXP X_r,
														SEXP y_r,
														SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
	return constrained_binomial_hessian_cpp_impl(X, y, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Weighted Log-Binomial Regression Score
//' @description Calculates the weighted score vector for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the weighted score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_log_binomial_regression_weighted_score_cpp(SEXP X_r,
                                                               SEXP y_r,
                                                               SEXP weights_r,
                                                               SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
  return constrained_binomial_weighted_score_cpp_impl(X, y, weights, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Weighted Log-Binomial Regression Hessian
//' @description Calculates the weighted Hessian matrix for a log-binomial regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the weighted Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_log_binomial_regression_weighted_hessian_cpp(SEXP X_r,
                                                                 SEXP y_r,
                                                                 SEXP weights_r,
                                                                 SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
  return constrained_binomial_weighted_hessian_cpp_impl(X, y, weights, beta, BinomialConstrainedLink::kLog);
}

//' @title Compute Identity-Binomial Regression Score
//' @description Calculates the score vector for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_identity_binomial_regression_score_cpp(SEXP X_r,
														   SEXP y_r,
														   SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
	return constrained_binomial_score_cpp_impl(X, y, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Compute Identity-Binomial Regression Hessian
//' @description Calculates the Hessian matrix for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_identity_binomial_regression_hessian_cpp(SEXP X_r,
															 SEXP y_r,
															 SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
	return constrained_binomial_hessian_cpp_impl(X, y, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Compute Weighted Identity-Binomial Regression Score
//' @description Calculates the weighted score vector for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the weighted score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_identity_binomial_regression_weighted_score_cpp(SEXP X_r,
                                                                    SEXP y_r,
                                                                    SEXP weights_r,
                                                                    SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
  return constrained_binomial_weighted_score_cpp_impl(X, y, weights, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Compute Weighted Identity-Binomial Regression Hessian
//' @description Calculates the weighted Hessian matrix for a binomial regression model with an identity link.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the weighted Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_identity_binomial_regression_weighted_hessian_cpp(SEXP X_r,
                                                                      SEXP y_r,
                                                                      SEXP weights_r,
                                                                      SEXP beta_sexp) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
    NumericVector beta_r(beta_sexp);
    Eigen::Map<const Eigen::VectorXd> beta(beta_r.begin(), beta_r.size());
  return constrained_binomial_weighted_hessian_cpp_impl(X, y, weights, beta, BinomialConstrainedLink::kIdentity);
}

//' @title Fast Log-Binomial Regression (C++)
//' @description High-performance log-binomial regression fitting using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients and fitted values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_log_binomial_regression_cpp(SEXP X_r,
                                      SEXP y_r,
                                      int maxit = 100,
                                      double tol = 1e-8,
                                      Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                      bool smart_cold_start = true,
                                      Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                      bool estimate_only = false) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kLog, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info, estimate_only);
}

//' @title Fast Log-Binomial Regression with Variance (C++)
//' @description Log-binomial regression with variance-covariance matrix and standard errors.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, vcov, and standard errors.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_log_binomial_regression_with_var_cpp(SEXP X_r,
                                               SEXP y_r,
                                               int j = 2,
                                               int maxit = 100,
                                               double tol = 1e-8,
                                               Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                               bool smart_cold_start = true,
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
  return fit_constrained_binomial_with_var_cpp_impl(X, y, BinomialConstrainedLink::kLog, j, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info);
}

//' @title Fast Weighted Log-Binomial Regression (C++)
//' @description High-performance weighted log-binomial regression fitting using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients and fitted values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_log_binomial_regression_weighted_cpp(SEXP X_r,
                                               SEXP y_r,
                                               SEXP weights_r,
                                               int maxit = 100,
                                               double tol = 1e-8,
                                               Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                               bool smart_cold_start = true,
                                               Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                               bool estimate_only = false) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
  return fit_constrained_binomial_weighted_cpp_impl(X, y, weights, BinomialConstrainedLink::kLog, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info, estimate_only);
}

//' @title Fast Identity-Binomial Regression (C++)
//' @description High-performance binomial regression with identity link using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients and fitted values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_identity_binomial_regression_cpp(SEXP X_r,
                                           SEXP y_r,
                                           int maxit = 100,
                                           double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = true,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
  return fit_constrained_binomial_cpp_impl(X, y, BinomialConstrainedLink::kIdentity, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info);
}

//' @title Fast Identity-Binomial Regression with Variance (C++)
//' @description Binomial regression with identity link, providing variance-covariance matrix and standard errors.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, vcov, and standard errors.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_identity_binomial_regression_with_var_cpp(SEXP X_r,
                                                    SEXP y_r,
                                                    int j = 2,
                                                    int maxit = 100,
                                                    double tol = 1e-8,
                                                    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                    bool smart_cold_start = true,
                                                    Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
  return fit_constrained_binomial_with_var_cpp_impl(X, y, BinomialConstrainedLink::kIdentity, j, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info);
}

//' @title Fast Weighted Identity-Binomial Regression (C++)
//' @description High-performance weighted binomial regression with identity link using Fisher scoring.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A nonnegative numeric vector of observation weights.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients and fitted values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_identity_binomial_regression_weighted_cpp(SEXP X_r,
                                                    SEXP y_r,
                                                    SEXP weights_r,
                                                    int maxit = 100,
                                                    double tol = 1e-8,
                                                    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                    bool smart_cold_start = true,
                                                    Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                                    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                                    bool estimate_only = false) {
    NumericMatrix X_mat(X_r);
    NumericVector y_vec(y_r);
    Eigen::Map<const Eigen::MatrixXd> X(X_mat.begin(), X_mat.nrow(), X_mat.ncol());
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    NumericVector weights_vec(weights_r);
    Eigen::Map<const Eigen::VectorXd> weights(weights_vec.begin(), weights_vec.size());
  return fit_constrained_binomial_weighted_cpp_impl(X, y, weights, BinomialConstrainedLink::kIdentity, maxit, tol, fixed_idx, fixed_values, warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info, estimate_only);
}
