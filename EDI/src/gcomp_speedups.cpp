// [[Rcpp::depends(RcppEigen)]]
#include "_helper_functions.h"
#include <unordered_map>
#include <utility>

using namespace Rcpp;

namespace {

inline double plogis_stable(double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

inline Eigen::ArrayXd plogis_array(const Eigen::ArrayXd& eta) {
  const Eigen::Array<bool, Eigen::Dynamic, 1> nonnegative = (eta >= 0.0);
  const Eigen::ArrayXd pos = 1.0 / (1.0 + (-eta).exp());
  const Eigen::ArrayXd neg_exp = eta.exp();
  const Eigen::ArrayXd neg = neg_exp / (1.0 + neg_exp);
  return nonnegative.select(pos, neg);
}

Eigen::MatrixXd cluster_meat(const Eigen::MatrixXd& X_fit,
                             const Eigen::VectorXd& resid,
                             const IntegerVector& cluster_id) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (cluster_id.size() != n) {
    stop("dimension mismatch in cluster_meat");
  }

  Eigen::MatrixXd meat = Eigen::MatrixXd::Zero(p, p);
  std::unordered_map<int, Eigen::VectorXd> cluster_scores;
  cluster_scores.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const int id = cluster_id[i];
    auto it = cluster_scores.find(id);
    if (it == cluster_scores.end()) {
      it = cluster_scores.emplace(id, Eigen::VectorXd::Zero(p)).first;
    }
    it->second.noalias() += X_fit.row(i).transpose() * resid[i];
  }
  for (const auto& entry : cluster_scores) {
    const auto& score_g = entry.second;
    meat.noalias() += score_g * score_g.transpose();
  }
  return meat;
}

}  // namespace

//' @title Fast G-Computation Point Estimate for Fractional Logit (C++)
//' @description Computes marginal mean difference under the fractional logit (quasi-binomial) model using G-computation.
//' @param X_fit Numeric matrix of predictors including intercept.
//' @param coef_hat Numeric vector of fitted coefficients.
//' @param j_treat 1-based column index of the treatment indicator in X_fit.
//' @return A list with elements \code{mean1}, \code{mean0}, and \code{md} (mean difference).
//' @export
//' @keywords internal
// [[Rcpp::export]]
List gcomp_fractional_logit_point_estimate_cpp(SEXP X_fit_sexp,
                                               SEXP coef_hat_sexp,
                                               int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  const int j_treat0 = j_treat - 1;

  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }

  Eigen::VectorXd eta = X_fit * coef_hat;
  Eigen::VectorXd eta_base = eta - coef_hat[j_treat0] * X_fit.col(j_treat0);

  Eigen::ArrayXd risk1_arr = plogis_array_safe((eta_base.array() + coef_hat[j_treat0]));
  Eigen::ArrayXd risk0_arr = plogis_array_safe(eta_base.array());

  double mean1 = risk1_arr.mean();
  double mean0 = risk0_arr.mean();

  return List::create(
    _["mean1"] = mean1,
    _["mean0"] = mean0,
    _["md"] = mean1 - mean0
  );
}

//' @title Fast G-Computation Point Estimate for Logistic Regression (C++)
//' @description Computes marginal risk difference and risk ratio under the logistic model using G-computation.
//' @param X_fit Numeric matrix of predictors including intercept.
//' @param coef_hat Numeric vector of fitted coefficients.
//' @param j_treat 1-based column index of the treatment indicator in X_fit.
//' @return A list with elements \code{mean1}, \code{mean0}, and \code{md} (mean difference).
//' @export
//' @keywords internal
// [[Rcpp::export]]
List gcomp_logistic_point_estimate_cpp(SEXP X_fit_sexp,
                                        SEXP coef_hat_sexp,
                                        int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  return gcomp_fractional_logit_point_estimate_cpp(X_fit_sexp, coef_hat_sexp, j_treat);
}

namespace {

struct LogisticPostFitResult {
  Eigen::MatrixXd vcov;
  Eigen::VectorXd std_err;
  Eigen::VectorXd z_vals;
  double risk1;
  double risk0;
  double rd;
  double se_rd;
  double log_rr;
  double rr;
  double se_log_rr;
};

LogisticPostFitResult compute_gcomp_logistic_post_fit(SEXP X_fit_sexp,
                                                       SEXP y_sexp,
                                                       SEXP coef_hat_sexp,
                                                       SEXP mu_hat_sexp,
                                                       int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericVector mu_hat_r(mu_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> mu_hat(mu_hat_r.begin(), mu_hat_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  const int j_treat0 = j_treat - 1;

  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n) {
    stop("dimension mismatch in gcomp_logistic_post_fit_cpp");
  }

  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
  }
  const Eigen::VectorXd W = (mu_hat.array() * (1.0 - mu_hat.array())).matrix();

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, W);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }

  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd score_resid = y - mu_hat;
  const Eigen::VectorXd resid_sq = score_resid.array().square().matrix();
  Eigen::MatrixXd meat = weighted_crossprod(X_fit, resid_sq);
  Eigen::MatrixXd vcov_robust = bread * meat * bread;
  vcov_robust = 0.5 * (vcov_robust + vcov_robust.transpose());

  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      if (!R_finite(vcov_robust(j, k))) {
        stop("non-finite robust covariance");
      }
    }
  }

  const double ssq_treat = vcov_robust(j_treat0, j_treat0);
  if (!R_finite(ssq_treat) || ssq_treat <= 0.0) {
    stop("non-positive treatment variance");
  }

  const Eigen::VectorXd eta = X_fit * coef_hat;
  const Eigen::VectorXd eta_base = eta - coef_hat[j_treat0] * X_fit.col(j_treat0);

  const Eigen::ArrayXd risk1_arr = plogis_array(eta_base.array() + coef_hat[j_treat0]);
  const Eigen::ArrayXd risk0_arr = plogis_array(eta_base.array());
  const Eigen::VectorXd risk1_i = risk1_arr.matrix();
  const Eigen::VectorXd risk0_i = risk0_arr.matrix();

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  const double inv_n = 1.0 / static_cast<double>(n);
  const Eigen::VectorXd wt1 = (risk1_arr * (1.0 - risk1_arr) * inv_n).matrix();
  const Eigen::VectorXd wt0 = (risk0_arr * (1.0 - risk0_arr) * inv_n).matrix();
  Eigen::VectorXd grad1 = X_fit.transpose() * wt1;
  Eigen::VectorXd grad0 = X_fit.transpose() * wt0;
  grad1[j_treat0] = wt1.sum();
  grad0[j_treat0] = 0.0;

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const Eigen::VectorXd bread_grad_rd = ldlt.solve(grad_rd);
  const double var_rd = (bread_grad_rd.transpose() * meat * bread_grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const Eigen::VectorXd bread_grad_log_rr = ldlt.solve(grad_log_rr);
    const double var_log_rr = (bread_grad_log_rr.transpose() * meat * bread_grad_log_rr)(0, 0);
    if (R_finite(var_log_rr) && var_log_rr >= 0.0) {
      se_log_rr = std::sqrt(var_log_rr);
    }
  }

  Eigen::VectorXd std_err(p);
  Eigen::VectorXd z_vals(p);
  for (int j = 0; j < p; ++j) {
    const double var_j = vcov_robust(j, j);
    std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
    z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? coef_hat[j] / std_err[j] : NA_REAL;
  }

  return {
    std::move(vcov_robust),
    std::move(std_err),
    std::move(z_vals),
    risk1,
    risk0,
    rd,
    se_rd,
    log_rr,
    rr,
    se_log_rr
  };
}

}  // namespace

// [[Rcpp::export]]
List gcomp_logistic_post_fit_cpp(SEXP X_fit_sexp,
                                 SEXP y_sexp,
                                 SEXP coef_hat_sexp,
                                 SEXP mu_hat_sexp,
                                 int j_treat) {
  LogisticPostFitResult result = compute_gcomp_logistic_post_fit(
    X_fit_sexp, y_sexp, coef_hat_sexp, mu_hat_sexp, j_treat
  );
  return List::create(
    _["vcov"] = result.vcov,
    _["std_err"] = result.std_err,
    _["z_vals"] = result.z_vals,
    _["risk1"] = result.risk1,
    _["risk0"] = result.risk0,
    _["rd"] = result.rd,
    _["se_rd"] = result.se_rd,
    _["log_rr"] = result.log_rr,
    _["rr"] = result.rr,
    _["se_log_rr"] = result.se_log_rr
  );
}

// [[Rcpp::export]]
List gcomp_fractional_logit_post_fit_cpp(SEXP X_fit_sexp,
                                         SEXP y_sexp,
                                         SEXP coef_hat_sexp,
                                         SEXP mu_hat_sexp,
                                         int j_treat) {
  LogisticPostFitResult result = compute_gcomp_logistic_post_fit(
    X_fit_sexp, y_sexp, coef_hat_sexp, mu_hat_sexp, j_treat
  );
  return List::create(
    _["vcov"] = result.vcov,
    _["std_err"] = result.std_err,
    _["z_vals"] = result.z_vals,
    _["mean1"] = result.risk1,
    _["mean0"] = result.risk0,
    _["md"] = result.rd,
    _["se_md"] = result.se_rd
  );
}

// [[Rcpp::export]]
List gcomp_logistic_cluster_post_fit_cpp(SEXP X_fit_sexp,
                                         SEXP y_sexp,
                                         SEXP coef_hat_sexp,
                                         SEXP mu_hat_sexp,
                                         const IntegerVector& cluster_id,
                                         int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericVector mu_hat_r(mu_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> mu_hat(mu_hat_r.begin(), mu_hat_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  const int j_treat0 = j_treat - 1;

  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n) {
    stop("dimension mismatch in gcomp_logistic_cluster_post_fit_cpp");
  }

  for (int i = 0; i < n; ++i) {
    const double mu_i = mu_hat[i];
    if (!R_finite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0) {
      stop("non-finite or boundary fitted values");
    }
  }
  const Eigen::VectorXd W = (mu_hat.array() * (1.0 - mu_hat.array())).matrix();

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, W);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }

  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd score_resid = y - mu_hat;
  Eigen::MatrixXd meat = cluster_meat(X_fit, score_resid, cluster_id);
  Eigen::MatrixXd vcov_robust = bread * meat * bread;
  vcov_robust = 0.5 * (vcov_robust + vcov_robust.transpose());

  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      if (!R_finite(vcov_robust(j, k))) {
        stop("non-finite robust covariance");
      }
    }
  }

  const double ssq_treat = vcov_robust(j_treat0, j_treat0);
  if (!R_finite(ssq_treat) || ssq_treat <= 0.0) {
    stop("non-positive treatment variance");
  }

  const Eigen::VectorXd eta = X_fit * coef_hat;
  const Eigen::VectorXd eta_base = eta - coef_hat[j_treat0] * X_fit.col(j_treat0);

  const Eigen::ArrayXd risk1_arr = plogis_array(eta_base.array() + coef_hat[j_treat0]);
  const Eigen::ArrayXd risk0_arr = plogis_array(eta_base.array());
  const Eigen::VectorXd risk1_i = risk1_arr.matrix();
  const Eigen::VectorXd risk0_i = risk0_arr.matrix();

  const double risk1 = risk1_i.mean();
  const double risk0 = risk0_i.mean();

  const double inv_n = 1.0 / static_cast<double>(n);
  const Eigen::VectorXd wt1 = (risk1_arr * (1.0 - risk1_arr) * inv_n).matrix();
  const Eigen::VectorXd wt0 = (risk0_arr * (1.0 - risk0_arr) * inv_n).matrix();
  Eigen::VectorXd grad1 = X_fit.transpose() * wt1;
  Eigen::VectorXd grad0 = X_fit.transpose() * wt0;
  grad1[j_treat0] = wt1.sum();
  grad0[j_treat0] = 0.0;

  const double rd = risk1 - risk0;
  const Eigen::VectorXd grad_rd = grad1 - grad0;
  const Eigen::VectorXd bread_grad_rd = ldlt.solve(grad_rd);
  const double var_rd = (bread_grad_rd.transpose() * meat * bread_grad_rd)(0, 0);
  const double se_rd = (R_finite(var_rd) && var_rd >= 0.0) ? std::sqrt(var_rd) : NA_REAL;

  double log_rr = NA_REAL;
  double rr = NA_REAL;
  double se_log_rr = NA_REAL;
  if (risk1 > 0.0 && risk0 > 0.0) {
    log_rr = std::log(risk1) - std::log(risk0);
    rr = std::exp(log_rr);
    const Eigen::VectorXd grad_log_rr = grad1 / risk1 - grad0 / risk0;
    const Eigen::VectorXd bread_grad_log_rr = ldlt.solve(grad_log_rr);
    const double var_log_rr = (bread_grad_log_rr.transpose() * meat * bread_grad_log_rr)(0, 0);
    if (R_finite(var_log_rr) && var_log_rr >= 0.0) {
      se_log_rr = std::sqrt(var_log_rr);
    }
  }

  Eigen::VectorXd std_err(p);
  Eigen::VectorXd z_vals(p);
  for (int j = 0; j < p; ++j) {
    const double var_j = vcov_robust(j, j);
    std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
    z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? coef_hat[j] / std_err[j] : NA_REAL;
  }

  return List::create(
    _["vcov"] = vcov_robust,
    _["std_err"] = std_err,
    _["z_vals"] = z_vals,
    _["risk1"] = risk1,
    _["risk0"] = risk0,
    _["rd"] = rd,
    _["se_rd"] = se_rd,
    _["log_rr"] = log_rr,
    _["rr"] = rr,
    _["se_log_rr"] = se_log_rr
  );
}

// [[Rcpp::export]]
List gcomp_ordinal_proportional_odds_post_fit_cpp(SEXP X_fit_sexp,
                                                  SEXP coef_hat_sexp,
                                                  SEXP alpha_hat_sexp,
                                                  int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericVector alpha_hat_r(alpha_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> alpha_hat(alpha_hat_r.begin(), alpha_hat_r.size());
  const int n = X_fit.rows();
  const int K_minus_1 = alpha_hat.size();
  const int j_treat0 = j_treat - 1;

  Eigen::VectorXd eta_base = X_fit * coef_hat - coef_hat[j_treat0] * X_fit.col(j_treat0);
  Eigen::VectorXd eta1 = eta_base.array() + coef_hat[j_treat0];
  Eigen::VectorXd eta0 = eta_base;

  auto compute_mean = [&](const Eigen::VectorXd& eta_vec) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      double m = 1.0;
      for (int k = 0; k < K_minus_1; ++k) {
        m += 1.0 - plogis_stable(alpha_hat[k] - eta_vec[i]);
      }
      sum += m;
    }
    return sum / n;
  };

  double mean1 = compute_mean(eta1);
  double mean0 = compute_mean(eta0);

  return List::create(
    _["mean1"] = mean1,
    _["mean0"] = mean0,
    _["md"] = mean1 - mean0
  );
}
