// [[Rcpp::depends(RcppEigen)]]
#include "_helper_functions.h"
#include <limits>
#include <unordered_map>

using namespace Rcpp;

namespace {

bool all_finite_mat(const Eigen::MatrixXd& M) {
  for (int j = 0; j < M.cols(); ++j) {
    for (int i = 0; i < M.rows(); ++i) {
      if (!R_finite(M(i, j))) {
        return false;
      }
    }
  }
  return true;
}

bool all_finite_vec(const Eigen::VectorXd& x) {
  for (int i = 0; i < x.size(); ++i) {
    if (!R_finite(x[i])) {
      return false;
    }
  }
  return true;
}

List summarize_with_vcov(const Eigen::VectorXd& coef_hat,
                         const Eigen::MatrixXd& vcov,
                         int j_treat) {
  const int p = coef_hat.size();
  const int j_treat0 = j_treat - 1;
  if (j_treat0 < 0 || j_treat0 >= p) {
    stop("treatment column index is out of bounds");
  }
  if (!all_finite_mat(vcov)) {
    stop("non-finite covariance matrix");
  }

  Eigen::VectorXd std_err(p);
  Eigen::VectorXd z_vals(p);
  for (int j = 0; j < p; ++j) {
    const double var_j = vcov(j, j);
    std_err[j] = (R_finite(var_j) && var_j >= 0.0) ? std::sqrt(var_j) : NA_REAL;
    z_vals[j] = (R_finite(std_err[j]) && std_err[j] > 0.0) ? coef_hat[j] / std_err[j] : NA_REAL;
  }

  const double ssq_hat = vcov(j_treat0, j_treat0);
  const double beta_hat = coef_hat[j_treat0];

  return List::create(
    _["beta_hat"] = beta_hat,
    _["ssq_hat"] = ssq_hat,
    _["se"] = (R_finite(ssq_hat) && ssq_hat >= 0.0) ? std::sqrt(ssq_hat) : NA_REAL,
    _["vcov"] = vcov,
    _["std_err"] = std_err,
    _["z_vals"] = z_vals
  );
}

Eigen::MatrixXd cluster_meat(const Eigen::MatrixXd& X_fit,
                             const Eigen::VectorXd& resid,
                             const IntegerVector& cluster_id) {
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (cluster_id.size() != n) {
    stop("dimension mismatch in cluster_meat");
  }

  std::unordered_map<int, int> cluster_lookup;
  std::vector<Eigen::VectorXd> cluster_scores;
  cluster_scores.reserve(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    const int id = cluster_id[i];
    auto it = cluster_lookup.find(id);
    int pos;
    if (it == cluster_lookup.end()) {
      pos = static_cast<int>(cluster_scores.size());
      cluster_lookup.emplace(id, pos);
      cluster_scores.push_back(Eigen::VectorXd::Zero(p));
    } else {
      pos = it->second;
    }
    cluster_scores[static_cast<std::size_t>(pos)].noalias() += X_fit.row(i).transpose() * resid[i];
  }

  Eigen::MatrixXd meat = Eigen::MatrixXd::Zero(p, p);
  for (const auto& score_g : cluster_scores) {
    meat.noalias() += score_g * score_g.transpose();
  }
  return meat;
}

}  // namespace

// [[Rcpp::export]]
List ols_hc2_setup_cpp(SEXP X_fit_sexp) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (!all_finite_mat(X_fit)) {
    stop("non-finite design matrix");
  }

  const Eigen::MatrixXd XtX = X_fit.transpose() * X_fit;
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'X");
  }
  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success || !all_finite_mat(bread)) {
    stop("failed to invert X'X");
  }

  const Eigen::VectorXd hat = (X_fit * bread).cwiseProduct(X_fit).rowwise().sum();

  return List::create(
    _["bread"] = bread,
    _["hat"] = hat
  );
}

// [[Rcpp::export]]
List ols_hc2_post_fit_precomputed_cpp(SEXP X_fit_sexp,
                                      SEXP y_sexp,
                                      SEXP coef_hat_sexp,
                                      SEXP bread_sexp,
                                      SEXP hat_sexp,
                                      int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericMatrix bread_r(bread_sexp);
  Rcpp::NumericVector hat_r(hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::MatrixXd> bread(bread_r.begin(), bread_r.nrow(), bread_r.ncol());
  Eigen::Map<const Eigen::VectorXd> hat(hat_r.begin(), hat_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (y.size() != n || coef_hat.size() != p || bread.rows() != p || bread.cols() != p || hat.size() != n) {
    stop("dimension mismatch in ols_hc2_post_fit_precomputed_cpp");
  }
  if (!all_finite_vec(coef_hat) || !all_finite_mat(bread) || !all_finite_vec(hat)) {
    stop("non-finite inputs");
  }

  const Eigen::VectorXd resid = y - X_fit * coef_hat;
  if (!all_finite_vec(resid)) {
    stop("non-finite residuals");
  }

  Eigen::VectorXd omega(n);
  for (int i = 0; i < n; ++i) {
    omega[i] = resid[i] * resid[i] / std::max(1.0 - hat[i], std::numeric_limits<double>::epsilon());
  }

  Eigen::MatrixXd meat = weighted_crossprod(X_fit, omega);
  Eigen::MatrixXd vcov = bread * meat * bread;
  vcov = 0.5 * (vcov + vcov.transpose());
  return summarize_with_vcov(coef_hat, vcov, j_treat);
}

// [[Rcpp::export]]
List ols_hc2_post_fit_cpp(SEXP X_fit_sexp,
                          SEXP y_sexp,
                          SEXP coef_hat_sexp,
                          int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  List setup = ols_hc2_setup_cpp(X_fit_sexp);
  return ols_hc2_post_fit_precomputed_cpp(
    X_fit_sexp,
    y_sexp,
    coef_hat_sexp,
    setup["bread"],
    setup["hat"],
    j_treat
  );
}

// [[Rcpp::export]]
List glm_sandwich_post_fit_cpp(SEXP X_fit_sexp,
                               SEXP y_sexp,
                               SEXP coef_hat_sexp,
                               SEXP mu_hat_sexp,
                               SEXP working_weights_sexp,
                               int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericVector mu_hat_r(mu_hat_sexp);
  Rcpp::NumericVector working_weights_r(working_weights_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> mu_hat(mu_hat_r.begin(), mu_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> working_weights(working_weights_r.begin(), working_weights_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n || working_weights.size() != n) {
    stop("dimension mismatch in glm_sandwich_post_fit_cpp");
  }
  if (!all_finite_vec(coef_hat) || !all_finite_vec(mu_hat) || !all_finite_vec(working_weights)) {
    stop("non-finite inputs");
  }

  for (int i = 0; i < n; ++i) {
    if (working_weights[i] <= 0.0) {
      stop("non-positive working weights");
    }
  }

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, working_weights);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }
  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success || !all_finite_mat(bread)) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd resid = y - mu_hat;
  const Eigen::VectorXd resid_sq = resid.array().square().matrix();
  Eigen::MatrixXd meat = weighted_crossprod(X_fit, resid_sq);
  Eigen::MatrixXd vcov = bread * meat * bread;
  vcov = 0.5 * (vcov + vcov.transpose());
  return summarize_with_vcov(coef_hat, vcov, j_treat);
}

// [[Rcpp::export]]
List glm_cluster_sandwich_post_fit_cpp(SEXP X_fit_sexp,
                                       SEXP y_sexp,
                                       SEXP coef_hat_sexp,
                                       SEXP mu_hat_sexp,
                                       SEXP working_weights_sexp,
                                       const IntegerVector& cluster_id,
                                       int j_treat) {
  Rcpp::NumericMatrix X_fit_r(X_fit_sexp);
  Rcpp::NumericVector y_r(y_sexp);
  Rcpp::NumericVector coef_hat_r(coef_hat_sexp);
  Rcpp::NumericVector mu_hat_r(mu_hat_sexp);
  Rcpp::NumericVector working_weights_r(working_weights_sexp);
  Eigen::Map<const Eigen::MatrixXd> X_fit(X_fit_r.begin(), X_fit_r.nrow(), X_fit_r.ncol());
  Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> coef_hat(coef_hat_r.begin(), coef_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> mu_hat(mu_hat_r.begin(), mu_hat_r.size());
  Eigen::Map<const Eigen::VectorXd> working_weights(working_weights_r.begin(), working_weights_r.size());
  const int n = X_fit.rows();
  const int p = X_fit.cols();
  if (y.size() != n || coef_hat.size() != p || mu_hat.size() != n || working_weights.size() != n) {
    stop("dimension mismatch in glm_cluster_sandwich_post_fit_cpp");
  }
  if (!all_finite_vec(coef_hat) || !all_finite_vec(mu_hat) || !all_finite_vec(working_weights)) {
    stop("non-finite inputs");
  }

  for (int i = 0; i < n; ++i) {
    if (working_weights[i] <= 0.0) {
      stop("non-positive working weights");
    }
  }

  const Eigen::MatrixXd XtWX = weighted_crossprod(X_fit, working_weights);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    stop("failed to factorize X'WX");
  }
  const Eigen::MatrixXd bread = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  if (ldlt.info() != Eigen::Success || !all_finite_mat(bread)) {
    stop("failed to invert X'WX");
  }

  const Eigen::VectorXd resid = y - mu_hat;
  Eigen::MatrixXd meat = cluster_meat(X_fit, resid, cluster_id);
  Eigen::MatrixXd vcov = bread * meat * bread;
  vcov = 0.5 * (vcov + vcov.transpose());
  return summarize_with_vcov(coef_hat, vcov, j_treat);
}
