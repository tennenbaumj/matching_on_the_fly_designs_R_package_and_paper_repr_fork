#include "_helper_functions.h"
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

struct LogChooseTable {
  std::vector<int> offsets;
  std::vector<double> values;
};

std::vector<double> build_log_factorials(int n) {
  std::vector<double> log_factorial(static_cast<std::size_t>(n + 1), 0.0);
  for (int i = 2; i <= n; ++i) {
    log_factorial[static_cast<std::size_t>(i)] =
      log_factorial[static_cast<std::size_t>(i - 1)] + std::log(static_cast<double>(i));
  }
  return log_factorial;
}

inline double log_choose_from_factorials(const std::vector<double>& log_factorial,
                                         int n,
                                         int k) {
  if (k < 0 || k > n) return R_NegInf;
  return log_factorial[static_cast<std::size_t>(n)] -
    log_factorial[static_cast<std::size_t>(k)] -
    log_factorial[static_cast<std::size_t>(n - k)];
}

LogChooseTable build_log_choose_table(const std::vector<int>& total_counts,
                                      const std::vector<double>& log_factorial) {
  const int K = static_cast<int>(total_counts.size());
  LogChooseTable table;
  table.offsets.resize(static_cast<std::size_t>(K + 1), 0);
  int total_size = 0;
  for (int k = 0; k < K; ++k) {
    table.offsets[static_cast<std::size_t>(k)] = total_size;
    total_size += total_counts[static_cast<std::size_t>(k)] + 1;
  }
  table.offsets[static_cast<std::size_t>(K)] = total_size;
  table.values.resize(static_cast<std::size_t>(total_size));

  for (int k = 0; k < K; ++k) {
    const int nk = total_counts[static_cast<std::size_t>(k)];
    const int offset = table.offsets[static_cast<std::size_t>(k)];
    for (int tk = 0; tk <= nk; ++tk) {
      table.values[static_cast<std::size_t>(offset + tk)] =
        log_choose_from_factorials(log_factorial, nk, tk);
    }
  }
  return table;
}

inline double log_choose_from_table(const LogChooseTable& table, int k, int t) {
  return table.values[static_cast<std::size_t>(table.offsets[static_cast<std::size_t>(k)] + t)];
}

int compute_jt_statistic2_from_counts(const std::vector<int>& treat_counts,
                                      const std::vector<int>& total_counts) {
  const int K = static_cast<int>(total_counts.size());
  int treated_seen = 0;
  int total_seen = 0;
  int stat2 = 0;

  for (int k = 0; k < K; ++k) {
    const int tk = treat_counts[static_cast<std::size_t>(k)];
    const int nk = total_counts[static_cast<std::size_t>(k)];
    const int lower_controls = total_seen - treated_seen;
    stat2 += tk * (2 * lower_controls + (nk - tk));
    treated_seen += tk;
    total_seen += nk;
  }
  return stat2;
}

void recurse_jt_distribution(int idx,
                             int remaining_treated,
                             const std::vector<int>& total_counts,
                             const LogChooseTable& log_choose_table,
                             int treated_seen,
                             int total_seen,
                             int stat2_so_far,
                             double log_weight,
                             std::vector<double>& stat_prob,
                             std::vector<unsigned char>& stat_active,
                             std::vector<int>& active_stats,
                             double log_norm) {
  const int K = static_cast<int>(total_counts.size());
  const int nk = total_counts[static_cast<std::size_t>(idx)];
  const int lower_controls = total_seen - treated_seen;

  if (idx == K - 1) {
    if (remaining_treated < 0 || remaining_treated > nk) return;
    const int stat2 = stat2_so_far +
      remaining_treated * (2 * lower_controls + (nk - remaining_treated));
    const double lw = log_weight + log_choose_from_table(log_choose_table, idx, remaining_treated) - log_norm;
    double& prob = stat_prob[static_cast<std::size_t>(stat2)];
    if (!stat_active[static_cast<std::size_t>(stat2)]) {
      stat_active[static_cast<std::size_t>(stat2)] = 1U;
      active_stats.push_back(stat2);
    }
    prob += std::exp(lw);
    return;
  }

  const int max_take = std::min(nk, remaining_treated);
  for (int tk = 0; tk <= max_take; ++tk) {
    const int stat2_next = stat2_so_far +
      tk * (2 * lower_controls + (nk - tk));
    recurse_jt_distribution(
      idx + 1,
      remaining_treated - tk,
      total_counts,
      log_choose_table,
      treated_seen + tk,
      total_seen + nk,
      stat2_next,
      log_weight + log_choose_from_table(log_choose_table, idx, tk),
      stat_prob,
      stat_active,
      active_stats,
      log_norm
    );
  }
}

}  // namespace

// [[Rcpp::export]]
List exact_jonckheere_terpstra_pval_cpp(SEXP y_sexp,
                                        SEXP w_sexp) {
  IntegerVector y_r(y_sexp);
  IntegerVector w_r(w_sexp);
  Eigen::Map<const Eigen::VectorXi> y(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXi> w(w_r.begin(), w_r.size());
  const int n = y.size();
  if (w.size() != n) stop("dimension mismatch in exact_jonckheere_terpstra_pval_cpp");
  if (n == 0) stop("empty input in exact_jonckheere_terpstra_pval_cpp");

  const int* y_ptr = y.data();
  const int* w_ptr = w.data();

  int n_treat = 0;
  int n_control = 0;
  for (int i = 0; i < n; ++i) {
    if (y_ptr[i] == NA_INTEGER || w_ptr[i] == NA_INTEGER) {
      stop("missing values are not allowed in exact_jonckheere_terpstra_pval_cpp");
    }
    if (w_ptr[i] == 1) {
      ++n_treat;
    } else if (w_ptr[i] == 0) {
      ++n_control;
    } else {
      stop("treatment assignments must be 0/1 in exact_jonckheere_terpstra_pval_cpp");
    }
  }
  if (n_treat == 0 || n_control == 0) {
    stop("both treatment arms must be present in exact_jonckheere_terpstra_pval_cpp");
  }

  std::vector<int> levels(n);
  for(int i=0; i<n; ++i) levels[i] = y_ptr[i];
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  const int K = static_cast<int>(levels.size());

  std::vector<int> total_counts(static_cast<std::size_t>(K), 0);
  std::vector<int> treat_counts_obs(static_cast<std::size_t>(K), 0);

  for (int i = 0; i < n; ++i) {
    const int yi = y_ptr[i];
    const int k = static_cast<int>(std::lower_bound(levels.begin(), levels.end(), yi) - levels.begin());
    ++total_counts[static_cast<std::size_t>(k)];
    if (w_ptr[i] == 1) ++treat_counts_obs[static_cast<std::size_t>(k)];
  }

  const int stat2_obs = compute_jt_statistic2_from_counts(treat_counts_obs, total_counts);
  const std::vector<double> log_factorial = build_log_factorials(n);
  const double log_norm = log_choose_from_factorials(log_factorial, n, n_treat);
  const LogChooseTable log_choose_table = build_log_choose_table(total_counts, log_factorial);
  const int max_stat2 = 2 * n_treat * n_control;
  thread_local std::vector<double> stat_prob;
  thread_local std::vector<unsigned char> stat_active;
  thread_local std::vector<int> active_stats;
  for (int stat2 : active_stats) {
    stat_prob[static_cast<std::size_t>(stat2)] = 0.0;
    stat_active[static_cast<std::size_t>(stat2)] = 0U;
  }
  active_stats.clear();
  if (stat_prob.size() < static_cast<std::size_t>(max_stat2 + 1)) {
    stat_prob.resize(static_cast<std::size_t>(max_stat2 + 1), 0.0);
  }
  if (stat_active.size() < static_cast<std::size_t>(max_stat2 + 1)) {
    stat_active.resize(static_cast<std::size_t>(max_stat2 + 1), 0U);
  }
  active_stats.reserve(static_cast<std::size_t>(std::min(max_stat2 + 1, 1024)));
  recurse_jt_distribution(0, n_treat, total_counts, log_choose_table, 0, 0, 0, 0.0, stat_prob, stat_active, active_stats, log_norm);

  double p_lower = 0.0;
  double p_upper = 0.0;
  for (int stat2 : active_stats) {
    const double prob = stat_prob[static_cast<std::size_t>(stat2)];
    if (stat2 <= stat2_obs) p_lower += prob;
    if (stat2 >= stat2_obs) p_upper += prob;
  }

  const double p_exact = std::min(1.0, 2.0 * std::min(p_lower, p_upper));
  const double superiority = static_cast<double>(stat2_obs) / (2.0 * n_treat * n_control);

  return List::create(
    _["stat2"] = stat2_obs,
    _["n_treat"] = n_treat,
    _["n_control"] = n_control,
    _["superiority"] = superiority,
    _["p_lower"] = p_lower,
    _["p_upper"] = p_upper,
    _["p_exact"] = p_exact
  );
}

// BRT variant of the two-group superiority statistic P(Y_T > Y_C) + 0.5 P(Y_T = Y_C) - 0.5
// (the Jonckheere-Terpstra effect estimate). Each replicate b resamples rows i_mat(., b)
// and pairs them with the fresh assignment w_mat(., b); no sharp-null shift is supported
// (ordinal responses; the R hook declines delta != 0).
// [[Rcpp::export]]
NumericVector compute_jt_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
  if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
  {
    std::vector<double> y_t, y_c;

#pragma omp for schedule(static)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      y_t.clear(); y_c.clear();
      for (int i = 0; i < n; ++i) {
        const double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
        if (!std::isfinite(yv)) continue;
        if (w_col[i] == 1) y_t.push_back(yv); else y_c.push_back(yv);
      }
      if (y_t.empty() || y_c.empty()) continue;
      std::sort(y_c.begin(), y_c.end());
      double num = 0.0;
      for (double yt : y_t) {
        const std::size_t less = std::lower_bound(y_c.begin(), y_c.end(), yt) - y_c.begin();
        const std::size_t leq = std::upper_bound(y_c.begin(), y_c.end(), yt) - y_c.begin();
        num += static_cast<double>(less) + 0.5 * static_cast<double>(leq - less);
      }
      res_ptr[b] = num / (static_cast<double>(y_t.size()) * static_cast<double>(y_c.size())) - 0.5;
    }
  }

  return wrap(results_vec);
}
