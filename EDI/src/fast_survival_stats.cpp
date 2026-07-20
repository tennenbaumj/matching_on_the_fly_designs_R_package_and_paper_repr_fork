#include <RcppEigen.h>
#include <algorithm> // for std::sort
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace Eigen;

// File-local struct and helper for the parallel BRT kernel.
// Avoids all Rcpp wrap/unwrap in the per-draw inner loop so OpenMP can be used safely.
namespace {

struct SurvEntry { double time; int status; };

// Compute KM median or RMST for one sorted group; utimes/sprobs are reused across calls.
inline double km_stat_inline(SurvEntry* grp, int ng, bool do_rmst,
                              std::vector<double>& utimes, std::vector<double>& sprobs) {
    if (ng == 0) return NA_REAL;
    std::sort(grp, grp + ng, [](const SurvEntry& a, const SurvEntry& b){ return a.time < b.time; });
    utimes.clear(); sprobs.clear();
    utimes.push_back(0.0); sprobs.push_back(1.0);
    double sp = 1.0;
    for (int i = 0; i < ng; ) {
        double ct = grp[i].time;
        int ar = ng - i, ev = 0, j = i;
        while (j < ng && grp[j].time == ct) { if (grp[j].status == 1) ev++; j++; }
        if (ev > 0) { sp *= 1.0 - (double)ev / ar; utimes.push_back(ct); sprobs.push_back(sp); }
        i = j;
    }
    if (!do_rmst) {
        // Median: first time KM drops below 0.5 with linear interpolation
        const int sz = (int)sprobs.size();
        for (int i = 0; i < sz; ++i) {
            if (sprobs[i] < 0.5) {
                if (i > 0) {
                    double p1 = sprobs[i-1], p2 = sprobs[i], t1 = utimes[i-1], t2 = utimes[i];
                    return t1 + (t2 - t1) * (0.5 - p1) / (p2 - p1);
                }
                return utimes[i];
            }
        }
        return R_PosInf;
    } else {
        // RMST: area under the KM curve (trapezoidal integration)
        double rmst = 0.0;
        const int sz = (int)utimes.size();
        for (int i = 0; i + 1 < sz; ++i)
            rmst += sprobs[i] * (utimes[i+1] - utimes[i]);
        if (sz > 1)
            rmst += sprobs.back() * (grp[ng-1].time - utimes.back());
        return rmst;
    }
}

} // namespace

//' Calculates the median or restricted mean survival time for a single group
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param requested_stat A string, either "median" or "restricted_mean".
//' @return The calculated statistic.
//' @keywords internal
// [[Rcpp::export]]
double get_survival_stat_for_group(SEXP y_sexp, SEXP dead_sexp, std::string requested_stat) {
    NumericVector y_vec(y_sexp);
    IntegerVector dead_int(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> dead(dead_int.begin(), dead_int.size());
    // Combine y and dead into a data frame-like structure for sorting
    int n = y.size();
    if (n == 0) {
        return NA_REAL;
    }

    struct Subject {
        double time;
        int status;
    };

    std::vector<Subject> subjects(n);
    for (int i = 0; i < n; ++i) {
        subjects[i] = {y[i], dead[i]};
    }

    // Sort subjects by time
    std::sort(subjects.begin(), subjects.end(), [](const Subject& a, const Subject& b) {
        return a.time < b.time;
    });

    // Calculate Kaplan-Meier survival probability
    double survival_prob = 1.0;
    std::vector<double> unique_times;
    std::vector<double> survival_probs;

    unique_times.push_back(0.0);
    survival_probs.push_back(1.0);

    double last_unique_time = -1.0;
    int at_risk = n;
    int event_count_at_time = 0;
    int at_risk_at_time = n;

    for (int i = 0; i < n; ) {
        double current_time = subjects[i].time;
        at_risk_at_time = n - i;
        event_count_at_time = 0;

        int j = i;
        while (j < n && subjects[j].time == current_time) {
            if (subjects[j].status == 1) {
                event_count_at_time++;
            }
            j++;
        }

        if (event_count_at_time > 0) {
            survival_prob *= (1.0 - (double)event_count_at_time / at_risk_at_time);
            unique_times.push_back(current_time);
            survival_probs.push_back(survival_prob);
        }

        i = j;
    }


    if (requested_stat == "median") {
        for (size_t i = 0; i < survival_probs.size(); ++i) {
            if (survival_probs[i] < 0.5) {
                // simple linear interpolation
                if (i > 0){
                    double p1 = survival_probs[i-1];
                    double p2 = survival_probs[i];
                    double t1 = unique_times[i-1];
                    double t2 = unique_times[i];
                    return t1 + (t2 - t1) * (0.5 - p1) / (p2 - p1);
                } else {
                    return unique_times[i]; // should be 0 if it happens at first obs
                }
            }
        }
        return R_PosInf; // Median is beyond the last observation time
    } else if (requested_stat == "restricted_mean") {
        double restricted_mean = 0.0;
        for (size_t i = 0; i < unique_times.size() - 1; ++i) {
            restricted_mean += survival_probs[i] * (unique_times[i+1] - unique_times[i]);
        }
        // Add the last interval
        if (unique_times.size() > 1){
             restricted_mean += survival_probs.back() * (subjects.back().time - unique_times.back());
        }

        return restricted_mean;
    }

    return NA_REAL; // Should not be reached
}


//' Calculates the difference in a survival statistic (median or restricted mean)
//' between two groups (treatment vs control)
//'
//' @param y_sexp Numeric vector of survival times.
//' @param dead_sexp Integer vector of event indicators (1=event, 0=censored).
//' @param w_sexp Integer vector of treatment assignments (1=treatment, 0=control).
//' @param requested_stat A string, either "median" or "restricted_mean".
//' @return The difference in the statistic (treatment - control).
//' @keywords internal
// [[Rcpp::export]]
double get_survival_stat_diff(SEXP y_sexp, SEXP dead_sexp, SEXP w_sexp, std::string requested_stat) {
    NumericVector y_vec(y_sexp);
    IntegerVector dead_int(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> dead(dead_int.begin(), dead_int.size());
    IntegerVector w_int(w_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_int.begin(), w_int.size());
    std::vector<int> control_indices_std, treatment_indices_std;
    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0) {
            control_indices_std.push_back(i);
        } else {
            treatment_indices_std.push_back(i);
        }
    }

    std::vector<double> y_control_std, y_treatment_std;
    std::vector<int> dead_control_std, dead_treatment_std;

    for (int idx : control_indices_std) {
        y_control_std.push_back(y[idx]);
        dead_control_std.push_back(dead[idx]);
    }
    for (int idx : treatment_indices_std) {
        y_treatment_std.push_back(y[idx]);
        dead_treatment_std.push_back(dead[idx]);
    }

    double stat_control = get_survival_stat_for_group(wrap(y_control_std), wrap(dead_control_std), requested_stat);
    double stat_treatment = get_survival_stat_for_group(wrap(y_treatment_std), wrap(dead_treatment_std), requested_stat);

    if (R_IsNA(stat_treatment) || R_IsNA(stat_control)) {
        return NA_REAL;
    }

    return stat_treatment - stat_control;
}


//' Calculates standard variance using the formula from Uno et al
//'
//' \eqn{Var(RMST) = \sum_j A(t_j)^2 d_j / (n_j (n_j - d_j))}
//' where \eqn{A(t_j) = \int_{t_j}^{\tau} S(u) du} is the remaining area under the KM
//' curve from event time \eqn{t_j} to the last observation \eqn{\tau}.
//' Here \eqn{d_j} is the number of events at \eqn{t_j}, and \eqn{n_j}
//' is the number at risk just before \eqn{t_j}.
//' Terms where n_j == d_j are omitted: S drops to 0 there, so A(t_j) = 0 and the
//' contribution is 0 in the limit regardless of the undefined Greenwood denominator.
//'
//' @param y_sexp Numeric vector of survival times.
//' @param dead_sexp Integer vector of event indicators (1=event, 0=censored).
//' @return The standard error of the restricted mean.
//' @keywords internal
// [[Rcpp::export]]
double get_restricted_mean_se_for_group(SEXP y_sexp, SEXP dead_sexp) {
    NumericVector y_vec(y_sexp);
    IntegerVector dead_int(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> dead(dead_int.begin(), dead_int.size());
    int n = y.size();
    if (n == 0) return NA_REAL;

    struct Subject { double time; int status; };
    std::vector<Subject> subjects(n);
    for (int i = 0; i < n; ++i) subjects[i] = {y[i], dead[i]};
    std::sort(subjects.begin(), subjects.end(), [](const Subject& a, const Subject& b) {
        return a.time < b.time;
    });

    double tau = subjects.back().time;

    // Build KM event table: one entry per unique event time
    struct EventInfo { double time; double S_after; int n_j; int d_j; };
    std::vector<EventInfo> events;
    double S = 1.0;
    for (int i = 0; i < n; ) {
        double t = subjects[i].time;
        int n_at_risk = n - i;
        int d = 0;
        int j = i;
        while (j < n && subjects[j].time == t) {
            if (subjects[j].status == 1) d++;
            j++;
        }
        if (d > 0) {
            S *= (1.0 - (double)d / n_at_risk);
            events.push_back({t, S, n_at_risk, d});
        }
        i = j;
    }

    int K = (int)events.size();
    if (K == 0) return 0.0;  // no events: RMST equals tau with zero variance

    // Compute A(t_j) = integral_{t_j}^{tau} S(u) du via suffix sums.
    // S(u) = events[k].S_after for u in [events[k].time, events[k+1].time);
    // the final interval extends to tau.
    std::vector<double> A(K);
    A[K - 1] = events[K - 1].S_after * (tau - events[K - 1].time);
    for (int k = K - 2; k >= 0; --k) {
        A[k] = A[k + 1] + events[k].S_after * (events[k + 1].time - events[k].time);
    }

    // Var(RMST) = sum_j A(t_j)^2 * d_j / (n_j * (n_j - d_j)), skipping n_j == d_j
    double rmst_var = 0.0;
    for (int k = 0; k < K; ++k) {
        int nj = events[k].n_j;
        int dj = events[k].d_j;
        if (nj > dj) {
            rmst_var += A[k] * A[k] * (double)dj / ((double)nj * (nj - dj));
        }
    }

    return sqrt(rmst_var);
}

//' Calculates the standard error of the difference in restricted mean survival times
//'
//' @param y_sexp Numeric vector of survival times.
//' @param dead_sexp Integer vector of event indicators (1=event, 0=censored).
//' @param w_sexp Integer vector of treatment assignments (1=treatment, 0=control).
//' @return The standard error of the difference.
//' @keywords internal
// [[Rcpp::export]]
double get_restricted_mean_se_diff(SEXP y_sexp, SEXP dead_sexp, SEXP w_sexp) {
    NumericVector y_vec(y_sexp);
    IntegerVector dead_int(dead_sexp);
    Eigen::Map<const Eigen::VectorXd> y(y_vec.begin(), y_vec.size());
    Eigen::Map<const Eigen::VectorXi> dead(dead_int.begin(), dead_int.size());
    IntegerVector w_int(w_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_int.begin(), w_int.size());
    std::vector<int> control_indices_std, treatment_indices_std;
    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0) {
            control_indices_std.push_back(i);
        } else {
            treatment_indices_std.push_back(i);
        }
    }

    std::vector<double> y_control_std, y_treatment_std;
    std::vector<int> dead_control_std, dead_treatment_std;

    for (int idx : control_indices_std) {
        y_control_std.push_back(y[idx]);
        dead_control_std.push_back(dead[idx]);
    }
    for (int idx : treatment_indices_std) {
        y_treatment_std.push_back(y[idx]);
        dead_treatment_std.push_back(dead[idx]);
    }

    double se_control = get_restricted_mean_se_for_group(wrap(y_control_std), wrap(dead_control_std));
    double se_treatment = get_restricted_mean_se_for_group(wrap(y_treatment_std), wrap(dead_treatment_std));

    if (R_IsNA(se_treatment) || R_IsNA(se_control)) {
        return NA_REAL;
    }

    return sqrt(pow(se_treatment, 2.0) + pow(se_control, 2.0));
}


//' Parallel BRT kernel for KM-diff (median) and RMST-diff.
//' Each replicate resamples rows i_mat(.,b) and pairs them with assignment w_mat(.,b).
//' Sharp-null shift is multiplicative on treated times (exp(delta)). Uses an inline
//' pure-C++ KM calculator — no R objects inside the loop, so OpenMP is safe.
//' @param do_rmst TRUE for RMST-diff, FALSE for median (KM-diff).
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    bool do_rmst,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const double mult = std::exp(delta);

  const bool has_noise = noise_mat.isNotNull();
  NumericMatrix noise_m;
  const double* noise_ptr = nullptr;
  if (has_noise) {
    noise_m = NumericMatrix(noise_mat);
    noise_ptr = noise_m.begin();
  }

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
  {
    // Per-thread reusable buffers: avoids heap allocation in the hot loop
    std::vector<SurvEntry> y_t(n), y_c(n);
    std::vector<double> utimes_t, utimes_c, sprobs_t, sprobs_c;
    utimes_t.reserve(n); sprobs_t.reserve(n);
    utimes_c.reserve(n); sprobs_c.reserve(n);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      int nt = 0, nc = 0;
      for (int i = 0; i < n; ++i) {
        const int row0 = i_col[i] - 1;  // i_mat is 1-based
        double yv = y0_ptr[row0];
        if (has_noise) yv += noise_ptr[(size_t)b * n + i];
        SurvEntry e;
        e.status = dead_ptr[row0];
        if (w_col[i] == 1) {
          e.time = (delta != 0.0) ? yv * mult : yv;
          y_t[nt++] = e;
        } else {
          e.time = yv;
          y_c[nc++] = e;
        }
      }
      if (nt == 0 || nc == 0) continue;
      double stat_t = km_stat_inline(y_t.data(), nt, do_rmst, utimes_t, sprobs_t);
      double stat_c = km_stat_inline(y_c.data(), nc, do_rmst, utimes_c, sprobs_c);
      if (std::isfinite(stat_t) && std::isfinite(stat_c))
        res_ptr[b] = stat_t - stat_c;
    }
  }
  return wrap(results_vec);
}


// BRT variant of the KM survival-statistic difference (median or restricted mean):
// each replicate b resamples rows i_mat(., b) and pairs them with the fresh assignment
// w_mat(., b); the sharp-null shift is multiplicative on the treated times (delta on the
// log scale). Serial: get_survival_stat_for_group allocates R objects, which is not
// thread-safe; the win over the R paths is eliminating per-replicate object duplication.
// [[Rcpp::export]]
NumericVector compute_survival_stat_diff_rand_bootstrap_serial_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    std::string requested_stat) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  NumericVector results(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  const double mult = std::exp(delta);

  std::vector<double> y_t, y_c;
  std::vector<int> d_t, d_c;

  for (int b = 0; b < nsim; ++b) {
    const int* i_col = i_ptr + (size_t)b * n;
    const int* w_col = w_ptr + (size_t)b * n;
    y_t.clear(); y_c.clear(); d_t.clear(); d_c.clear();
    for (int i = 0; i < n; ++i) {
      const int row0 = i_col[i] - 1; // i_mat is 1-based
      if (w_col[i] == 1) {
        y_t.push_back(delta != 0.0 ? y0_ptr[row0] * mult : y0_ptr[row0]);
        d_t.push_back(dead_ptr[row0]);
      } else {
        y_c.push_back(y0_ptr[row0]);
        d_c.push_back(dead_ptr[row0]);
      }
    }
    if (y_t.empty() || y_c.empty()) continue;
    const double stat_t = get_survival_stat_for_group(wrap(y_t), wrap(d_t), requested_stat);
    const double stat_c = get_survival_stat_for_group(wrap(y_c), wrap(d_c), requested_stat);
    if (std::isfinite(stat_t) && std::isfinite(stat_c)) results[b] = stat_t - stat_c;
  }

  return results;
}
