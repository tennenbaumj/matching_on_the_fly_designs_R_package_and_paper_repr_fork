// Helpers for DesignFixedRerandomization.
//
// complete_randomization_forced_balanced_cpp : r balanced draws (n/2 treated), returns r×n
// complete_randomization_imbalanced_cpp      : r draws with nT treated, returns r×n
// compute_objective_vals_cpp                 : per-row objective value for an r×n indicTs matrix
// rerandomization_search_cpp                 : parallel rejection sampler; returns n×k (k<=r)
//   matrix of allocations whose objective value is <= cutoff, using at most max_draws total
//   random draws. Uses RcppEigen M-matrix trick for SIMD-optimized objective evaluation and
//   OpenMP for parallel draw generation.

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <random>
#include <chrono>
#include <cmath>
#include <atomic>
#include <numeric>
#include <algorithm>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::LLT;
using Eigen::Success;

// ── Helpers shared by the simple randomization functions ────────────────────

static std::default_random_engine make_rng(int seed) {
    unsigned int s = (seed == NA_INTEGER)
        ? static_cast<unsigned int>(
              std::chrono::system_clock::now().time_since_epoch().count())
        : static_cast<unsigned int>(seed);
    return std::default_random_engine(s);
}

static void shuffle_vec(std::vector<int>& v, std::default_random_engine& rng) {
    for (int i = static_cast<int>(v.size()) - 1; i > 0; i--) {
        std::uniform_int_distribution<int> d(0, i);
        std::swap(v[i], v[d(rng)]);
    }
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix complete_randomization_forced_balanced_cpp(int n, int r, int seed) {
    std::default_random_engine rng = make_rng(seed);
    Rcpp::IntegerMatrix out(r, n);
    std::vector<int> v(n, 0);
    int nT = n / 2;
    for (int i = 0; i < nT; i++) v[i] = 1;
    for (int row = 0; row < r; row++) {
        shuffle_vec(v, rng);
        for (int j = 0; j < n; j++) out(row, j) = v[j];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix complete_randomization_imbalanced_cpp(int n, int nT, int r, int seed) {
    std::default_random_engine rng = make_rng(seed);
    Rcpp::IntegerMatrix out(r, n);
    std::vector<int> v(n, 0);
    for (int i = 0; i < nT; i++) v[i] = 1;
    for (int row = 0; row < r; row++) {
        shuffle_vec(v, rng);
        for (int j = 0; j < n; j++) out(row, j) = v[j];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_objective_vals_cpp(
    const Rcpp::NumericMatrix& X,
    const Rcpp::IntegerMatrix& indicTs,
    std::string         objective,
    Rcpp::Nullable<Rcpp::NumericMatrix> inv_cov_X = R_NilValue
) {
    int n = X.nrow();
    int p = X.ncol();
    int r = indicTs.nrow();
    const bool abs_mode = (objective == "abs_sum_diff");
    const bool mahal_mode = (objective == "mahal_dist");

    if (indicTs.ncol() != n) {
        Rcpp::stop("indicTs must have ncol(indicTs) == nrow(X)");
    }
    if (!abs_mode && !mahal_mode) {
        Rcpp::stop("objective must be 'abs_sum_diff' or 'mahal_dist'");
    }

    std::vector<double> sum_all(p, 0.0), sumsq_all(p, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++) {
            double x = X(i, j);
            sum_all[j]   += x;
            sumsq_all[j] += x * x;
        }

    std::vector<double> sd_all;
    if (abs_mode) {
        sd_all.resize(p);
        for (int j = 0; j < p; j++) {
            double var = (sumsq_all[j] - sum_all[j] * sum_all[j] / n) / (n - 1);
            sd_all[j] = std::sqrt(std::max(var, 0.0));
        }
    }

    Rcpp::NumericMatrix Sinv;
    if (mahal_mode) {
        if (inv_cov_X.isNull()) Rcpp::stop("inv_cov_X required for mahal_dist");
        Sinv = Rcpp::NumericMatrix(inv_cov_X);
    }

    Rcpp::NumericVector vals(r);
    std::vector<double> sum_T(static_cast<std::size_t>(r) * p, 0.0);
    std::vector<int> nT_by_row(r, 0);
    std::vector<double> x_row(p), diff(p);
    const int* indic_ptr = indicTs.begin();
    const double* x_ptr = X.begin();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            x_row[j] = x_ptr[i + static_cast<std::size_t>(j) * n];
        }
        const int* indic_col = indic_ptr + static_cast<std::size_t>(i) * r;
        for (int row = 0; row < r; row++) {
            if (indic_col[row] == 1) {
                ++nT_by_row[row];
                double* row_sum_T = sum_T.data() + static_cast<std::size_t>(row) * p;
                for (int j = 0; j < p; j++) row_sum_T[j] += x_row[j];
            }
        }
    }

    for (int row = 0; row < r; row++) {
        const int nT = nT_by_row[row];
        int nC = n - nT;
        const double* row_sum_T = sum_T.data() + static_cast<std::size_t>(row) * p;
        for (int j = 0; j < p; j++)
            diff[j] = row_sum_T[j] / nT - (sum_all[j] - row_sum_T[j]) / nC;

        if (abs_mode) {
            double v = 0.0;
            for (int j = 0; j < p; j++) v += std::fabs(diff[j] / sd_all[j]);
            vals[row] = v;
        } else {
            double v = 0.0;
            for (int i = 0; i < p; i++)
                for (int k = 0; k < p; k++) v += diff[i] * Sinv(i, k) * diff[k];
            vals[row] = v;
        }
    }
    return vals;
}

// ── Parallel rejection-sampling rerandomization search ──────────────────────
//
// Generates balanced random allocations in parallel (OpenMP), evaluates the
// objective using the same M-matrix / Eigen SIMD path as the greedy search,
// and collects those whose objective value <= cutoff until r are found or
// max_draws is exhausted.
//
// Returns an n×k IntegerMatrix (k <= r).  The caller handles recycling when k < r.

// [[Rcpp::export]]
Rcpp::IntegerMatrix rerandomization_search_cpp(
    const Eigen::Map<Eigen::MatrixXd> X_raw,
    const int                         r,
    const std::string&                objective,
    const double                      cutoff,
    const int                         max_draws
) {
    const int n  = static_cast<int>(X_raw.rows());
    const int p  = static_cast<int>(X_raw.cols());
    const int nt = n / 2;

    // ── Build M (p×n): f = ‖M*(2w−1)‖₁  or  ‖M*(2w−1)‖² ───────────────────
    RowVectorXd col_means = X_raw.colwise().mean();
    MatrixXd    X         = X_raw.rowwise() - col_means;

    MatrixXd M(p, n);
    bool abs_mode = true;

    if (objective == "mahal_dist") {
        MatrixXd cov_mat = (X.transpose() * X) / std::max(1, n - 1);
        LLT<MatrixXd> llt(cov_mat);
        if (llt.info() == Success) {
            M        = llt.matrixL().solve(X.transpose()) / static_cast<double>(n);
            abs_mode = false;
        }
        // singular covariance → fall through to abs_sum_diff
    }

    if (abs_mode) {
        VectorXd inv_sd(p);
        for (int j = 0; j < p; j++) {
            double var = X.col(j).squaredNorm() / std::max(1, n - 1);
            inv_sd[j]  = (var < 1e-24) ? 1.0 : 1.0 / std::sqrt(var);
        }
        M = (X * inv_sd.asDiagonal()).transpose() / static_cast<double>(n);
    }

    // ── Seed per-thread RNGs from R's RNG ────────────────────────────────────
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<uint32_t> seeds(static_cast<std::size_t>(nthreads));
    GetRNGstate();
    for (int t = 0; t < nthreads; t++)
        seeds[static_cast<std::size_t>(t)] =
            static_cast<uint32_t>(::unif_rand() * 4294967295.0);
    PutRNGstate();

    // ── Shared result storage ────────────────────────────────────────────────
    // Pre-allocate n×r; only the first `found` columns are valid on return.
    Rcpp::IntegerMatrix result(n, r);
    const bool use_sequential = (nthreads <= 1);

    if (use_sequential) {
        std::mt19937 rng(seeds[0]);
        std::vector<int> w(n), order(n);
        VectorXd dbl_w(n), d(p);
        int found = 0;

        for (int draw = 0; draw < max_draws && found < r; draw++) {
            std::iota(order.begin(), order.end(), 0);
            for (int i = n - 1; i > 0; i--)
                std::swap(order[static_cast<std::size_t>(i)],
                          order[static_cast<std::size_t>(
                              std::uniform_int_distribution<int>(0, i)(rng))]);
            std::fill(w.begin(), w.end(), 0);
            for (int i = 0; i < nt; i++) w[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = 1;

            for (int i = 0; i < n; i++) dbl_w[i] = 2.0 * static_cast<double>(w[static_cast<std::size_t>(i)]) - 1.0;
            d.noalias() = M * dbl_w;
            double f = abs_mode ? d.lpNorm<1>() : d.squaredNorm();

            if (f <= cutoff) {
                for (int i = 0; i < n; i++) result(i, found) = w[static_cast<std::size_t>(i)];
                ++found;
            }
        }

        if (found == r) return result;
        Rcpp::IntegerMatrix trimmed(n, found);
        for (int j = 0; j < found; j++)
            for (int i = 0; i < n; i++) trimmed(i, j) = result(i, j);
        return trimmed;
    }

    std::atomic<int> found(0);
    std::atomic<int> next_draw(0);
    const int chunk_size = 64;

    // ── Parallel rejection loop ──────────────────────────────────────────────
#pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::mt19937 rng(seeds[static_cast<std::size_t>(tid)]);

        // thread-local buffers
        std::vector<int> w(n), order(n);
        VectorXd dbl_w(n), d(p);

        while (found.load(std::memory_order_relaxed) < r) {
            const int start = next_draw.fetch_add(chunk_size, std::memory_order_relaxed);
            if (start >= max_draws) break;
            const int end = std::min(max_draws, start + chunk_size);

            for (int draw = start; draw < end && found.load(std::memory_order_relaxed) < r; draw++) {
                // Fisher-Yates balanced init
                std::iota(order.begin(), order.end(), 0);
                for (int i = n - 1; i > 0; i--)
                    std::swap(order[static_cast<std::size_t>(i)],
                              order[static_cast<std::size_t>(
                                  std::uniform_int_distribution<int>(0, i)(rng))]);
                std::fill(w.begin(), w.end(), 0);
                for (int i = 0; i < nt; i++) w[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = 1;

                for (int i = 0; i < n; i++) dbl_w[i] = 2.0 * static_cast<double>(w[static_cast<std::size_t>(i)]) - 1.0;
                d.noalias() = M * dbl_w;
                double f = abs_mode ? d.lpNorm<1>() : d.squaredNorm();

                if (f <= cutoff) {
                    int slot = found.fetch_add(1, std::memory_order_relaxed);
                    if (slot < r) {
                        for (int i = 0; i < n; i++) result(i, slot) = w[static_cast<std::size_t>(i)];
                    }
                }
            }
        }
    }

    int k = std::min(found.load(), r);
    if (k == r) return result;
    // Return only the filled columns
    Rcpp::IntegerMatrix trimmed(n, k);
    for (int j = 0; j < k; j++)
        for (int i = 0; i < n; i++) trimmed(i, j) = result(i, j);
    return trimmed;
}
