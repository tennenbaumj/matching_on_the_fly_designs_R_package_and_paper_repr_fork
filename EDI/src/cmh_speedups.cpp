#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
double compute_cmh_block_se_cpp(const NumericVector& y,
                                   const IntegerVector& m_vec,
                                   int n_total) {
  if (y.size() != m_vec.size()) {
    stop("compute_cmh_block_se_cpp: y and m_vec must have the same length.");
  }
  if (n_total <= 0) {
    return NA_REAL;
  }

  // Find max block ID (skipping NA/invalid entries) to size the flat accumulator.
  int max_block_id = 0;
  for (int i = 0; i < (int)m_vec.size(); ++i) {
    const int mid = m_vec[i];
    if (mid != NA_INTEGER && mid > max_block_id) max_block_id = mid;
  }
  if (max_block_id <= 0) return NA_REAL;

  // y is binary (0/1), so y[i]^2 == y[i] and sum_sq == sum; store integer sum only.
  // -1 sentinel: block not yet seen (distinguishes "unseen" from "seen but all zeros").
  std::vector<int> block_sums(max_block_id + 1, -1);
  int B = 0;
  int n_included = 0;
  for (int i = 0; i < (int)y.size(); ++i) {
    const int match_id = m_vec[i];
    if (match_id == NA_INTEGER || match_id <= 0 || !R_finite(y[i])) {
      continue;
    }
    int& bs = block_sums[match_id];
    if (bs < 0) { bs = 0; ++B; }
    bs += static_cast<int>(y[i]);
    ++n_included;
  }

  if (B <= 0) {
    return NA_REAL;
  }

  // Equal block sizes guaranteed by assert_equal_block_sizes on the R side.
  // Derive block size once and cast to double once.
  const int    n_b   = n_included / B;
  if (n_b <= 1) {
    return NA_REAL;
  }
  const double n_b_d = static_cast<double>(n_b);
  const double denom = n_b_d - 1.0;

  double var_cmh = 0.0;
  for (int bid = 1; bid <= max_block_id; ++bid) {
    if (block_sums[bid] < 0) continue;  // unseen block
    const double s = static_cast<double>(block_sums[bid]);
    // ss  = sum_sq - sum^2/n_b = s - s^2/n_b = s*(n_b - s)/n_b
    // contribution = (n_b/(n_b-1)) * ss = s*(n_b - s)/(n_b - 1)
    var_cmh += s * (n_b_d - s) / denom;
  }

  return (2.0 / static_cast<double>(n_total)) * std::sqrt(var_cmh);
}

// [[Rcpp::export]]
double compute_extended_robins_block_se_cpp(const NumericVector& y,
                                            const NumericVector& w,
                                            const IntegerVector& m_vec,
                                            int n_total) {
  if (y.size() != m_vec.size() || y.size() != w.size()) {
    stop("compute_extended_robins_block_se_cpp: y, w, and m_vec must have the same length.");
  }
  if (n_total <= 0) {
    return NA_REAL;
  }

  // Find max block ID to size the flat accumulator.
  int max_block_id = 0;
  for (int i = 0; i < (int)m_vec.size(); ++i) {
    const int mid = m_vec[i];
    if (mid != NA_INTEGER && mid > max_block_id) max_block_id = mid;
  }
  if (max_block_id <= 0) return NA_REAL;

  struct RobinsBlockAccumulator {
    int n = 0;
    double sum_t = 0.0;
    double sum_c = 0.0;
  };

  std::vector<RobinsBlockAccumulator> flat_blocks(max_block_id + 1);
  int B = 0;

  int total_t = 0;
  int total_c = 0;
  double total_sum_t = 0.0;
  double total_sum_c = 0.0;

  for (int i = 0; i < (int)y.size(); ++i) {
    const int match_id = m_vec[i];
    if (match_id == NA_INTEGER || match_id <= 0) {
      continue;
    }
    if (!R_finite(y[i]) || !R_finite(w[i]) || (w[i] != -1.0 && w[i] != 1.0)) {
      return NA_REAL;
    }

    RobinsBlockAccumulator& block = flat_blocks[match_id];
    if (block.n == 0) ++B;
    block.n += 1;
    if (w[i] == 1.0) {
      ++total_t;
      block.sum_t += y[i];
      total_sum_t += y[i];
    } else {
      ++total_c;
      block.sum_c += y[i];
      total_sum_c += y[i];
    }
  }

  if (B <= 0) {
    return NA_REAL;
  }

  if (total_t <= 0 || total_c <= 0) {
    return NA_REAL;
  }

  double variance_tot = 0.0;
  for (int bid = 1; bid <= max_block_id; ++bid) {
    const RobinsBlockAccumulator& block = flat_blocks[bid];
    if (block.n == 0) continue;  // unseen block
    if (block.n <= 1) {
      return NA_REAL;
    }

    const double n_b = static_cast<double>(block.n);
    const double n_b_over_two = n_b / 2.0;
    const double p_hat_T_b = block.sum_t / n_b_over_two;
    const double p_hat_C_b = block.sum_c / n_b_over_two;
    const double m_1_b = std::max(p_hat_T_b, p_hat_C_b);
    const double m_0_b = std::min(p_hat_T_b, p_hat_C_b);

    variance_tot +=
      m_1_b * (1.0 - m_1_b) / n_b_over_two +
      m_0_b * (1.0 - m_0_b) / n_b_over_two +
      ((2.0 * m_0_b - m_1_b) * (1.0 - m_1_b) - m_0_b * (1.0 - m_0_b)) / n_b;
  }
  //now divide by number of blocks squared
  variance_tot /= static_cast<double>(B * B);

  const double p_hat_T = total_sum_t / static_cast<double>(total_t);
  const double p_hat_C = total_sum_c / static_cast<double>(total_c);
  const double var_robbins_ext =
    (p_hat_T * (1.0 - p_hat_T) + p_hat_C * (1.0 - p_hat_C)) /
    static_cast<double>(n_total);

  return std::sqrt(variance_tot + var_robbins_ext);
}
