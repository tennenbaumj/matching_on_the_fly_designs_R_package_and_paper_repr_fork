#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <cstdint>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

// Helper: which_cols_vary_subset
std::vector<int> which_cols_vary_subset_int(const MatrixXd& X, int n_rows) {
	int p = X.cols();
	std::vector<int> var_cols;
	for (int j = 0; j < p; ++j) {
		if (n_rows <= 1) continue;
		double first_val = X(0, j);
		bool varies = false;
		for (int i = 1; i < n_rows; ++i) {
			if (X(i, j) != first_val) {
				varies = true;
				break;
			}
		}
		if (varies) var_cols.push_back(j);
	}
	return var_cols;
}

// Helper: extract_submatrix
MatrixXd extract_submatrix_int(const MatrixXd& X, int n_rows, const std::vector<int>& cols) {
	MatrixXd res(n_rows, cols.size());
	for (int i = 0; i < n_rows; ++i) {
		for (size_t j = 0; j < cols.size(); ++j) {
			res(i, j) = X(i, cols[j]);
		}
	}
	return res;
}

// Helper: find_independent_cols
std::vector<int> find_independent_cols_int(const MatrixXd& X) {
	int p = X.cols();
	if (p == 0) return std::vector<int>();
	FullPivHouseholderQR<MatrixXd> qr(X);
	int rank = qr.rank();
	std::vector<int> indep_cols;
	if (rank == 0) return indep_cols;
	MatrixXi P = qr.colsPermutation().indices();
	for (int i = 0; i < rank; ++i) indep_cols.push_back(P(i));
	std::sort(indep_cols.begin(), indep_cols.end());
	return indep_cols;
}

// Helper: compute_atkinson_weight_internal
double compute_atkinson_weight_internal(const VectorXd& w_prev, const MatrixXd& X_prev, const VectorXd& xt_prev, int t) {
	int rows = w_prev.size();
	int p = X_prev.cols();
	int cols = p + 2;
	if (rows == 0 || cols < 2) return 0.5;

	MatrixXd XprevWT(rows, cols);
	XprevWT.col(0) = w_prev;
	XprevWT.col(1).setOnes();
	XprevWT.rightCols(p) = X_prev;

	MatrixXd XwtXw = XprevWT.transpose() * XprevWT;
	FullPivLU<MatrixXd> lu(XwtXw);
	if (!lu.isInvertible()) return 0.5;

	MatrixXd M = static_cast<double>(t - 1) * lu.inverse();
	VectorXd row_segment = M.row(0).segment(1, p + 1);
	VectorXd xt(p + 1);
	xt(0) = 1.0; xt.tail(p) = xt_prev;

	double A = row_segment.dot(xt);
	if (A == 0 || !std::isfinite(A)) return 0.5;

	double val = M(0, 0) / A + 1.0;
	double s_over_A_plus_one_sq = val * val;
	double prob = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1.0);
	return std::max(0.0, std::min(1.0, prob));
}

// Local RNG helpers — seed from R once per call, then run independently
inline double u01(std::mt19937_64& rng) {
	return (rng() >> 11) * 0x1.0p-53;
}

inline std::mt19937_64 make_local_rng() {
	return std::mt19937_64(static_cast<uint64_t>(
		R::unif_rand() * static_cast<double>(std::numeric_limits<uint64_t>::max())));
}

} // namespace

// [[Rcpp::export]]
List generate_permutations_matching_cpp(const IntegerVector& m_vec, int nsim, double prob_T) {
  auto rng = make_local_rng();
  int n = m_vec.size();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  const int* m_orig_ptr = m_vec.begin();

  int max_m = 0;
  for (int i = 0; i < n; ++i) if (m_orig_ptr[i] > max_m) max_m = m_orig_ptr[i];

  std::vector<std::vector<int>> pairs(max_m);
  std::vector<int> reservoir;
  for (int i = 0; i < n; ++i) {
    int m = m_orig_ptr[i];
    if (m > 0) pairs[m-1].push_back(i);
    else reservoir.push_back(i);
  }

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int i : reservoir) w_col[i] = (u01(rng) < prob_T) ? 1 : 0;
    for (int m = 0; m < max_m; ++m) {
      if (pairs[m].size() == 2) {
        int first_is_T = (u01(rng) < prob_T) ? 1 : 0;
        w_col[pairs[m][0]] = first_is_T;
        w_col[pairs[m][1]] = 1 - first_is_T;
      } else if (pairs[m].size() == 1) {
        w_col[pairs[m][0]] = (u01(rng) < prob_T) ? 1 : 0;
      }
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue); // Assume m_vec is fixed for KK randomization test
}

// [[Rcpp::export]]
List generate_permutations_bernoulli_cpp(int n, int nsim, double prob_T) {
  auto rng = make_local_rng();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int i = 0; i < n; ++i) w_col[i] = (u01(rng) < prob_T) ? 1 : 0;
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_ibcrd_cpp(int n, int nsim, double prob_T) {
  auto rng = make_local_rng();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  int n_T = std::round(n * prob_T);
  std::vector<int> w_base(n);
  for (int i = 0; i < n; ++i) w_base[i] = (i < n_T) ? 1 : 0;
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    std::vector<int> w_shuffled = w_base;
    std::shuffle(w_shuffled.begin(), w_shuffled.end(), rng);
    for (int i = 0; i < n; ++i) w_col[i] = w_shuffled[i];
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_blocking_cpp(int n, int nsim, double prob_T, List strata_indices) {
  auto rng = make_local_rng();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  int num_strata = strata_indices.size();
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int s = 0; s < num_strata; ++s) {
      IntegerVector idxs = strata_indices[s];
      int m = idxs.size();
      int n_T = std::round(m * prob_T);
      std::vector<int> w_stratum(m);
      for (int i = 0; i < m; ++i) w_stratum[i] = (i < n_T) ? 1 : 0;
      std::shuffle(w_stratum.begin(), w_stratum.end(), rng);
      for (int i = 0; i < m; ++i) w_col[idxs[i] - 1] = w_stratum[i];
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_efron_cpp(int n, int nsim, double prob_T, double weighted_coin_prob) {
  auto rng = make_local_rng();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    int n_T = 0, n_C = 0;
    for (int i = 0; i < n; ++i) {
      double p;
      double sT = n_T * prob_T, sC = n_C * (1.0 - prob_T);
      if (sT > sC) p = 1.0 - weighted_coin_prob;
      else if (sT < sC) p = weighted_coin_prob;
      else p = prob_T;
      if (u01(rng) < p) { w_col[i] = 1; n_T++; }
      else { w_col[i] = 0; n_C++; }
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_atkinson_cpp(SEXP X_sexp, int n, int p_raw, double prob_T, int nsim) {
  auto rng = make_local_rng();
  Rcpp::NumericMatrix X_r(X_sexp);
  Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());

  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  int bernoulli_threshold = p_raw + 2 + 1;
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int t = 1; t <= n; ++t) {
      if (t <= bernoulli_threshold) {
        w_col[t - 1] = (u01(rng) < prob_T) ? 1 : 0;
        continue;
      }
      std::vector<int> var_cols = which_cols_vary_subset_int(X, t);
      if (var_cols.empty()) { w_col[t - 1] = (u01(rng) < prob_T) ? 1 : 0; continue; }
      MatrixXd X_var = extract_submatrix_int(X, t - 1, var_cols);
      std::vector<int> indep_cols = find_independent_cols_int(X_var);
      if (indep_cols.empty()) { w_col[t - 1] = (u01(rng) < prob_T) ? 1 : 0; continue; }
      MatrixXd X_prev(t - 1, indep_cols.size());
      for (int i = 0; i < t - 1; ++i) {
        for (size_t j = 0; j < indep_cols.size(); ++j) X_prev(i, j) = X_var(i, indep_cols[j]);
      }
      VectorXd xt_prev(indep_cols.size());
      for (size_t j = 0; j < indep_cols.size(); ++j) xt_prev(j) = X(t - 1, var_cols[indep_cols[j]]);
      VectorXd w_prev_e(t - 1);
      for (int i = 0; i < t - 1; ++i) w_prev_e(i) = static_cast<double>(w_col[i]);
      double p = compute_atkinson_weight_internal(w_prev_e, X_prev, xt_prev, t);
      w_col[t - 1] = (u01(rng) < p) ? 1 : 0;
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_pocock_simon_cpp(const IntegerMatrix& x_levels_matrix, int num_levels_total, const NumericVector& weights, double p_best, double prob_T, int nsim) {
  auto rng = make_local_rng();
  int n = x_levels_matrix.nrow();
  int num_covs = x_levels_matrix.ncol();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    NumericMatrix counts(num_levels_total, 2); // 2 treatments

    for (int i = 0; i < n; ++i) {
      IntegerVector subject_levels = x_levels_matrix.row(i);

      // Inline the assign logic to avoid R API contention from IntegerVector
      std::vector<double> G(2, 0.0);
      for (int k = 0; k < 2; ++k) {
        double G_k = 0.0;
        for (int j = 0; j < num_covs; ++j) {
          int row_idx = subject_levels[j] - 1;
          double c0 = counts(row_idx, 0) + (k == 0 ? 1 : 0);
          double c1 = counts(row_idx, 1) + (k == 1 ? 1 : 0);
          double m = (c0 + c1) / 2.0;
          double var = (c0-m)*(c0-m) + (c1-m)*(c1-m); // Variance proportional
          G_k += weights[j] * var;
        }
        G[k] = G_k;
      }

      int best_trt = (G[1] < G[0]) ? 1 : (G[0] < G[1] ? 0 : ((u01(rng) < prob_T) ? 1 : 0));
      int assigned_w = (u01(rng) < p_best) ? best_trt : (1 - best_trt);

      w_col[i] = assigned_w;
      for (int j = 0; j < num_covs; ++j) {
        counts(subject_levels[j] - 1, assigned_w)++;
      }
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_cluster_cpp(int n, int nsim, double prob_T, List cluster_indices) {
  auto rng = make_local_rng();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  int num_clusters = cluster_indices.size();
  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int c = 0; c < num_clusters; ++c) {
      IntegerVector idxs = cluster_indices[c];
      int w_cluster = (u01(rng) < prob_T) ? 1 : 0;
      int m = idxs.size();
      for (int i = 0; i < m; ++i) w_col[idxs[i] - 1] = w_cluster;
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}

// [[Rcpp::export]]
List generate_permutations_spbr_cpp(const CharacterVector& strata_keys, int block_size, double prob_T, int nsim) {
  auto rng = make_local_rng();
  int n = strata_keys.size();
  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();

  int n_T_block = std::round(block_size * prob_T);

  // Pre-convert string keys to dense integer IDs once — eliminates per-sim string map lookups
  std::unordered_map<std::string, int> key_to_id;
  key_to_id.reserve(64);
  std::vector<int> strata_ids(n);
  int num_strata = 0;
  for (int i = 0; i < n; ++i) {
    std::string key = as<std::string>(strata_keys[i]);
    auto result = key_to_id.emplace(key, num_strata);
    if (result.second) ++num_strata;
    strata_ids[i] = result.first->second;
  }

  // Base block shuffled per new block; reuse allocation across simulations
  std::vector<int> base_block(block_size);
  for (int k = 0; k < block_size; ++k) base_block[k] = (k < n_T_block) ? 1 : 0;

  // Persistent strata state — clear() at top of each simulation keeps capacity
  std::vector<std::vector<int>> strata_states(num_strata);

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (auto& v : strata_states) v.clear();

    for (int i = 0; i < n; ++i) {
      int sid = strata_ids[i];
      if (strata_states[sid].empty()) {
        strata_states[sid] = base_block;
        std::shuffle(strata_states[sid].begin(), strata_states[sid].end(), rng);
      }
      w_col[i] = strata_states[sid].back();
      strata_states[sid].pop_back();
    }
  }
  return List::create(_["w_mat"] = w_mat, _["m_mat"] = R_NilValue);
}
