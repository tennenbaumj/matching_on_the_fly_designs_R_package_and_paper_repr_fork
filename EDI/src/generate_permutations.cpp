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
std::vector<int> find_independent_cols_int(
	const MatrixXd& X,
	FullPivHouseholderQR<MatrixXd>& qr
) {
	int p = X.cols();
	if (p == 0) return std::vector<int>();
	qr.compute(X);
	int rank = qr.rank();
	std::vector<int> indep_cols;
	if (rank == 0) return indep_cols;
	MatrixXi P = qr.colsPermutation().indices();
	for (int i = 0; i < rank; ++i) indep_cols.push_back(P(i));
	std::sort(indep_cols.begin(), indep_cols.end());
	return indep_cols;
}

struct AtkinsonStepData {
	bool usable = false;
	MatrixXd X_prev;
	VectorXd xt_prev;
};

// Helper: compute_atkinson_weight_internal
double compute_atkinson_weight_internal(
	const int* w_prev,
	const MatrixXd& X_prev,
	const VectorXd& xt_prev,
	int t,
	FullPivLU<MatrixXd>& lu,
	MatrixXd& design_workspace,
	MatrixXd& crossprod_workspace,
	MatrixXd& inverse_workspace,
	VectorXd& xt_workspace
) {
	int rows = t - 1;
	int p = X_prev.cols();
	int cols = p + 2;
	if (rows == 0 || cols < 2) return 0.5;

	auto XprevWT = design_workspace.topLeftCorner(rows, cols);
	for (int i = 0; i < rows; ++i) XprevWT(i, 0) = static_cast<double>(w_prev[i]);
	XprevWT.col(1).setOnes();
	XprevWT.rightCols(p) = X_prev;

	auto XwtXw = crossprod_workspace.topLeftCorner(cols, cols);
	XwtXw.noalias() = XprevWT.transpose() * XprevWT;
	lu.compute(XwtXw);
	if (!lu.isInvertible()) return 0.5;

	auto M = inverse_workspace.topLeftCorner(cols, cols);
	M.noalias() = static_cast<double>(t - 1) * lu.inverse();
	auto xt = xt_workspace.head(p + 1);
	xt(0) = 1.0;
	xt.tail(p) = xt_prev;

	double A = M.row(0).segment(1, p + 1).dot(xt);
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

	// The varying/independent covariate columns depend only on X and t, not on
	// the simulated assignments. Compute each QR and processed design once.
	std::vector<AtkinsonStepData> step_data(static_cast<std::size_t>(n));
	FullPivHouseholderQR<MatrixXd> qr_workspace;
	for (int t = bernoulli_threshold + 1; t <= n; ++t) {
		std::vector<int> var_cols = which_cols_vary_subset_int(X, t);
		if (var_cols.empty()) continue;
		MatrixXd X_var = extract_submatrix_int(X, t - 1, var_cols);
		std::vector<int> indep_cols = find_independent_cols_int(X_var, qr_workspace);
		if (indep_cols.empty()) continue;

		AtkinsonStepData& step = step_data[static_cast<std::size_t>(t - 1)];
		step.X_prev.resize(t - 1, static_cast<int>(indep_cols.size()));
		step.xt_prev.resize(static_cast<int>(indep_cols.size()));
		for (int i = 0; i < t - 1; ++i) {
			for (std::size_t j = 0; j < indep_cols.size(); ++j) {
				step.X_prev(i, static_cast<int>(j)) = X_var(i, indep_cols[j]);
			}
		}
		for (std::size_t j = 0; j < indep_cols.size(); ++j) {
			step.xt_prev(static_cast<int>(j)) = X(t - 1, var_cols[indep_cols[j]]);
		}
		step.usable = true;
	}

	const int max_cols = X.cols() + 2;
	FullPivLU<MatrixXd> lu_workspace;
	MatrixXd design_workspace(n, max_cols);
	MatrixXd crossprod_workspace(max_cols, max_cols);
	MatrixXd inverse_workspace(max_cols, max_cols);
	VectorXd xt_workspace(max_cols - 1);

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int t = 1; t <= n; ++t) {
      if (t <= bernoulli_threshold) {
        w_col[t - 1] = (u01(rng) < prob_T) ? 1 : 0;
        continue;
      }
      const AtkinsonStepData& step = step_data[static_cast<std::size_t>(t - 1)];
      if (!step.usable) { w_col[t - 1] = (u01(rng) < prob_T) ? 1 : 0; continue; }
      double p = compute_atkinson_weight_internal(
		w_col, step.X_prev, step.xt_prev, t, lu_workspace,
		design_workspace, crossprod_workspace, inverse_workspace, xt_workspace);
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
  if (num_levels_total <= 0) stop("num_levels_total must be positive");
  if (weights.size() != num_covs) stop("weights length must match the number of covariates");

  // Convert the R matrix's one-based global level rows to a cache-friendly,
  // zero-based row-major buffer once, validating before any pointer indexing.
  std::vector<int> level_rows(static_cast<std::size_t>(n) * num_covs);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < num_covs; ++j) {
      const int row_idx = x_levels_matrix(i, j) - 1;
      if (row_idx < 0 || row_idx >= num_levels_total) {
        stop("x_levels_matrix contains a level index outside 1..num_levels_total");
      }
      level_rows[static_cast<std::size_t>(i) * num_covs + j] = row_idx;
    }
  }

  IntegerMatrix w_mat(n, nsim);
  int* w_ptr = w_mat.begin();
  const double* weights_ptr = weights.begin();
  std::vector<int> counts(static_cast<std::size_t>(num_levels_total) * 2);

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    std::fill(counts.begin(), counts.end(), 0);

    for (int i = 0; i < n; ++i) {
      const int* subject_levels = level_rows.data() + static_cast<std::size_t>(i) * num_covs;

      double G[2] = {0.0, 0.0};
      for (int k = 0; k < 2; ++k) {
        double G_k = 0.0;
        for (int j = 0; j < num_covs; ++j) {
          const int row_idx = subject_levels[j];
          double c0 = counts[static_cast<std::size_t>(row_idx) * 2] + (k == 0 ? 1 : 0);
          double c1 = counts[static_cast<std::size_t>(row_idx) * 2 + 1] + (k == 1 ? 1 : 0);
          double m = (c0 + c1) / 2.0;
          double var = (c0-m)*(c0-m) + (c1-m)*(c1-m); // Variance proportional
          G_k += weights_ptr[j] * var;
        }
        G[k] = G_k;
      }

      int best_trt = (G[1] < G[0]) ? 1 : (G[0] < G[1] ? 0 : ((u01(rng) < prob_T) ? 1 : 0));
      int assigned_w = (u01(rng) < p_best) ? best_trt : (1 - best_trt);

      w_col[i] = assigned_w;
      for (int j = 0; j < num_covs; ++j) {
        counts[static_cast<std::size_t>(subject_levels[j]) * 2 + assigned_w]++;
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

  // Flatten the R list once. Offsets retain cluster order and allow overlapping
  // clusters to preserve the original last-write-wins behavior.
  std::vector<std::size_t> cluster_offsets(static_cast<std::size_t>(num_clusters) + 1);
  std::vector<int> subject_indices;
  subject_indices.reserve(static_cast<std::size_t>(n));
  for (int c = 0; c < num_clusters; ++c) {
    cluster_offsets[static_cast<std::size_t>(c)] = subject_indices.size();
    IntegerVector idxs = cluster_indices[c];
    for (int i = 0; i < idxs.size(); ++i) {
      const int subject = idxs[i] - 1;
      if (subject < 0 || subject >= n) stop("cluster_indices contains an index outside 1..n");
      subject_indices.push_back(subject);
    }
  }
  cluster_offsets[static_cast<std::size_t>(num_clusters)] = subject_indices.size();

  for (int b = 0; b < nsim; ++b) {
    int* w_col = w_ptr + (size_t)b * n;
    for (int c = 0; c < num_clusters; ++c) {
      int w_cluster = (u01(rng) < prob_T) ? 1 : 0;
      const std::size_t begin = cluster_offsets[static_cast<std::size_t>(c)];
      const std::size_t end = cluster_offsets[static_cast<std::size_t>(c + 1)];
      for (std::size_t i = begin; i < end; ++i) w_col[subject_indices[i]] = w_cluster;
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
