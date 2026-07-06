#include <RcppEigen.h>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <cstdint>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

namespace {

inline bool kk_is_success(double y_i) {
	return R_finite(y_i) && y_i != 0.0;
}

inline uint64_t bounded_rand(std::mt19937_64& rng, uint64_t s) {
	uint64_t x = rng();
	__uint128_t m = static_cast<__uint128_t>(x) * s;
	uint64_t l = static_cast<uint64_t>(m);
	if (l < s) {
		uint64_t t = (-s) % s;
		while (l < t) {
			x = rng();
			m = static_cast<__uint128_t>(x) * s;
			l = static_cast<uint64_t>(m);
		}
	}
	return static_cast<uint64_t>(m >> 64);
}

inline std::mt19937_64 make_local_rng() {
	return std::mt19937_64(static_cast<uint64_t>(
		R::unif_rand() * static_cast<double>(std::numeric_limits<uint64_t>::max())));
}

} // namespace

// [[Rcpp::export]]
IntegerMatrix bootstrap_m_indices_cpp(
	const IntegerVector& m_vec,
	const IntegerVector& i_reservoir,
	int n_reservoir,
	int m,
	int B
) {
	auto rng = make_local_rng();
	int row_length = n_reservoir + 2 * m;
	IntegerMatrix result(B, row_length);

	std::vector< std::array<int, 2> > match_pairs(m);
	std::vector<int> count(m, 0);
	for (int idx = 0; idx < m_vec.size(); ++idx) {
	int match_id = m_vec[idx];
	if (match_id > 0 && match_id <= m) {
		int pos = count[match_id - 1]++;
		if (pos < 2) {
		match_pairs[match_id - 1][pos] = idx + 1;
		}
	}
	}

	const uint64_t u_res = static_cast<uint64_t>(n_reservoir);
	const uint64_t u_m   = static_cast<uint64_t>(m);
	for (int row = 0; row < B; ++row) {
	if (n_reservoir > 0) {
		for (int j = 0; j < n_reservoir; ++j) {
		result(row, j) = i_reservoir[static_cast<int>(bounded_rand(rng, u_res))];
		}
	}

	for (int k = 0; k < m; ++k) {
		int match_id = static_cast<int>(bounded_rand(rng, u_m));
		auto pair = match_pairs[match_id];
		result(row, n_reservoir + 2 * k) = pair[0];
		result(row, n_reservoir + 2 * k + 1) = pair[1];
	}
	}
	return result;
}

List compute_zhang_match_data_cpp(const NumericMatrix& X,
								  const NumericVector& y,
								  const IntegerVector& w,
								  const IntegerVector& m_vec);

// [[Rcpp::export]]
List draw_matching_bootstrap_sample_cpp(
	const IntegerVector& i_reservoir,
	const IntegerMatrix& pair_rows,
	int n_reservoir
) {
	const int m = pair_rows.nrow();
	const int out_n = n_reservoir + 2 * m;
	IntegerVector i_b(out_n);
	IntegerVector m_vec_b(out_n);
	Function sample_fn("sample");

	for (int j = 0; j < n_reservoir; ++j) {
		m_vec_b[j] = 0;
	}
	if (n_reservoir > 0) {
		IntegerVector i_reservoir_b = sample_fn(i_reservoir, Named("size", n_reservoir), Named("replace", true));
		for (int j = 0; j < n_reservoir; ++j) {
			i_b[j] = i_reservoir_b[j];
		}
	}

	if (m > 0) {
		IntegerVector pair_ids = seq_len(m);
		IntegerVector pairs_to_include = sample_fn(pair_ids, Named("size", m), Named("replace", true));
		for (int pair_idx = 0; pair_idx < m; ++pair_idx) {
			const int sampled_pair = pairs_to_include[pair_idx] - 1;
			const int out_idx = n_reservoir + 2 * pair_idx;
			i_b[out_idx] = pair_rows(sampled_pair, 0);
			i_b[out_idx + 1] = pair_rows(sampled_pair, 1);
			m_vec_b[out_idx] = pair_idx + 1;
			m_vec_b[out_idx + 1] = pair_idx + 1;
		}
	}

	return List::create(
		_["i_b"] = i_b,
		_["m_vec_b"] = m_vec_b
	);
}

// [[Rcpp::export]]
List compute_bootstrap_matching_stats_cpp(
	const NumericMatrix& X,
	const NumericVector& y,
	const IntegerVector& w,
	const IntegerVector& i_b,
	int n_reservoir
) {
	const int n_rows = i_b.size();
	const int p = X.ncol();
	const int matched_rows = std::max(0, n_rows - n_reservoir);
	const int m = matched_rows / 2;

	NumericVector yTs_matched(m, NA_REAL);
	NumericVector yCs_matched(m, NA_REAL);
	NumericVector y_matched_diffs(m, NA_REAL);
	NumericMatrix X_matched_diffs_full(m, p);
	NumericMatrix X_matched_means_full(m, p);
	std::vector<int> found_t(static_cast<std::size_t>(m), 0);
	std::vector<int> found_c(static_cast<std::size_t>(m), 0);

	NumericMatrix X_reservoir(n_reservoir, p);
	NumericVector y_reservoir(n_reservoir);
	IntegerVector w_reservoir(n_reservoir);
	int nRT = 0;
	int nRC = 0;
	int n11 = 0;
	int n10 = 0;
	int n01 = 0;
	int n00 = 0;

	for (int res_idx = 0; res_idx < n_reservoir; ++res_idx) {
		const int src_idx = i_b[res_idx] - 1;
		const int w_i = w[src_idx];
		const double y_i = y[src_idx];
		const bool success = kk_is_success(y_i);

		y_reservoir[res_idx] = y_i;
		w_reservoir[res_idx] = w_i;
		for (int j = 0; j < p; ++j) {
			X_reservoir(res_idx, j) = X(src_idx, j);
		}

		if (w_i == 1) {
			++nRT;
			if (success) {
				++n11;
			} else {
				++n10;
			}
		} else {
			++nRC;
			if (success) {
				++n01;
			} else {
				++n00;
			}
		}
	}

	for (int pair_idx = 0; pair_idx < m; ++pair_idx) {
		const int row1 = i_b[n_reservoir + 2 * pair_idx] - 1;
		const int row2 = i_b[n_reservoir + 2 * pair_idx + 1] - 1;
		const int rows[2] = {row1, row2};

		for (int rr = 0; rr < 2; ++rr) {
			const int src_idx = rows[rr];
			const int w_i = w[src_idx];
			const double y_i = y[src_idx];

			if (w_i == 1) {
				yTs_matched[pair_idx] = y_i;
				found_t[static_cast<std::size_t>(pair_idx)] = 1;
				for (int j = 0; j < p; ++j) {
					X_matched_diffs_full(pair_idx, j) += X(src_idx, j);
					X_matched_means_full(pair_idx, j) += X(src_idx, j) / 2.0;
				}
			} else {
				yCs_matched[pair_idx] = y_i;
				found_c[static_cast<std::size_t>(pair_idx)] = 1;
				for (int j = 0; j < p; ++j) {
					X_matched_diffs_full(pair_idx, j) -= X(src_idx, j);
					X_matched_means_full(pair_idx, j) += X(src_idx, j) / 2.0;
				}
			}
		}
	}

	int d_plus = 0;
	int d_minus = 0;
	std::vector<int> keep_cols;
	keep_cols.reserve(static_cast<std::size_t>(p));
	for (int pair_idx = 0; pair_idx < m; ++pair_idx) {
		if (found_t[static_cast<std::size_t>(pair_idx)] &&
			found_c[static_cast<std::size_t>(pair_idx)]) {
			y_matched_diffs[pair_idx] = yTs_matched[pair_idx] - yCs_matched[pair_idx];
			const bool y_t_success = kk_is_success(yTs_matched[pair_idx]);
			const bool y_c_success = kk_is_success(yCs_matched[pair_idx]);
			if (y_t_success && !y_c_success) {
				++d_plus;
			} else if (!y_t_success && y_c_success) {
				++d_minus;
			}
		}
	}

	for (int j = 0; j < p; ++j) {
		bool nonzero = false;
		for (int pair_idx = 0; pair_idx < m; ++pair_idx) {
			if (X_matched_diffs_full(pair_idx, j) != 0.0) {
				nonzero = true;
				break;
			}
		}
		if (nonzero) {
			keep_cols.push_back(j);
		}
	}

	NumericMatrix X_matched_diffs(m, keep_cols.size());
	for (int pair_idx = 0; pair_idx < m; ++pair_idx) {
		for (int j = 0; j < static_cast<int>(keep_cols.size()); ++j) {
			X_matched_diffs(pair_idx, j) = X_matched_diffs_full(pair_idx, keep_cols[static_cast<std::size_t>(j)]);
		}
	}

	return List::create(
		_["X_matched_diffs"] = X_matched_diffs,
		_["X_matched_diffs_full"] = X_matched_diffs_full,
		_["X_matched_means_full"] = X_matched_means_full,
		_["yTs_matched"] = yTs_matched,
		_["yCs_matched"] = yCs_matched,
		_["y_matched_diffs"] = y_matched_diffs,
		_["X_reservoir"] = X_reservoir,
		_["y_reservoir"] = y_reservoir,
		_["w_reservoir"] = w_reservoir,
		_["nRT"] = nRT,
		_["nRC"] = nRC,
		_["m"] = m,
		_["d_plus"] = d_plus,
		_["d_minus"] = d_minus,
		_["n11"] = n11,
		_["n10"] = n10,
		_["n01"] = n01,
		_["n00"] = n00
	);
}

// [[Rcpp::export]]
List match_stats_from_indices_cpp(
	const NumericMatrix& X,
	const NumericVector& y,
	const NumericVector& w,
	const IntegerVector& original_m_vec,
	const IntegerVector& i_b,
	int m
) {
	int n_rows = i_b.size();
	int p = X.ncol();
	NumericVector y_sample(n_rows);
	NumericVector w_sample(n_rows);
	NumericMatrix X_sample(n_rows, p);
	IntegerVector m_vec_sample(n_rows); // This will hold the resampled m_vec

	for (int i = 0; i < n_rows; ++i) {
	int idx = i_b[i] - 1;
	y_sample[i] = y[idx];
	w_sample[i] = w[idx];
	m_vec_sample[i] = original_m_vec[idx]; // Corrected resampling
	for (int j = 0; j < p; ++j) {
		X_sample(i, j) = X(idx, j);
	}
	}

	List match_data = compute_zhang_match_data_cpp(
		X_sample,
		y_sample,
		as<IntegerVector>(w_sample),
		m_vec_sample
	);

	return match_data;
}
