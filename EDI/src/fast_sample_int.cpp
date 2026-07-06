#include <Rcpp.h>
#include <random>
#include <cstdint>
#include <limits>
using namespace Rcpp;

namespace {

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
IntegerVector sample_int_replace_cpp(int n, int size) {
	auto rng = make_local_rng();
	const uint64_t un = static_cast<uint64_t>(n);
	IntegerVector result(size);
	for (int i = 0; i < size; ++i) {
		result[i] = 1 + static_cast<int>(bounded_rand(rng, un));
	}
	return result;
}

// [[Rcpp::export]]
IntegerVector resample_group_rows_cpp(const IntegerVector& group_id, int sample_size) {
	const int n = group_id.size();
	if (sample_size < 0) {
		stop("sample_size must be non-negative.");
	}
	if (n == 0 || sample_size == 0) {
		return IntegerVector(0);
	}

	int max_group = 0;
	for (int i = 0; i < n; ++i) {
		const int g = group_id[i];
		if (g == NA_INTEGER || g <= 0) {
			stop("group_id must contain only positive integers.");
		}
		if (g > max_group) {
			max_group = g;
		}
	}

	std::vector<int> counts(max_group, 0);
	for (int i = 0; i < n; ++i) {
		counts[group_id[i] - 1]++;
	}
	for (int g = 0; g < max_group; ++g) {
		if (counts[g] == 0) {
			stop("group_id must be consecutive positive integers starting at 1.");
		}
	}

	std::vector< std::vector<int> > rows_by_group(max_group);
	for (int g = 0; g < max_group; ++g) {
		rows_by_group[g].reserve(counts[g]);
	}
	for (int i = 0; i < n; ++i) {
		rows_by_group[group_id[i] - 1].push_back(i + 1);
	}

	auto rng = make_local_rng();
	const uint64_t u_mg = static_cast<uint64_t>(max_group);

	IntegerVector sampled_groups(sample_size);
	int out_size = 0;
	for (int draw = 0; draw < sample_size; ++draw) {
		const int sampled_group = 1 + static_cast<int>(bounded_rand(rng, u_mg));
		sampled_groups[draw] = sampled_group;
		out_size += rows_by_group[sampled_group - 1].size();
	}

	IntegerVector out(out_size);
	int out_idx = 0;
	for (int draw = 0; draw < sample_size; ++draw) {
		const std::vector<int>& rows = rows_by_group[sampled_groups[draw] - 1];
		for (size_t j = 0; j < rows.size(); ++j) {
			out[out_idx++] = rows[j];
		}
	}
	return out;
}
