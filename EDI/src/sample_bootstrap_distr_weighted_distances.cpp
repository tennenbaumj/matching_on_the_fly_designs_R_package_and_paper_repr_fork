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
NumericVector compute_bootstrapped_weighted_sqd_distances_cpp(
	NumericMatrix X_all_scaled_col_subset,
	NumericVector covariate_weights,
	int t, // self$t
	int B) { // private$other_params$num_boot

	int d = covariate_weights.size();
	NumericVector bootstrapped_weighted_sqd_distances(B);

	auto rng = make_local_rng();
	const uint64_t ut = static_cast<uint64_t>(t);

	for (int b = 0; b < B; ++b) {
		int i1 = static_cast<int>(bounded_rand(rng, ut));
		int i2 = static_cast<int>(bounded_rand(rng, ut));
		if (t > 1) {
			while (i1 == i2) {
				i2 = static_cast<int>(bounded_rand(rng, ut));
			}
		}

		double sqd_weighted_sum = 0.0;
		for (int j = 0; j < d; ++j) {
			double delta = X_all_scaled_col_subset(i1, j) - X_all_scaled_col_subset(i2, j);
			sqd_weighted_sum += delta * delta * covariate_weights[j];
		}

		bootstrapped_weighted_sqd_distances[b] = sqd_weighted_sum;
	}

	return bootstrapped_weighted_sqd_distances;
}
