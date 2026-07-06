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
IntegerMatrix bootstrap_indices_cpp(int n, int B) {
	auto rng = make_local_rng();
	const uint64_t un = static_cast<uint64_t>(n);
	IntegerMatrix idx(B, n);
	for (int i = 0; i < B; ++i) {
		for (int j = 0; j < n; ++j) {
			idx(i, j) = 1 + static_cast<int>(bounded_rand(rng, un));
		}
	}
	return idx;
}
