library(testthat)
library(EDI)

test_that("cluster permutation flattening preserves cluster assignments", {
	n_clusters <- 12L
	cluster_size <- 8L
	n <- n_clusters * cluster_size
	clusters <- lapply(seq_len(n_clusters), function(i) {
		as.integer(((i - 1L) * cluster_size + 1L):(i * cluster_size))
	})

	set.seed(82)
	first <- EDI:::generate_permutations_cluster_cpp(n, 400L, 0.4, clusters)
	set.seed(82)
	second <- EDI:::generate_permutations_cluster_cpp(n, 400L, 0.4, clusters)

	expect_identical(first$w_mat, second$w_mat)
	expect_equal(dim(first$w_mat), c(n, 400L))
	expect_true(all(first$w_mat %in% 0:1))
	expect_null(first$m_mat)
	for (indices in clusters) {
		expect_true(all(apply(first$w_mat[indices, , drop = FALSE], 2L, function(x) length(unique(x))) == 1L))
	}
	expect_lt(abs(mean(first$w_mat) - 0.4), 0.03)
})

test_that("cluster permutation indices are validated before simulation", {
	expect_error(
		EDI:::generate_permutations_cluster_cpp(5L, 2L, 0.5, list(c(0L, 1L))),
		"outside 1..n"
	)
	expect_error(
		EDI:::generate_permutations_cluster_cpp(5L, 2L, 0.5, list(c(1L, 6L))),
		"outside 1..n"
	)
})
