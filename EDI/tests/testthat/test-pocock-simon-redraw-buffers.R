library(testthat)
library(EDI)

pocock_simon_redraw_reference <- function(level_rows, num_levels_total, weights,
		p_best, prob_treatment) {
	counts <- matrix(0, num_levels_total, 2L)
	assignments <- integer(nrow(level_rows))
	for (i in seq_len(nrow(level_rows))) {
		imbalance <- vapply(0:1, function(treatment) {
			sum(vapply(seq_len(ncol(level_rows)), function(j) {
				row <- level_rows[i, j]
				after <- counts[row, ] + as.numeric(0:1 == treatment)
				weights[j] * sum((after - mean(after))^2)
			}, numeric(1L)))
		}, numeric(1L))
		if (imbalance[1L] == imbalance[2L]) {
			assigned <- as.integer(runif(1L) < prob_treatment)
		} else {
			best <- if (imbalance[2L] < imbalance[1L]) 1L else 0L
			assigned <- if (runif(1L) < p_best) best else 1L - best
		}
		assignments[i] <- assigned
		for (row in level_rows[i, ]) counts[row, assigned + 1L] <- counts[row, assigned + 1L] + 1
	}
	assignments
}

test_that("Pocock-Simon redraw buffers match an independent R reference", {
	set.seed(84)
	n <- 80L
	level_rows <- cbind(
		sample(1:2, n, TRUE),
		sample(3:5, n, TRUE),
		sample(6:7, n, TRUE)
	)
	storage.mode(level_rows) <- "integer"
	weights <- c(1, 0.75, 1.25)

	set.seed(841)
	expected <- pocock_simon_redraw_reference(level_rows, 7L, weights, 0.7, 0.5)
	set.seed(841)
	actual <- EDI:::pocock_simon_redraw_w_cpp(level_rows, 7L, weights, 0.7, 0.5)

	expect_identical(as.integer(actual), expected)
	expect_length(actual, n)
	expect_true(all(actual %in% 0:1))
})

test_that("Pocock-Simon redraw input rows are validated", {
	valid <- matrix(as.integer(c(1, 2, 3, 4)), nrow = 2L)
	expect_error(
		EDI:::pocock_simon_redraw_w_cpp(valid - 1L, 4L, c(1, 1), 0.7, 0.5),
		"outside 1..num_levels_total"
	)
	expect_error(
		EDI:::pocock_simon_redraw_w_cpp(valid, 4L, 1, 0.7, 0.5),
		"weights length"
	)
})
