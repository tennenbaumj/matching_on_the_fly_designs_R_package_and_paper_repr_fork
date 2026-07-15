library(EDI)

expect_error_contains = function(expr, pattern) {
	err = tryCatch(
		{
			force(expr)
			NULL
		},
		error = function(e) conditionMessage(e)
	)
	if (is.null(err) || !grepl(pattern, err, fixed = TRUE)) {
		stop("Expected error containing '", pattern, "'. Got: ", err %||% "<no error>", call. = FALSE)
	}
	invisible(TRUE)
}

build_design = function(response_type, design_cls, X, y_fun) {
	des = design_cls$new(response_type = response_type, n = nrow(X), design_formula = ~ .)
	for (i in seq_len(nrow(X))) {
		w_i = des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		y_i = y_fun(i, w_i)
		des$add_one_subject_response(i, y_i, dead = 1)
	}
	des
}

set.seed(20260423)
n = 120
p = 3
X = as.data.frame(matrix(rnorm(n * p), n, p))
names(X) = paste0("x", seq_len(p))

# Wald-only path should reject both non-auto preferences.
des_cont = build_design("continuous", DesignSeqOneByOneBernoulli, X, function(i, w_i) {
	0.4 * w_i + X$x1[i] + rnorm(1, sd = 0.2)
})
inf_cont = InferenceAllSimpleMeanDiff$new(des_cont)
stopifnot(identical(inf_cont$get_supported_information_preferences(), "auto"))
expect_error_contains(inf_cont$set_information_preference("observed"), "does not support information_preference")
expect_error_contains(inf_cont$set_information_preference("fisher"), "does not support information_preference")

# Observed-only likelihood-backed path should reject Fisher but allow observed.
des_kk = build_design("incidence", DesignSeqOneByOneKK14, X, function(i, w_i) {
	eta = -0.1 + 0.35 * w_i + 0.2 * X$x1[i]
	rbinom(1, 1, plogis(eta))
})
inf_kk = InferenceIncidKKCondLogitPlusGLMMOneLik$new(des_kk, model_formula = ~ .)
stopifnot(identical(inf_kk$get_supported_information_preferences(), c("auto", "observed")))
expect_error_contains(inf_kk$set_information_preference("fisher"), "does not support information_preference")
inf_kk$set_information_preference("observed")
stopifnot(identical(inf_kk$get_information_preference(), "observed"))
stopifnot(is.null(inf_kk$get_information_source_used()))

# Fisher-capable path should allow all three preferences.
des_logit = build_design("incidence", DesignSeqOneByOneBernoulli, X, function(i, w_i) {
	eta = -0.2 + 0.5 * w_i + 0.15 * X$x1[i]
	rbinom(1, 1, plogis(eta))
})
inf_logit = InferenceIncidLogRegr$new(des_logit, model_formula = ~ .)
stopifnot(identical(inf_logit$get_supported_information_preferences(), c("auto", "observed", "fisher")))
inf_logit$set_information_preference("fisher")
stopifnot(identical(inf_logit$get_information_preference(), "fisher"))
stopifnot(is.null(inf_logit$get_information_source_used()))
inf_logit$set_information_preference("observed")
stopifnot(identical(inf_logit$get_information_preference(), "observed"))
stopifnot(is.null(inf_logit$get_information_source_used()))

cat("information_preference_tests: ok\n")
