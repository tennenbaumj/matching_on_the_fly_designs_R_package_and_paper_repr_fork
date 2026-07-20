make_completed_fixed_design <- function(response_type, x, w, y, dead = NULL) {
	des <- EDI:::DesignFixed$new(n = length(y), response_type = response_type, verbose = FALSE)
	des$add_all_subjects_to_experiment(data.frame(x = x))
	des$overwrite_all_subject_assignments(w)
	if (is.null(dead)) {
		des$add_all_subject_responses(y)
	} else {
		des$add_all_subject_responses(y, dead)
	}
	des
}

test_that("smart_cold_start_default TRUE and FALSE agree across core optimization families", {
	set.seed(101)
	n <- 80
	x <- rnorm(n)
	w <- rep(c(1, -1), length.out = n)
	w01 <- (w + 1) / 2

	y_logit <- rbinom(n, 1, plogis(-0.3 + 0.8 * w01 + 0.5 * x))
	des_logit <- make_completed_fixed_design("incidence", x, w, y_logit)
	inf_logit_s <- InferenceIncidLogRegr$new(des_logit, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_logit_l <- InferenceIncidLogRegr$new(des_logit, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_logit_s$compute_estimate(), inf_logit_l$compute_estimate(), tolerance = 5e-3)

	y_probit <- rbinom(n, 1, pnorm(-0.3 + 0.8 * w01 + 0.5 * x))
	des_probit <- make_completed_fixed_design("incidence", x, w, y_probit)
	inf_probit_inc_s <- InferenceIncidProbitRegr$new(des_probit, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_probit_inc_l <- InferenceIncidProbitRegr$new(des_probit, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_probit_inc_s$compute_estimate(), inf_probit_inc_l$compute_estimate(), tolerance = 5e-3)

	y_pois <- rpois(n, lambda = exp(0.2 + 0.3 * w01 - 0.2 * x))
	des_pois <- make_completed_fixed_design("count", x, w, y_pois)
	inf_pois_s <- InferenceCountPoisson$new(des_pois, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_pois_l <- InferenceCountPoisson$new(des_pois, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_pois_s$compute_estimate(), inf_pois_l$compute_estimate(), tolerance = 5e-3)

	y_nb <- MASS::rnegbin(n, mu = exp(0.4 + 0.25 * w01 + 0.1 * x), theta = 2)
	des_nb <- make_completed_fixed_design("count", x, w, y_nb)
	inf_nb_s <- InferenceCountNegBin$new(des_nb, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_nb_l <- InferenceCountNegBin$new(des_nb, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_nb_s$compute_estimate(), inf_nb_l$compute_estimate(), tolerance = 5e-3)

	y_surv <- exp(1 + 0.3 * w01 - 0.15 * x + rnorm(n, sd = 0.15))
	dead <- rbinom(n, 1, 0.85)
	des_surv <- make_completed_fixed_design("survival", x, w, y_surv, dead = dead)
	inf_surv_s <- InferenceSurvivalWeibullRegr$new(des_surv, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_surv_l <- InferenceSurvivalWeibullRegr$new(des_surv, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_surv_s$compute_estimate(), inf_surv_l$compute_estimate(), tolerance = 5e-3)

	eta_ord <- 0.6 * w01 - 0.25 * x + rnorm(n)
	y_ord <- as.numeric(cut(eta_ord, breaks = c(-Inf, -0.4, 0.5, Inf), labels = FALSE))
	des_ord <- make_completed_fixed_design("ordinal", x, w, y_ord)
	inf_ord_s <- InferenceOrdinalPropOddsRegr$new(des_ord, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_ord_l <- InferenceOrdinalPropOddsRegr$new(des_ord, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_ord_s$compute_estimate(), inf_ord_l$compute_estimate(), tolerance = 5e-3)

	inf_probit_s <- InferenceOrdinalOrderedProbitRegr$new(des_ord, verbose = FALSE, smart_cold_start_default = TRUE)
	inf_probit_l <- InferenceOrdinalOrderedProbitRegr$new(des_ord, verbose = FALSE, smart_cold_start_default = FALSE)
	expect_equal(inf_probit_s$compute_estimate(), inf_probit_l$compute_estimate(), tolerance = 5e-3)
})
