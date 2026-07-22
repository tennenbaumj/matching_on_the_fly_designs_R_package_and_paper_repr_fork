#' Likelihood-Backed Asymptotic Inference
#'
#' @name InferenceAsympLik
#' @description Intermediate base class for asymptotic inference families that expose
#' likelihood / partial-likelihood / working-likelihood test paths in addition
#' to Wald inference. The term "likelihood" is used broadly: subclasses may be
#' backed by a true full likelihood, a partial likelihood (e.g. Cox PH), a
#' quasi-likelihood (e.g. GEE, quasi-Poisson), or a composite/combined
#' likelihood. Classes requiring a full generative likelihood — i.e. those
#' supporting parametric-bootstrap LR calibration — inherit instead from
#' \code{InferenceParamBootstrap}.
#'
#' @keywords internal
InferenceAsympLik = R6::R6Class("InferenceAsympLik",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMSummaryTable,
	public = list(
		#' @description Computes an asymptotic confidence interval using the configured test.
		#'
		#' @param alpha Significance level 1 - \code{alpha}. Default 0.05.
		#'
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_confidence_interval_impl(alpha),
				score = private$compute_score_confidence_interval_impl(alpha),
				gradient = private$compute_gradient_confidence_interval_impl(alpha),
				lik_ratio = private$compute_lik_ratio_confidence_interval_impl(alpha),
				lik_ratio_bartlett_approx = private$compute_lik_ratio_bartlett_approx_confidence_interval_impl(alpha),
				lik_ratio_bartlett_exact = private$compute_lik_ratio_bartlett_exact_confidence_interval_impl(alpha)
			)
		},
		#' @description Computes an asymptotic two-sided p-value using the configured test.
		#'
		#' @param delta Null treatment effect to test against. Default 0.
		#'
		#' @return The asymptotic p-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			switch(
				private$testing_type,
				wald = private$compute_wald_two_sided_pval_impl(delta),
				score = private$compute_score_two_sided_pval_impl(delta),
				gradient = private$compute_gradient_two_sided_pval_impl(delta),
				lik_ratio = private$compute_lik_ratio_two_sided_pval_impl(delta),
				lik_ratio_bartlett_approx = private$compute_lik_ratio_bartlett_approx_two_sided_pval_impl(delta),
				lik_ratio_bartlett_exact = private$compute_lik_ratio_bartlett_exact_two_sided_pval_impl(delta)
			)
		},
		#' @description Sets the asymptotic testing method used by p-values and CIs.
		#'
		#' @param testing_type One of \code{"wald"}, \code{"score"}, \code{"gradient"}, \code{"lik_ratio"}, \code{"lik_ratio_bartlett_approx"}, or \code{"lik_ratio_bartlett_exact"}.
		#'
		#' @return The inference object, invisibly.
		set_testing_type = function(testing_type = c("wald", "score", "gradient", "lik_ratio", "lik_ratio_bartlett_approx", "lik_ratio_bartlett_exact")){
			testing_type = private$normalize_testing_type(testing_type)
			supported = private$get_supported_testing_types_with_bartlett()
			if (!testing_type %in% supported) {
				stop(
					class(self)[1], " does not support testing_type = \"", testing_type,
					"\". Supported values are: ", paste(supported, collapse = ", "),
					call. = FALSE
				)
			}
			private$testing_type = testing_type
			invisible(self)
		},
		#' @description Sets the information matrix preference used by score-test dispatch.
		#'
		#' @param information_preference One of \code{"auto"}, \code{"fisher"}, or \code{"observed"}.
		#'
		#' @return The inference object, invisibly.
		set_information_preference = function(information_preference = c("auto", "fisher", "observed")){
			information_preference = private$normalize_information_preference(information_preference)
			supported = private$get_supported_information_preferences_impl()
			if (!information_preference %in% supported) {
				stop(
					class(self)[1], " does not support information_preference = \"", information_preference,
					"\". Supported values are: ", paste(supported, collapse = ", "),
					call. = FALSE
				)
			}
			private$information_preference = information_preference
			private$information_source_used = NULL
			private$clear_likelihood_test_eval_cache()
			invisible(self)
		},
		#' @description Gets the asymptotic testing method used by p-values and CIs.
		get_testing_type = function(){
			private$testing_type
		},
		#' @description Gets the score-test information matrix preference.
		get_information_preference = function(){
			private$information_preference
		},
		#' @description Gets the actual information source used by the most recent information-backed computation.
		get_information_source_used = function(){
			private$information_source_used
		},
		#' @description Gets the asymptotic testing methods supported by this inference object.
		get_supported_testing_types = function(){
			private$get_supported_testing_types_with_bartlett()
		},
		#' @description Gets the score-test information matrix preferences supported by this inference object.
		get_supported_information_preferences = function(){
			private$get_supported_information_preferences_impl()
		},
		#' @description Computes the score two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_score_two_sided_pval = function(delta = 0){
			private$compute_score_two_sided_pval_impl(delta)
		},
		#' @description Computes the score confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_score_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_score_confidence_interval_impl(alpha)
		},

		#' @description Computes the likelihood-ratio two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			private$compute_lik_ratio_two_sided_pval_impl(delta)
		},

		#' @description Computes the likelihood-ratio confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_lik_ratio_confidence_interval_impl(alpha)
		},

		#' @description Computes the approximate (Monte-Carlo) Bartlett-corrected likelihood-ratio
		#' two-sided p-value regardless of configured testing type. Returns \code{NA_real_}
		#' for subclasses that do not implement an approximate Bartlett correction factor.
		#'
		#' The approximate Bartlett factor is estimated by Monte Carlo (e.g. the generic
		#' \code{InferenceParamBootstrap} factor): \code{B} datasets are simulated under
		#' the null-restricted fit at \code{delta} and refit to approximate
		#' \code{E[LR | H0]}, the quantity a classical analytic Bartlett correction
		#' targets exactly. The Monte-Carlo draws are seeded from this object's own
		#' \code{seed} (see \code{set_seed()}), so repeated calls at the same
		#' \code{delta} with the same \code{B} are reproducible; there is no separate
		#' seed argument here.
		#'
		#' See \code{compute_lik_ratio_bartlett_exact_two_sided_pval()} for the
		#' closed-form analytic counterpart (no simulation, no \code{B}).
		#'
		#' @param delta Null treatment effect. Default 0.
		#' @param B Number of Monte-Carlo replicates used to estimate the Bartlett
		#'   factor. Default 99.
		compute_lik_ratio_bartlett_approx_two_sided_pval = function(delta = 0, B = 99){
			private$compute_lik_ratio_bartlett_approx_two_sided_pval_impl(delta, B = B)
		},

		#' @description Computes the approximate (Monte-Carlo) Bartlett-corrected likelihood-ratio
		#' confidence interval regardless of configured testing type. Returns
		#' \code{c(NA_real_, NA_real_)} for subclasses that do not implement an
		#' approximate Bartlett correction factor.
		#'
		#' See \code{compute_lik_ratio_bartlett_approx_two_sided_pval()} for what
		#' \code{B} controls and how the Monte-Carlo seed is inherited from this
		#' object's own \code{seed}. Each p-value evaluation during the
		#' confidence-interval search re-simulates \code{B} replicates, so this can
		#' be substantially more expensive than the p-value alone.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of Monte-Carlo replicates used to estimate the Bartlett
		#'   factor. Default 99.
		compute_lik_ratio_bartlett_approx_confidence_interval = function(alpha = 0.05, B = 99){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_lik_ratio_bartlett_approx_confidence_interval_impl(alpha, B = B)
		},

		#' @description Computes the exact (closed-form analytic) Bartlett-corrected
		#' likelihood-ratio two-sided p-value regardless of configured testing type.
		#' Returns \code{NA_real_} for subclasses that do not implement an exact,
		#' bespoke analytic Bartlett correction factor (no family does yet; see
		#' \code{get_bartlett_factor_exact()}).
		#'
		#' Unlike \code{compute_lik_ratio_bartlett_approx_two_sided_pval()}, this path
		#' involves no simulation and no Monte-Carlo replicate count.
		#'
		#' @param delta Null treatment effect. Default 0.
		compute_lik_ratio_bartlett_exact_two_sided_pval = function(delta = 0){
			private$compute_lik_ratio_bartlett_exact_two_sided_pval_impl(delta)
		},

		#' @description Computes the exact (closed-form analytic) Bartlett-corrected
		#' likelihood-ratio confidence interval regardless of configured testing type.
		#' Returns \code{c(NA_real_, NA_real_)} for subclasses that do not implement an
		#' exact, bespoke analytic Bartlett correction factor.
		#'
		#' @param alpha Significance level. Default 0.05.
		compute_lik_ratio_bartlett_exact_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_lik_ratio_bartlett_exact_confidence_interval_impl(alpha)
		},

		#' @description Computes "the best available" Bartlett-corrected likelihood-ratio
		#' two-sided p-value regardless of configured testing type: uses the exact
		#' (closed-form analytic) factor if this class implements one, otherwise
		#' falls back to the approximate (Monte-Carlo) factor. Errors if the class
		#' supports neither (see \code{supports_bartlett_likelihood_ratio_exact()}/
		#' \code{supports_bartlett_likelihood_ratio_approx()}).
		#'
		#' This is a convenience entry point for callers who want a Bartlett-corrected
		#' p-value without caring which mechanism produced it. \strong{Because exact and
		#' approximate factors are computed differently} (deterministic closed form vs.
		#' seeded Monte-Carlo simulation), \strong{the same call can silently start
		#' returning different numeric results on a future package version} once a
		#' family gains an exact implementation where previously only the
		#' approximate path existed. Callers who need results stable across package
		#' versions (e.g. for reproducibility or regression tests) should call
		#' \code{compute_lik_ratio_bartlett_approx_two_sided_pval()} or
		#' \code{compute_lik_ratio_bartlett_exact_two_sided_pval()} directly instead.
		#'
		#' @param delta Null treatment effect. Default 0.
		#' @param B Number of Monte-Carlo replicates, used only when the exact factor
		#'   is unavailable and the approximate factor is used instead. If explicitly
		#'   supplied but the exact factor is used (so \code{B} has no effect), a
		#'   warning is issued; \code{B} left at its default is silently ignored in
		#'   that case. Default 99.
		compute_lik_ratio_bartlett_two_sided_pval = function(delta = 0, B = 99){
			private$compute_lik_ratio_bartlett_two_sided_pval_impl(delta, B = B, B_missing = missing(B))
		},

		#' @description Computes "the best available" Bartlett-corrected likelihood-ratio
		#' confidence interval regardless of configured testing type: uses the exact
		#' (closed-form analytic) factor if this class implements one, otherwise
		#' falls back to the approximate (Monte-Carlo) factor. Errors if the class
		#' supports neither.
		#'
		#' See \code{compute_lik_ratio_bartlett_two_sided_pval()} for the exact-over-approx
		#' selection rule, the \code{B}-ignored warning behavior, and why callers who
		#' need version-to-version reproducibility should prefer the explicit
		#' \code{_approx}/\code{_exact} methods instead.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of Monte-Carlo replicates, used only when the exact factor
		#'   is unavailable and the approximate factor is used instead. Default 99.
		compute_lik_ratio_bartlett_confidence_interval = function(alpha = 0.05, B = 99){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_lik_ratio_bartlett_confidence_interval_impl(alpha, B = B, B_missing = missing(B))
		},

		#' @description Computes the gradient two-sided p-value regardless of configured testing type.
		#' @param delta Null treatment effect.
		compute_gradient_two_sided_pval = function(delta = 0){
			private$compute_gradient_two_sided_pval_impl(delta)
		},

		#' @description Computes the gradient confidence interval regardless of configured testing type.
		#' @param alpha Significance level. Default 0.05.
		compute_gradient_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$compute_gradient_confidence_interval_impl(alpha)
		}
	),
	private = c(InferenceMixinCIInversion$private, InferenceMixinInformationMatrix$private, InferenceMixinLikelihoodTestMemoization$private, list(
		likelihood_ci_max_abs = 10,
		testing_type = "wald",
		information_preference = "auto",
		information_source_used = NULL,

		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_score_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			private$invert_gradient_ci_uniroot(alpha)
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			private$invert_lik_ratio_ci_newton(alpha)
		},
		compute_lik_ratio_bartlett_approx_two_sided_pval_impl = function(delta, B = 99){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio_bartlett_approx", bartlett_B = B)
		},
		compute_lik_ratio_bartlett_approx_confidence_interval_impl = function(alpha, B = 99){
			private$invert_test_pval_confidence_interval(alpha, testing_type = "lik_ratio_bartlett_approx", bartlett_B = B)
		},
		compute_lik_ratio_bartlett_exact_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio_bartlett_exact")
		},
		compute_lik_ratio_bartlett_exact_confidence_interval_impl = function(alpha){
			private$invert_test_pval_confidence_interval(alpha, testing_type = "lik_ratio_bartlett_exact")
		},
		warn_bartlett_B_ignored_by_exact = function(){
			warning(
				class(self)[1], " has an exact Bartlett correction factor for this class; ",
				"B is ignored (no simulation is performed).",
				call. = FALSE
			)
		},
		stop_bartlett_unsupported = function(){
			stop(
				class(self)[1], " does not support Bartlett-corrected likelihood-ratio inference ",
				"(neither an exact nor an approximate factor is implemented). ",
				"See supports_bartlett_likelihood_ratio_exact() / supports_bartlett_likelihood_ratio_approx().",
				call. = FALSE
			)
		},
		compute_lik_ratio_bartlett_two_sided_pval_impl = function(delta, B = 99, B_missing = TRUE){
			if (isTRUE(private$supports_bartlett_likelihood_ratio_exact())) {
				if (!B_missing) private$warn_bartlett_B_ignored_by_exact()
				return(private$compute_lik_ratio_bartlett_exact_two_sided_pval_impl(delta))
			}
			if (isTRUE(private$supports_bartlett_likelihood_ratio_approx())) {
				return(private$compute_lik_ratio_bartlett_approx_two_sided_pval_impl(delta, B = B))
			}
			private$stop_bartlett_unsupported()
		},
		compute_lik_ratio_bartlett_confidence_interval_impl = function(alpha, B = 99, B_missing = TRUE){
			if (isTRUE(private$supports_bartlett_likelihood_ratio_exact())) {
				if (!B_missing) private$warn_bartlett_B_ignored_by_exact()
				return(private$compute_lik_ratio_bartlett_exact_confidence_interval_impl(alpha))
			}
			if (isTRUE(private$supports_bartlett_likelihood_ratio_approx())) {
				return(private$compute_lik_ratio_bartlett_approx_confidence_interval_impl(alpha, B = B))
			}
			private$stop_bartlett_unsupported()
		},


		get_likelihood_test_spec = function(){
			NULL
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		supports_information_preference = function(){
			isTRUE(private$supports_likelihood_tests())
		},

		supports_observed_information = function(){
			isTRUE(private$supports_information_preference())
		},

		supports_fisher_information = function(){
			FALSE
		},

		#' Whether this class exposes an approximate (Monte-Carlo) Bartlett-corrected
		#' likelihood-ratio test/CI via get_bartlett_factor_approx(). Default FALSE;
		#' InferenceParamBootstrap overrides this to delegate to
		#' supports_lik_ratio_param_bootstrap().
		supports_bartlett_likelihood_ratio_approx = function(){
			FALSE
		},

		#' Monte-Carlo estimate of the Bartlett correction factor c(delta), such that
		#' LR_B(delta) = LR(delta) / c(delta) is referred to chi-square(1). B controls
		#' the number of null-simulated replicates (ignored by classes that don't
		#' simulate). Default NULL (unimplemented).
		get_bartlett_factor_approx = function(spec, delta, full_fit, null_fit, B = 99){
			NULL
		},

		#' Whether this class exposes an exact (closed-form analytic) Bartlett-corrected
		#' likelihood-ratio test/CI via get_bartlett_factor_exact(). Default FALSE;
		#' no family implements this yet -- individual families opt in once their
		#' bespoke analytic factor is derived.
		supports_bartlett_likelihood_ratio_exact = function(){
			FALSE
		},

		#' Closed-form analytic Bartlett correction factor c(delta) (e.g. a
		#' Cordeiro-style GLM correction). No simulation, no replicate count.
		#' Default NULL (unimplemented).
		get_bartlett_factor_exact = function(spec, delta, full_fit, null_fit){
			NULL
		},

		get_supported_testing_types_impl = function(){
			if (isTRUE(private$supports_likelihood_tests())) {
				c("wald", "score", "gradient", "lik_ratio")
			} else {
				"wald"
			}
		},

		#' Wraps get_supported_testing_types_impl() to conditionally append the
		#' Bartlett testing types. Kept as a separate wrapper (rather than folded
		#' into get_supported_testing_types_impl() itself) because many concrete
		#' and abstract subclasses across the codebase override
		#' get_supported_testing_types_impl() directly with a hard-coded vector;
		#' appending here means those overrides don't each need to be touched to
		#' pick up Bartlett support.
		get_supported_testing_types_with_bartlett = function(){
			types = private$get_supported_testing_types_impl()
			if (isTRUE(private$supports_bartlett_likelihood_ratio_approx())) {
				types = c(types, "lik_ratio_bartlett_approx")
			}
			if (isTRUE(private$supports_bartlett_likelihood_ratio_exact())) {
				types = c(types, "lik_ratio_bartlett_exact")
			}
			unique(types)
		},

		get_supported_information_preferences_impl = function(){
			if (isTRUE(private$supports_information_preference())) {
				c("auto", "observed")
			} else {
				"auto"
			}
		},

		get_default_information_source = function(){
			if (isTRUE(private$supports_fisher_information())) return("fisher")
			if (isTRUE(private$supports_observed_information())) return("observed")
			"legacy"
		},

		normalize_testing_type = function(testing_type){
			if (length(testing_type) != 1L) testing_type = testing_type[1L]
			testing_type = tolower(as.character(testing_type))
			switch(
				testing_type,
				wald = "wald",
				score = "score",
				gradient = "gradient",
				lr = "lik_ratio",
				lrt = "lik_ratio",
				lik_ratio = "lik_ratio",
				likelihood_ratio = "lik_ratio",
				lik_ratio_bartlett_approx = "lik_ratio_bartlett_approx",
				lr_bartlett_approx = "lik_ratio_bartlett_approx",
				lrb_approx = "lik_ratio_bartlett_approx",
				bartlett_approx = "lik_ratio_bartlett_approx",
				lik_ratio_bartlett_exact = "lik_ratio_bartlett_exact",
				lr_bartlett_exact = "lik_ratio_bartlett_exact",
				lrb_exact = "lik_ratio_bartlett_exact",
				bartlett_exact = "lik_ratio_bartlett_exact",
				stop("testing_type must be one of: wald, score, gradient, lik_ratio, lik_ratio_bartlett_approx, lik_ratio_bartlett_exact", call. = FALSE)
			)
		},
		normalize_information_preference = function(information_preference){
			if (length(information_preference) != 1L) information_preference = information_preference[1L]
			information_preference = tolower(as.character(information_preference))
			switch(
				information_preference,
				auto = "auto",
				fisher = "fisher",
				observed = "observed",
				obs = "observed",
				stop("information_preference must be one of: auto, fisher, observed", call. = FALSE)
			)
		}
	))
)
