context("KK GEE parity against backend implementations")

make_kk_design_for_gee_test <- function(response_type, n_pairs = 80L, n_single = 40L) {
	n <- 2L * n_pairs + n_single
	X <- data.frame(
		x1 = rnorm(n),
		x2 = rnorm(n),
		x3 = rbinom(n, 1L, 0.5)
	)
	des <- DesignSeqOneByOneKK14$new(n = n, response_type = response_type, verbose = FALSE)
	for (i in seq_len(n)) {
		des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	}
	des$.__enclos_env__$private$m <- c(rep(seq_len(n_pairs), each = 2L), rep(0L, n_single))
	des
}

simulate_kk_response_for_gee_test <- function(des, family = c("binomial", "poisson", "proportion", "ordinal")) {
	family <- match.arg(family)
	priv <- des$.__enclos_env__$private
	w <- priv$w
	X <- as.data.frame(priv$X)
	m <- priv$m
	group_id <- ifelse(m > 0L, m, max(m, 0L) + seq_along(m))
	b <- rnorm(max(group_id), sd = 0.35)
	eta <- 0.25 - 0.4 * w + 0.15 * X[[1]] - 0.1 * X[[2]] + 0.2 * X[[3]] + b[group_id]
	if (family == "binomial") {
		return(rbinom(length(eta), size = 1L, prob = plogis(eta)))
	}
	if (family == "poisson") {
		return(rpois(length(eta), lambda = exp(pmin(eta, 5))))
	}
	if (family == "proportion") {
		return(as.numeric(rbinom(length(eta), size = 1L, prob = plogis(eta))))
	}
	latent <- eta + rnorm(length(eta), sd = 0.8)
	cuts <- stats::quantile(latent, probs = c(0.33, 0.66))
	as.integer(cut(latent, breaks = c(-Inf, cuts, Inf), labels = FALSE, right = TRUE))
}

compare_kk_gee_wrapper_paths <- function(class_name, des, use_rcpp_tolerance_est, use_rcpp_tolerance_se) {
	inf_fast <- get(class_name)$new(des, use_rcpp = TRUE, verbose = FALSE)
	inf_ref  <- get(class_name)$new(des, use_rcpp = FALSE, verbose = FALSE)

	est_fast <- inf_fast$compute_estimate()
	est_ref  <- inf_ref$compute_estimate()
	p_fast <- inf_fast$compute_asymp_two_sided_pval()
	p_ref  <- inf_ref$compute_asymp_two_sided_pval()
	se_fast <- inf_fast$.__enclos_env__$private$cached_values$s_beta_hat_T
	se_ref  <- inf_ref$.__enclos_env__$private$cached_values$s_beta_hat_T

	expect_equal(est_fast, est_ref, tolerance = use_rcpp_tolerance_est)
	if (is.finite(se_ref) && is.finite(p_ref)) {
		expect_equal(se_fast, se_ref, tolerance = use_rcpp_tolerance_se)
		expect_equal(p_fast, p_ref, tolerance = 5e-3)
	} else {
		expect_true(inf_ref$is_nonestimable("se"))
		expect_true(is.finite(se_fast))
		expect_true(is.finite(p_fast))
	}
}

test_that("fast KK GEE direct solver matches geepack for binomial and Poisson", {
	skip_if_not_installed("geepack")
	set.seed(20260423)

	generate_direct_data <- function(family = c("binomial", "poisson"), n_pairs = 120L, n_single = 60L) {
		family <- match.arg(family)
		group_id <- c(rep(seq_len(n_pairs), each = 2L), rep(n_pairs + seq_len(n_single), each = 1L))
		n <- length(group_id)
		X <- cbind("(Intercept)" = 1, w = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
		b <- rnorm(max(group_id), sd = 0.45)
		eta <- drop(X %*% c(0.4, -0.35, 0.2, -0.15) + b[group_id])
		y <- if (family == "binomial") {
			rbinom(n, size = 1L, prob = plogis(eta))
		} else {
			rpois(n, lambda = exp(pmin(eta, 5)))
		}
		list(X = X, y = y, group_id = group_id)
	}

	for (family in c("binomial", "poisson")) {
		dat <- generate_direct_data(family)
		df <- data.frame(y = dat$y, w = dat$X[, "w"], x1 = dat$X[, "x1"], x2 = dat$X[, "x2"])
		mod_ref <- geepack::geeglm(
			y ~ w + x1 + x2,
			data = df,
			id = dat$group_id,
			family = if (family == "binomial") stats::binomial("logit") else stats::poisson("log"),
			corstr = "exchangeable",
			std.err = "san.se"
		)
		mod_fast <- EDI:::gee_pairs_singletons_cpp(dat$X, dat$y, dat$group_id, family_str = family)

		beta_tol <- if (family == "binomial") 1e-4 else 6e-3
		se_tol <- if (family == "binomial") 1e-5 else 6e-4
		alpha_tol <- if (family == "binomial") 1e-3 else 6e-2

		expect_equal(as.numeric(mod_fast$beta), as.numeric(stats::coef(mod_ref)), tolerance = beta_tol)
		expect_equal(as.numeric(sqrt(diag(mod_fast$vcov))), as.numeric(sqrt(diag(stats::vcov(mod_ref)))), tolerance = se_tol)
		expect_equal(as.numeric(mod_fast$alpha), as.numeric(mod_ref$geese$alpha[1]), tolerance = alpha_tol)
	}
})

test_that("incidence, count, and proportion KK GEE wrappers match backend fits", {
	skip_if_not_installed("geepack")
	set.seed(20260423)

	des_incid <- make_kk_design_for_gee_test("incidence")
	des_incid$add_all_subject_responses(simulate_kk_response_for_gee_test(des_incid, "binomial"))
	compare_kk_gee_wrapper_paths("InferenceIncidKKGEE", des_incid, use_rcpp_tolerance_est = 1e-4, use_rcpp_tolerance_se = 1e-5)

	des_count <- make_kk_design_for_gee_test("count")
	des_count$add_all_subject_responses(simulate_kk_response_for_gee_test(des_count, "poisson"))
	compare_kk_gee_wrapper_paths("InferenceCountPoissonKKGEE", des_count, use_rcpp_tolerance_est = 2e-2, use_rcpp_tolerance_se = 2e-3)

	des_prop <- make_kk_design_for_gee_test("proportion")
	des_prop$add_all_subject_responses(simulate_kk_response_for_gee_test(des_prop, "proportion"))
	compare_kk_gee_wrapper_paths("InferencePropKKGEE", des_prop, use_rcpp_tolerance_est = 2e-2, use_rcpp_tolerance_se = 2e-2)
})

test_that("ordinal KK GEE wrapper matches direct multgee backend fit", {
	skip_if_not_installed("multgee")
	set.seed(20260423)

	des <- make_kk_design_for_gee_test("ordinal")
	des$add_all_subject_responses(simulate_kk_response_for_gee_test(des, "ordinal"))

	inf <- InferenceOrdinalKKGEE$new(des, verbose = FALSE)
	est <- inf$compute_estimate()
	inf$compute_asymp_two_sided_pval()
	se <- inf$.__enclos_env__$private$cached_values$s_beta_hat_T

	priv <- inf$.__enclos_env__$private
	m_vec <- priv$m
	m_vec[is.na(m_vec)] <- 0L
	group_id <- m_vec
	reservoir_idx <- which(group_id == 0L)
	if (length(reservoir_idx) > 0L) {
		group_id[reservoir_idx] <- max(group_id) + seq_along(reservoir_idx)
	}

	pred_df <- priv$gee_predictors_df()
	dat <- data.frame(y = factor(priv$y, ordered = TRUE), pred_df, group_id = group_id)
	dat <- dat[order(dat$group_id), , drop = FALSE]
	id_sorted <- dat$group_id
	formula_gee <- stats::as.formula(paste("y ~", paste(setdiff(colnames(dat), c("y", "group_id")), collapse = " + ")))
	mod_ref <- suppressWarnings(multgee::ordLORgee(formula_gee, data = dat, id = id_sorted, LORstr = "uniform", link = "logit"))
	beta_ref <- stats::coef(mod_ref)
	j_treat <- priv$gee_treatment_index(beta_ref)
	se_ref <- sqrt(as.numeric(stats::vcov(mod_ref)[j_treat, j_treat]))

	expect_equal(est, as.numeric(beta_ref[j_treat]), tolerance = 1e-8)
	expect_equal(se, se_ref, tolerance = 1e-8)
})
