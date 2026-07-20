# BRT "smoothed" Fast-Kernel Noise Support — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend every `compute_fast_rand_bootstrap_distr` C++ kernel that operates on real-valued (continuous/survival) responses to accept an optional per-draw, per-row noise matrix, so that `type = "smoothed"` bootstrap-randomization p-values/CIs use the fast C++ path instead of the ~7x-slower R-level per-draw fallback (`InferenceRandBootstrap$approximate_rand_bootstrap_distribution_beta_hat_T`'s reused-worker / per-iteration path).

**Architecture:** The R-level dispatcher (`inference_all_abstract_rand_bootstrap.R`) currently disables *every* class's fast kernel unconditionally whenever a draw carries `smooth_noise`, because none of the C++ kernels can consume it. We move that decision down into each class: `rand_bootstrap_draw_matrices()` now also packs a `noise_mat` (n × B) whenever draws carry `smooth_noise`, and each capable kernel adds `noise_mat[i, b]` to the resampled response before applying the sharp-null shift — exactly mirroring what the R fallback already does in `load_rand_bootstrap_assignment_into_worker()` (`y_sim = y0_full[draw$i_b] + smooth_noise`, shift applied after). Two ordinal classes (Ridit, Jonckheere-Terpstra) operate on integer category codes, where continuous noise isn't semantically meaningful (per the existing roxygen caveat, "smoothed" is "only meaningful for continuous responses"); they explicitly decline noise-carrying draws in their own `compute_fast_rand_bootstrap_distr` and keep using the slow fallback, unchanged from today.

**Tech Stack:** R6 (EDI package), Rcpp/RcppEigen, OpenMP, testthat.

## Global Constraints

- Package root: `EDI/`. Sources: `EDI/src/*.cpp`. R6 classes: `EDI/R/*.R`. Tests: `EDI/tests/testthat/*.R`.
- After changing any `// [[Rcpp::export]]` function signature, regenerate bindings with `Rscript -e 'Rcpp::compileAttributes("EDI")'` **before** installing — otherwise `RcppExports.cpp`/`RcppExports.R` are stale and the build either fails or silently calls the wrong arity.
- Install with `R CMD INSTALL --no-docs EDI` (matches this repo's established convention; see `package_metadata/perf_experiments_final.md`). Do **not** run `roxygen2::roxygenize()` directly for `.Rd`/`NAMESPACE` regeneration — use `Rscript fast_roxygenize.R` from the repo root (only needed in Task 12, which touches roxygen text).
- Every new/changed C++ parameter is `Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat` with **no default value** — it is always passed explicitly from R (either `mats$noise_mat`, which may itself be R `NULL`, or `R_NilValue` from the one internal C++-to-C++ call site in Task 3). C++ forbids a non-defaulted parameter after a defaulted one, and `num_cores` (the trailing parameter in every kernel here) has no default either, so giving `noise_mat` one would gain nothing while risking silent positional mistakes.
- Preserve execution order exactly as today: noise is added to the resampled response **before** the delta (sharp-null) shift is applied, for **every** row (treated and control) — never conditionally on `w`. This matches `load_rand_bootstrap_assignment_into_worker()` in `EDI/R/inference_all_abstract_rand_bootstrap.R`.
- Every kernel task's correctness test uses the same *identity-resampling equivalence* trick (see Task 1's "Interfaces" section and Task 2 for the full pattern) — no per-statistic reference math needed, because `y0[i] + noise[i, b]` under an identity `i_mat` column is bitwise interchangeable with pre-adding the noise to `y0` and calling the kernel with `noise_mat = NULL`.

---

## Survey: the 8 kernel files in scope (and the 2 explicitly out of scope)

| # | R class | R file | C++ export | C++ file |
|---|---|---|---|---|
| 1 | `InferenceAllSimpleWilcox` | `inference_all_simple_wilcox.R` | `compute_wilcox_hl_rand_bootstrap_parallel_cpp` | `fast_wilcox_hl.cpp` |
| 2 | `InferenceAllSimpleMeanDiff` | `inference_all_mean_diff.R` | `compute_rand_bootstrap_mean_diff_parallel_cpp` | `rand_bootstrap_mean_diff_parallel.cpp` |
| 3 | `InferenceContinOLS` | `inference_continuous_ols.R` | `compute_rand_bootstrap_ols_parallel_cpp` | `rand_bootstrap_ols_parallel.cpp` |
| 4 | `InferenceContinRobustRegr` | `inference_continuous_robust_regr.R` | `compute_robust_rand_bootstrap_parallel_cpp` | `fast_robust_regression.cpp` |
| 5 | `InferenceSurvivalCoxPHRegr` | `inference_survival_coxph.R` | `compute_coxph_rand_bootstrap_parallel_cpp` | `fast_coxph_regression.cpp` |
| 6 | `InferenceSurvivalKKWeibullMarginal` | `inference_survival_KK_weibull_marginal.R` | `compute_weibull_rand_bootstrap_parallel_cpp` | `fast_weibull_regression.cpp` |
| 7 | `InferenceSurvivalLogRank` | `inference_survival_log_rank.R` | `compute_logrank_rand_bootstrap_parallel_cpp` | `fast_logrank.cpp` |
| 8 | `InferenceSurvivalRestrictedMeanDiff` + `InferenceSurvivalKMDiff` (share one kernel, `do_rmst` flag) | `inference_survival_rmst.R`, `inference_survival_km_diff.R` | `compute_survival_stat_diff_rand_bootstrap_parallel_cpp` | `fast_survival_stats.cpp` |

**Out of scope (documented, not silently skipped):** `InferenceOrdinalRidit` (`compute_ridit_rand_bootstrap_parallel_cpp`) and `InferenceOrdinalJonckheereTerpstraTest` (`compute_jt_rand_bootstrap_parallel_cpp`) operate on integer ordinal category codes. Adding continuous Gaussian noise to a category code is not statistically meaningful (this matches the existing roxygen: "`smoothed`... Only meaningful for continuous responses"). Task 1 makes both classes explicitly decline noise-carrying draws (rather than relying on the blanket dispatcher check being removed), so they keep using today's slow-but-correct fallback with no behavior change.

---

### Task 1: Foundation — `noise_mat` plumbing, dispatcher cleanup, ordinal opt-out

**Files:**
- Modify: `EDI/R/inference_all_abstract_rand_bootstrap.R:122-129` (site A, inside `approximate_rand_bootstrap_distribution_beta_hat_T`)
- Modify: `EDI/R/inference_all_abstract_rand_bootstrap.R:654-657` (site B, inside `get_brt_distribution_prefix`)
- Modify: `EDI/R/inference_all_abstract_rand_bootstrap.R:482-495` (`rand_bootstrap_draw_matrices`)
- Modify: `EDI/R/inference_ordinal_ridit.R:233-242` (`compute_fast_rand_bootstrap_distr`)
- Modify: `EDI/R/inference_ordinal_jonckheere_terpstra_test.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-noise-mat-plumbing.R`

**Interfaces:**
- Produces: `private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)` now returns `list(i_mat = <n x B IntegerMatrix>, w_mat = <n x B IntegerMatrix>, noise_mat = <n x B NumericMatrix> | NULL)`. `noise_mat` is `NULL` unless `rand_bootstrap_draws[[1]][["smooth_noise"]]` is non-NULL, in which case **every** draw must carry a length-`n` `smooth_noise` vector or the whole call returns `NULL` (same all-or-nothing contract the function already uses for `w_b`).
- Consumes (by Tasks 2-9): each class's `compute_fast_rand_bootstrap_distr` reads `mats$noise_mat` and passes it as the new `noise_mat` argument to its C++ kernel — `NULL` from R becomes `R_NilValue`/`Nullable::isNotNull() == false` in C++ automatically, no special-casing needed on the R side.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-noise-mat-plumbing.R`:

```r
library(EDI)

test_that("rand_bootstrap_draw_matrices packs a noise_mat only when draws carry smooth_noise", {
	set.seed(20260721)
	n = 10
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))
	inf = InferenceAllSimpleWilcox$new(des)
	priv = inf$.__enclos_env__$private

	draws = priv$generate_rand_bootstrap_draws(B = 4L, materialize_w = TRUE)
	mats_no_noise = priv$rand_bootstrap_draw_matrices(draws)
	expect_null(mats_no_noise$noise_mat)

	for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = rep(b, n)
	mats_noise = priv$rand_bootstrap_draw_matrices(draws)
	expect_true(is.matrix(mats_noise$noise_mat))
	expect_equal(dim(mats_noise$noise_mat), c(n, 4L))
	expect_equal(mats_noise$noise_mat[, 3], rep(3, n))

	# all-or-nothing contract: one draw missing smooth_noise -> NULL, like w_b today
	draws2 = draws
	draws2[[2]][["smooth_noise"]] = NULL
	expect_null(priv$rand_bootstrap_draw_matrices(draws2)$noise_mat)
})

test_that("Ridit and Jonckheere-Terpstra fast kernels decline smoothed (noise) draws", {
	set.seed(20260722)
	n = 12
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "ordinal")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(sample(1:4, n, replace = TRUE))

	inf_ridit = InferenceOrdinalRidit$new(des)
	priv_ridit = inf_ridit$.__enclos_env__$private
	priv_ridit$shared()
	draws = priv_ridit$generate_rand_bootstrap_draws(B = 3L, materialize_w = TRUE)
	for (b in seq_along(draws)) draws[[b]][["smooth_noise"]] = rnorm(n)
	expect_null(priv_ridit$compute_fast_rand_bootstrap_distr(as.numeric(priv_ridit$y), draws, 0, "none"))

	inf_jt = InferenceOrdinalJonckheereTerpstraTest$new(des)
	priv_jt = inf_jt$.__enclos_env__$private
	draws_jt = priv_jt$generate_rand_bootstrap_draws(B = 3L, materialize_w = TRUE)
	for (b in seq_along(draws_jt)) draws_jt[[b]][["smooth_noise"]] = rnorm(n)
	expect_null(priv_jt$compute_fast_rand_bootstrap_distr(as.numeric(priv_jt$y), draws_jt, 0, "none"))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-noise-mat-plumbing.R")'`
Expected: FAIL — `rand_bootstrap_draw_matrices` does not return a `noise_mat` element yet (first `expect_true(is.matrix(mats_noise$noise_mat))` fails since the whole `noise_mat` name doesn't exist), and the Ridit/JT `expect_null(...)` calls fail because today the fast kernels don't know about `smooth_noise` at all and would try (and fail differently, e.g. wrong-arity or non-NULL) — the important thing is the test does not pass as-is.

- [ ] **Step 3: Implement `rand_bootstrap_draw_matrices`**

In `EDI/R/inference_all_abstract_rand_bootstrap.R`, replace:

```r
		rand_bootstrap_draw_matrices = function(rand_bootstrap_draws){
			n = as.integer(private$n)
			B = length(rand_bootstrap_draws)
			if (B == 0L) return(NULL)
			i_mat = matrix(NA_integer_, nrow = n, ncol = B)
			w_mat = matrix(NA_integer_, nrow = n, ncol = B)
			for (b in seq_len(B)) {
				draw = rand_bootstrap_draws[[b]]
				if (is.null(draw$w_b) || length(draw$i_b) != n || length(draw$w_b) != n) return(NULL)
				i_mat[, b] = as.integer(draw$i_b)
				w_mat[, b] = as.integer(draw$w_b)
			}
			list(i_mat = i_mat, w_mat = w_mat)
		},
```

with:

```r
		rand_bootstrap_draw_matrices = function(rand_bootstrap_draws){
			n = as.integer(private$n)
			B = length(rand_bootstrap_draws)
			if (B == 0L) return(NULL)
			i_mat = matrix(NA_integer_, nrow = n, ncol = B)
			w_mat = matrix(NA_integer_, nrow = n, ncol = B)
			has_noise = !is.null(rand_bootstrap_draws[[1L]][["smooth_noise"]])
			noise_mat = if (has_noise) matrix(0, nrow = n, ncol = B) else NULL
			for (b in seq_len(B)) {
				draw = rand_bootstrap_draws[[b]]
				if (is.null(draw$w_b) || length(draw$i_b) != n || length(draw$w_b) != n) return(NULL)
				i_mat[, b] = as.integer(draw$i_b)
				w_mat[, b] = as.integer(draw$w_b)
				if (has_noise) {
					if (is.null(draw[["smooth_noise"]]) || length(draw[["smooth_noise"]]) != n) return(NULL)
					noise_mat[, b] = as.numeric(draw[["smooth_noise"]])
				}
			}
			list(i_mat = i_mat, w_mat = w_mat, noise_mat = noise_mat)
		},
```

- [ ] **Step 4: Remove the blanket fast-kernel bypass at site A**

In `EDI/R/inference_all_abstract_rand_bootstrap.R`, inside `approximate_rand_bootstrap_distribution_beta_hat_T`, replace:

```r
			has_custom_randomization_statistic =
				!is.null(private[["custom_randomization_statistic_function"]]) ||
				!is.null(private[["compiled_cpp_stat_fn"]])
			has_fast_kernel = !has_custom_randomization_statistic &&
				private$has_private_method("compute_fast_rand_bootstrap_distr")
			# Smoothed draws carry per-draw noise the C++ kernel cannot see; force worker/iteration path.
			if (has_fast_kernel && draws_supplied && length(rand_bootstrap_draws) > 0L && !is.null(rand_bootstrap_draws[[1L]][["smooth_noise"]]))
				has_fast_kernel = FALSE
```

with:

```r
			has_custom_randomization_statistic =
				!is.null(private[["custom_randomization_statistic_function"]]) ||
				!is.null(private[["compiled_cpp_stat_fn"]])
			has_fast_kernel = !has_custom_randomization_statistic &&
				private$has_private_method("compute_fast_rand_bootstrap_distr")
```

(Capability for noise-carrying draws is now decided per-class inside `compute_fast_rand_bootstrap_distr` itself — see Tasks 2-9 and the Ridit/JT change below.)

- [ ] **Step 5: Remove the blanket fast-kernel bypass at site B**

In the same file, inside `get_brt_distribution_prefix`, replace:

```r
			has_custom = !is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])
			has_fast_kernel = !has_custom && private$has_private_method("compute_fast_rand_bootstrap_distr")
			if (has_fast_kernel && length(draws) > 0L && !is.null(draws[[1L]][["smooth_noise"]]))
				has_fast_kernel = FALSE
```

with:

```r
			has_custom = !is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])
			has_fast_kernel = !has_custom && private$has_private_method("compute_fast_rand_bootstrap_distr")
```

- [ ] **Step 6: Make Ridit explicitly decline noise-carrying draws**

In `EDI/R/inference_ordinal_ridit.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# ordinal: no sharp-null shift supported
			if (delta != 0) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_ridit_rand_bootstrap_parallel_cpp(
				as.integer(y0_full), mats$i_mat, mats$w_mat, private$reference,
				private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# ordinal: no sharp-null shift supported
			if (delta != 0) return(NULL)
			# "smoothed" adds continuous Gaussian noise, which is not meaningful for integer
			# ordinal category codes; decline and let the R-level fallback (which truncates
			# via as.integer()) handle it, unchanged from before this kernel existed.
			if (length(rand_bootstrap_draws) > 0L && !is.null(rand_bootstrap_draws[[1L]][["smooth_noise"]])) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_ridit_rand_bootstrap_parallel_cpp(
				as.integer(y0_full), mats$i_mat, mats$w_mat, private$reference,
				private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 7: Make Jonckheere-Terpstra explicitly decline noise-carrying draws**

In `EDI/R/inference_ordinal_jonckheere_terpstra_test.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# ordinal: no sharp-null shift supported
			if (delta != 0) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_jt_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), mats$i_mat, mats$w_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# ordinal: no sharp-null shift supported
			if (delta != 0) return(NULL)
			# "smoothed" adds continuous Gaussian noise, which is not meaningful for integer
			# ordinal category codes; decline and let the R-level fallback (which truncates
			# via as.integer()) handle it, unchanged from before this kernel existed.
			if (length(rand_bootstrap_draws) > 0L && !is.null(rand_bootstrap_draws[[1L]][["smooth_noise"]])) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_jt_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), mats$i_mat, mats$w_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 8: Run the test again to verify it passes**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-noise-mat-plumbing.R")'`
Expected: PASS, 2 tests, 0 failures. (No C++ recompilation needed for this task — pure R changes.)

- [ ] **Step 9: Commit**

```bash
git add EDI/R/inference_all_abstract_rand_bootstrap.R EDI/R/inference_ordinal_ridit.R EDI/R/inference_ordinal_jonckheere_terpstra_test.R EDI/tests/testthat/test-brt-smoothed-noise-mat-plumbing.R
git commit -m "BRT smoothed: thread noise_mat through rand_bootstrap_draw_matrices; move fast-kernel capability check into each class"
```

---

### Task 2: Wilcox kernel — `compute_wilcox_hl_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_wilcox_hl.cpp:539-590`
- Modify: `EDI/R/inference_all_simple_wilcox.R:180-190` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-wilcox-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1's `rand_bootstrap_draw_matrices`.
- Produces: `compute_wilcox_hl_rand_bootstrap_parallel_cpp(y0_sexp, i_mat_sexp, w_mat_sexp, delta, transform_code, zero_one_logit_clamp, noise_mat, num_cores)` — same as today plus one new `noise_mat` argument (a `Rcpp::Nullable<Rcpp::NumericMatrix>`) inserted before `num_cores`.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-wilcox-kernel.R`. This uses the *identity-resampling equivalence* trick: with `i_mat` set to the identity permutation in every column, adding `noise_mat[, b]` inside the kernel must give exactly the same result as pre-adding that column's noise to `y0` and calling the kernel with `noise_mat = NULL`.

```r
library(EDI)

test_that("compute_wilcox_hl_rand_bootstrap_parallel_cpp adds noise before the sharp-null shift", {
	set.seed(20260723)
	n = 40; B = 5
	y0 = rnorm(n, sd = 3)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	# guarantee both arms present in every column
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.7), nrow = n, ncol = B)
	delta = 0.3
	transform_code = 0L
	clamp = .Machine$double.eps

	res_noisy = EDI:::compute_wilcox_hl_rand_bootstrap_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, noise_mat, 1L
	)

	res_ref = numeric(B)
	for (b in seq_len(B)) {
		y0_b = y0 + noise_mat[, b]
		res_ref[b] = EDI:::compute_wilcox_hl_rand_bootstrap_parallel_cpp(
			as.numeric(y0_b), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, transform_code, clamp, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-10)

	# noise_mat = NULL must be a no-op relative to an explicit zero matrix
	res_null = EDI:::compute_wilcox_hl_rand_bootstrap_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, NULL, 1L
	)
	res_zero = EDI:::compute_wilcox_hl_rand_bootstrap_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, matrix(0, n, B), 1L
	)
	expect_equal(res_null, res_zero, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-wilcox-kernel.R")'`
Expected: FAIL — `unused argument` / wrong number of arguments, since the kernel does not accept `noise_mat` yet.

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_wilcox_hl.cpp`, replace:

```cpp
NumericVector compute_wilcox_hl_rand_bootstrap_parallel_cpp(
    SEXP y0_sexp,
    SEXP i_mat_sexp,
    SEXP w_mat_sexp,
    double delta,
    int transform_code,
    double zero_one_logit_clamp,
    int num_cores) {

	NumericVector y0_vec(y0_sexp);
	IntegerMatrix i_int_mat(i_mat_sexp);
	IntegerMatrix w_int_mat(w_mat_sexp);
    Eigen::Map<const Eigen::VectorXd> y0(y0_vec.begin(), y0_vec.size());
    Eigen::Map<const Eigen::MatrixXi> i_mat(i_int_mat.begin(), i_int_mat.nrow(), i_int_mat.ncol());
    Eigen::Map<const Eigen::MatrixXi> w_mat(w_int_mat.begin(), w_int_mat.nrow(), w_int_mat.ncol());

    const int n = i_mat.rows();
    const int nsim = i_mat.cols();
    std::vector<double> results_vec(nsim, NA_REAL);
    const double* y0_ptr = y0.data();
    const int* i_ptr = i_mat.data();
    const int* w_ptr = w_mat.data();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        const int* i_col = i_ptr + (size_t)b * n;
        const int* w_col = w_ptr + (size_t)b * n;
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
            if (!std::isfinite(yv)) continue;
            if (w_col[i] == 1) {
                y_t.push_back(delta != 0.0 ? apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv);
            } else if (w_col[i] == 0) {
                y_c.push_back(yv);
            }
        }

        res_ptr[b] = hl_from_groups(std::move(y_t), std::move(y_c));
    }

    return wrap(results_vec);
}
```

with:

```cpp
NumericVector compute_wilcox_hl_rand_bootstrap_parallel_cpp(
    SEXP y0_sexp,
    SEXP i_mat_sexp,
    SEXP w_mat_sexp,
    double delta,
    int transform_code,
    double zero_one_logit_clamp,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores) {

	NumericVector y0_vec(y0_sexp);
	IntegerMatrix i_int_mat(i_mat_sexp);
	IntegerMatrix w_int_mat(w_mat_sexp);
    Eigen::Map<const Eigen::VectorXd> y0(y0_vec.begin(), y0_vec.size());
    Eigen::Map<const Eigen::MatrixXi> i_mat(i_int_mat.begin(), i_int_mat.nrow(), i_int_mat.ncol());
    Eigen::Map<const Eigen::MatrixXi> w_mat(w_int_mat.begin(), w_int_mat.nrow(), w_int_mat.ncol());

    const int n = i_mat.rows();
    const int nsim = i_mat.cols();
    std::vector<double> results_vec(nsim, NA_REAL);
    const double* y0_ptr = y0.data();
    const int* i_ptr = i_mat.data();
    const int* w_ptr = w_mat.data();
    double* res_ptr = results_vec.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        const int* i_col = i_ptr + (size_t)b * n;
        const int* w_col = w_ptr + (size_t)b * n;
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
            if (has_noise) yv += noise_ptr[(size_t)b * n + i];
            if (!std::isfinite(yv)) continue;
            if (w_col[i] == 1) {
                y_t.push_back(delta != 0.0 ? apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv);
            } else if (w_col[i] == 0) {
                y_c.push_back(yv);
            }
        }

        res_ptr[b] = hl_from_groups(std::move(y_t), std::move(y_c));
    }

    return wrap(results_vec);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_all_simple_wilcox.R`, replace:

```r
			compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
				if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
				transform_code = private$rand_bootstrap_transform_code(transform_responses)
				if (is.null(transform_code)) return(NULL)
				mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
				if (is.null(mats)) return(NULL)
				compute_wilcox_hl_rand_bootstrap_parallel_cpp(
					as.numeric(y0_full), mats$i_mat, mats$w_mat, as.numeric(delta),
					transform_code, as.numeric(zero_one_logit_clamp), private$n_cpp_threads(ncol(mats$w_mat))
				)
			},
```

with:

```r
			compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
				if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
				transform_code = private$rand_bootstrap_transform_code(transform_responses)
				if (is.null(transform_code)) return(NULL)
				mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
				if (is.null(mats)) return(NULL)
				compute_wilcox_hl_rand_bootstrap_parallel_cpp(
					as.numeric(y0_full), mats$i_mat, mats$w_mat, as.numeric(delta),
					transform_code, as.numeric(zero_one_logit_clamp), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
				)
			},
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

Run:
```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-wilcox-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_wilcox_hl.cpp EDI/R/inference_all_simple_wilcox.R EDI/tests/testthat/test-brt-smoothed-wilcox-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for Wilcox HL rand-bootstrap"
```

---

### Task 3: Mean-difference kernel — `compute_rand_bootstrap_mean_diff_parallel_cpp`

**Files:**
- Modify: `EDI/src/rand_bootstrap_mean_diff_parallel.cpp:53-118` (7-arg exported function + its 5-arg internal-only overload)
- Modify: `EDI/R/inference_all_mean_diff.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-mean-diff-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_rand_bootstrap_mean_diff_parallel_cpp(y0, i_mat, w_mat, delta, transform_code, zero_one_logit_clamp, noise_mat, num_cores)`.
- Note: this file also has a **plain C++ overload** (not `Rcpp::export`ed, unused from R — verified via `grep -rn "compute_rand_bootstrap_mean_diff_parallel_cpp" EDI/R/*.R EDI/src/*.cpp`, only the 7-arg exported version is called) that forwards to the 7-arg version with fixed `transform_code = 0, zero_one_logit_clamp = NA_REAL`. It must be updated to compile (pass `R_NilValue` for the new argument) even though nothing calls it today.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-mean-diff-kernel.R`:

```r
library(EDI)

test_that("compute_rand_bootstrap_mean_diff_parallel_cpp adds noise before the sharp-null shift", {
	set.seed(20260724)
	n = 30; B = 5
	y0 = rnorm(n, sd = 2)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.5), nrow = n, ncol = B)
	delta = -0.2
	transform_code = 0L
	clamp = .Machine$double.eps

	res_noisy = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, transform_code, clamp, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-10)

	res_null = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, NULL, 1L
	)
	res_zero = EDI:::compute_rand_bootstrap_mean_diff_parallel_cpp(
		as.numeric(y0), i_mat, w_mat, delta, transform_code, clamp, matrix(0, n, B), 1L
	)
	expect_equal(res_null, res_zero, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-mean-diff-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel and fix the internal overload**

In `EDI/src/rand_bootstrap_mean_diff_parallel.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int transform_code,
	double zero_one_logit_clamp,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = w_mat.rows();

	std::vector<double> results_vec(nsim);

	const double* y0_ptr = y0.begin();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;
		double sum_T = 0, sum_C = 0;
		int n_T = 0;

		for (int i = 0; i < n; ++i) {
			const double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
			const int is_t = (w_col[i] == 1);
			if (is_t) {
				sum_T += (delta != 0.0) ? brt_apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv;
				++n_T;
			} else {
				sum_C += yv;
			}
		}

		const int n_C = n - n_T;
		if (n_T == 0 || n_C == 0) {
			res_ptr[b] = NA_REAL;
		} else {
			res_ptr[b] = (sum_T / n_T) - (sum_C / n_C);
		}
	}

	return wrap(results_vec);
}

// Compatibility entry point for the generated Rcpp wrapper.  The R API still
// uses the original five-argument additive-scale interface.
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {
	return compute_rand_bootstrap_mean_diff_parallel_cpp(
		y0, i_mat, w_mat, delta, 0, NA_REAL, num_cores
	);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int transform_code,
	double zero_one_logit_clamp,
	Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = w_mat.rows();

	std::vector<double> results_vec(nsim);

	const double* y0_ptr = y0.begin();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

	const bool has_noise = noise_mat.isNotNull();
	NumericMatrix noise_m;
	const double* noise_ptr = nullptr;
	if (has_noise) {
		noise_m = NumericMatrix(noise_mat);
		noise_ptr = noise_m.begin();
	}

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;
		double sum_T = 0, sum_C = 0;
		int n_T = 0;

		for (int i = 0; i < n; ++i) {
			double yv = y0_ptr[i_col[i] - 1]; // i_mat is 1-based
			if (has_noise) yv += noise_ptr[(size_t)b * n + i];
			const int is_t = (w_col[i] == 1);
			if (is_t) {
				sum_T += (delta != 0.0) ? brt_apply_shift(yv, delta, transform_code, zero_one_logit_clamp) : yv;
				++n_T;
			} else {
				sum_C += yv;
			}
		}

		const int n_C = n - n_T;
		if (n_T == 0 || n_C == 0) {
			res_ptr[b] = NA_REAL;
		} else {
			res_ptr[b] = (sum_T / n_T) - (sum_C / n_C);
		}
	}

	return wrap(results_vec);
}

// Compatibility entry point for the generated Rcpp wrapper.  The R API still
// uses the original five-argument additive-scale interface.
NumericVector compute_rand_bootstrap_mean_diff_parallel_cpp(
	const NumericVector& y0,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {
	return compute_rand_bootstrap_mean_diff_parallel_cpp(
		y0, i_mat, w_mat, delta, 0, NA_REAL, R_NilValue, num_cores
	);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_all_mean_diff.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			transform_code = private$rand_bootstrap_transform_code(transform_responses)
			if (is.null(transform_code)) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_rand_bootstrap_mean_diff_parallel_cpp(
				as.numeric(y0_full), mats$i_mat, mats$w_mat, as.numeric(delta),
				transform_code, as.numeric(zero_one_logit_clamp), private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps) {
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			transform_code = private$rand_bootstrap_transform_code(transform_responses)
			if (is.null(transform_code)) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_rand_bootstrap_mean_diff_parallel_cpp(
				as.numeric(y0_full), mats$i_mat, mats$w_mat, as.numeric(delta),
				transform_code, as.numeric(zero_one_logit_clamp), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-mean-diff-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/rand_bootstrap_mean_diff_parallel.cpp EDI/R/inference_all_mean_diff.R EDI/tests/testthat/test-brt-smoothed-mean-diff-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for simple mean-difference rand-bootstrap"
```

---

### Task 4: OLS kernel — `compute_rand_bootstrap_ols_parallel_cpp`

**Files:**
- Modify: `EDI/src/rand_bootstrap_ols_parallel.cpp:22-82`
- Modify: `EDI/R/inference_continuous_ols.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-ols-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_rand_bootstrap_ols_parallel_cpp(y0, Xc, i_mat, w_mat, delta, noise_mat, num_cores)`.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-ols-kernel.R`:

```r
library(EDI)

test_that("compute_rand_bootstrap_ols_parallel_cpp adds noise before the sharp-null shift", {
	set.seed(20260725)
	n = 40; B = 5; p_cov = 2
	y0 = rnorm(n, sd = 2)
	Xc = matrix(rnorm(n * p_cov), nrow = n, ncol = p_cov)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.6), nrow = n, ncol = B)
	delta = 0.4

	res_noisy = EDI:::compute_rand_bootstrap_ols_parallel_cpp(
		as.numeric(y0), Xc, i_mat, w_mat, delta, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_rand_bootstrap_ols_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), Xc, i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-8)

	res_null = EDI:::compute_rand_bootstrap_ols_parallel_cpp(as.numeric(y0), Xc, i_mat, w_mat, delta, NULL, 1L)
	res_zero = EDI:::compute_rand_bootstrap_ols_parallel_cpp(as.numeric(y0), Xc, i_mat, w_mat, delta, matrix(0, n, B), 1L)
	expect_equal(res_null, res_zero, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-ols-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/rand_bootstrap_ols_parallel.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_ols_parallel_cpp(
	const NumericVector& y0,
	const NumericMatrix& Xc,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {

	const int nsim = w_mat.cols();
	const int n = w_mat.rows();
	const int p_cov = Xc.ncol();
	const int p = 2 + p_cov;

	std::vector<double> results_vec(nsim, NA_REAL);

	const double* y0_ptr = y0.begin();
	const double* xc_ptr = (p_cov > 0) ? Xc.begin() : nullptr;
	const int n_full = y0.size();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n * p * p, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;

		Eigen::MatrixXd M(n, p);
		Eigen::VectorXd yb(n);
		int n_T = 0;
		for (int i = 0; i < n; ++i) {
			const int row0 = i_col[i] - 1; // i_mat is 1-based
			const int is_t = (w_col[i] == 1);
			M(i, 0) = 1.0;
			M(i, 1) = static_cast<double>(is_t);
			for (int j = 0; j < p_cov; ++j) {
				M(i, 2 + j) = xc_ptr[(size_t)j * n_full + row0];
			}
			yb(i) = y0_ptr[row0] + (is_t ? delta : 0.0);
			n_T += is_t;
		}
		if (n_T == 0 || n_T == n || n <= p) continue; // leaves NA_REAL

		const Eigen::MatrixXd MtM = M.transpose() * M;
		const Eigen::VectorXd Mty = M.transpose() * yb;
		Eigen::LDLT<Eigen::MatrixXd> ldlt(MtM);
		if (ldlt.info() != Eigen::Success) continue;
		const Eigen::VectorXd beta = ldlt.solve(Mty);
		// guard against a numerically singular design (LDLT does not fail on rank deficiency)
		const double resid_norm = (MtM * beta - Mty).norm();
		if (!std::isfinite(beta(1)) || resid_norm > 1e-6 * (1.0 + Mty.norm())) continue;
		res_ptr[b] = beta(1);
	}

	return wrap(results_vec);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_rand_bootstrap_ols_parallel_cpp(
	const NumericVector& y0,
	const NumericMatrix& Xc,
	const IntegerMatrix& i_mat,
	const IntegerMatrix& w_mat,
	double delta,
	Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
	int num_cores) {

	const int nsim = w_mat.cols();
	const int n = w_mat.rows();
	const int p_cov = Xc.ncol();
	const int p = 2 + p_cov;

	std::vector<double> results_vec(nsim, NA_REAL);

	const double* y0_ptr = y0.begin();
	const double* xc_ptr = (p_cov > 0) ? Xc.begin() : nullptr;
	const int n_full = y0.size();
	const int* i_ptr = i_mat.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n * p * p, num_cores);

	const bool has_noise = noise_mat.isNotNull();
	NumericMatrix noise_m;
	const double* noise_ptr = nullptr;
	if (has_noise) {
		noise_m = NumericMatrix(noise_mat);
		noise_ptr = noise_m.begin();
	}

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* i_col = i_ptr + (size_t)b * n;
		const int* w_col = w_ptr + (size_t)b * n;

		Eigen::MatrixXd M(n, p);
		Eigen::VectorXd yb(n);
		int n_T = 0;
		for (int i = 0; i < n; ++i) {
			const int row0 = i_col[i] - 1; // i_mat is 1-based
			const int is_t = (w_col[i] == 1);
			M(i, 0) = 1.0;
			M(i, 1) = static_cast<double>(is_t);
			for (int j = 0; j < p_cov; ++j) {
				M(i, 2 + j) = xc_ptr[(size_t)j * n_full + row0];
			}
			double yv = y0_ptr[row0];
			if (has_noise) yv += noise_ptr[(size_t)b * n + i];
			yb(i) = yv + (is_t ? delta : 0.0);
			n_T += is_t;
		}
		if (n_T == 0 || n_T == n || n <= p) continue; // leaves NA_REAL

		const Eigen::MatrixXd MtM = M.transpose() * M;
		const Eigen::VectorXd Mty = M.transpose() * yb;
		Eigen::LDLT<Eigen::MatrixXd> ldlt(MtM);
		if (ldlt.info() != Eigen::Success) continue;
		const Eigen::VectorXd beta = ldlt.solve(Mty);
		// guard against a numerically singular design (LDLT does not fail on rank deficiency)
		const double resid_norm = (MtM * beta - Mty).norm();
		if (!std::isfinite(beta(1)) || resid_norm > 1e-6 * (1.0 + Mty.norm())) continue;
		res_ptr[b] = beta(1);
	}

	return wrap(results_vec);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_continuous_ols.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# the OLS kernel only implements the additive sharp-null shift
			if (delta != 0 && !identical(transform_responses, "none")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_full = private$create_design_matrix()
			Xc = if (ncol(X_full) > 2L) {
				as.matrix(X_full[, -(1:2), drop = FALSE])
			} else {
				matrix(numeric(0), nrow = as.integer(private$n), ncol = 0L)
			}
			compute_rand_bootstrap_ols_parallel_cpp(
				as.numeric(y0_full), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# the OLS kernel only implements the additive sharp-null shift
			if (delta != 0 && !identical(transform_responses, "none")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_full = private$create_design_matrix()
			Xc = if (ncol(X_full) > 2L) {
				as.matrix(X_full[, -(1:2), drop = FALSE])
			} else {
				matrix(numeric(0), nrow = as.integer(private$n), ncol = 0L)
			}
			compute_rand_bootstrap_ols_parallel_cpp(
				as.numeric(y0_full), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-ols-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/rand_bootstrap_ols_parallel.cpp EDI/R/inference_continuous_ols.R EDI/tests/testthat/test-brt-smoothed-ols-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for OLS rand-bootstrap"
```

---

### Task 5: Robust regression kernel — `compute_robust_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_robust_regression.cpp:317-381`
- Modify: `EDI/R/inference_continuous_robust_regr.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_robust_rand_bootstrap_parallel_cpp(y0, Xc, i_mat, w_mat, delta, method, noise_mat, num_cores)`.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R`:

```r
library(EDI)

test_that("compute_robust_rand_bootstrap_parallel_cpp adds noise before the sharp-null shift", {
	set.seed(20260726)
	n = 50; B = 4; p_cov = 2
	y0 = rnorm(n, sd = 2)
	Xc = matrix(rnorm(n * p_cov), nrow = n, ncol = p_cov)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.5), nrow = n, ncol = B)
	delta = 0.25
	method = "huber"

	res_noisy = EDI:::compute_robust_rand_bootstrap_parallel_cpp(
		as.numeric(y0), Xc, i_mat, w_mat, delta, method, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_robust_rand_bootstrap_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), Xc, i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, method, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-8)

	res_null = EDI:::compute_robust_rand_bootstrap_parallel_cpp(as.numeric(y0), Xc, i_mat, w_mat, delta, method, NULL, 1L)
	res_zero = EDI:::compute_robust_rand_bootstrap_parallel_cpp(as.numeric(y0), Xc, i_mat, w_mat, delta, method, matrix(0, n, B), 1L)
	expect_equal(res_null, res_zero, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_robust_regression.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_robust_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    std::string method,
    int num_cores)
{
    const int n      = i_mat.nrow();
    const int nsim   = i_mat.ncol();
    const int n_full = y0.size();
    const int p_cov  = Xc.ncol();
    const int p      = 2 + p_cov;  // intercept + treatment + covariates

    const double* y0_ptr = y0.begin();
    const double* xc_ptr = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr  = i_mat.begin();
    const int*    w_ptr  = w_mat.begin();

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = 1.0;
                X_b(i, 1) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 2) = xc_ptr[(size_t)j * n_full + r];
                y_b(i) = y0_ptr[r] + delta * wt;
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            RobustModelResult res = fast_robust_regression_internal(
                X_b, y_b,
                R_NilValue, true, method,
                1.345, 4.685, 50, 1e-7, -1.0,
                R_NilValue, R_NilValue,
                R_NilValue, R_NilValue,
                true, 0);
            if (res.b.size() > 1 && std::isfinite(res.b[1]))
                res_ptr[b] = res.b[1];
        }
    }
    return wrap(results);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_robust_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    std::string method,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores)
{
    const int n      = i_mat.nrow();
    const int nsim   = i_mat.ncol();
    const int n_full = y0.size();
    const int p_cov  = Xc.ncol();
    const int p      = 2 + p_cov;  // intercept + treatment + covariates

    const double* y0_ptr = y0.begin();
    const double* xc_ptr = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr  = i_mat.begin();
    const int*    w_ptr  = w_mat.begin();

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = 1.0;
                X_b(i, 1) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 2) = xc_ptr[(size_t)j * n_full + r];
                double yv = y0_ptr[r];
                if (has_noise) yv += noise_ptr[(size_t)b * n + i];
                y_b(i) = yv + delta * wt;
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            RobustModelResult res = fast_robust_regression_internal(
                X_b, y_b,
                R_NilValue, true, method,
                1.345, 4.685, 50, 1e-7, -1.0,
                R_NilValue, R_NilValue,
                R_NilValue, R_NilValue,
                true, 0);
            if (res.b.size() > 1 && std::isfinite(res.b[1]))
                res_ptr[b] = res.b[1];
        }
    }
    return wrap(results);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_continuous_robust_regr.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_data = private$get_X()
			Xc = if (!is.null(private$best_X_colnames) && length(private$best_X_colnames) > 0L) {
				cols = intersect(private$best_X_colnames, colnames(X_data))
				if (length(cols) > 0L) as.matrix(X_data[, cols, drop = FALSE]) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			} else if (!is.null(X_data) && ncol(as.matrix(X_data)) > 0L) {
				as.matrix(X_data)
			} else {
				matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_robust_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$rlm_method, private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_data = private$get_X()
			Xc = if (!is.null(private$best_X_colnames) && length(private$best_X_colnames) > 0L) {
				cols = intersect(private$best_X_colnames, colnames(X_data))
				if (length(cols) > 0L) as.matrix(X_data[, cols, drop = FALSE]) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			} else if (!is.null(X_data) && ncol(as.matrix(X_data)) > 0L) {
				as.matrix(X_data)
			} else {
				matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_robust_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$rlm_method, mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_robust_regression.cpp EDI/R/inference_continuous_robust_regr.R EDI/tests/testthat/test-brt-smoothed-robust-regr-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for robust regression rand-bootstrap"
```

---

### Task 6: CoxPH kernel — `compute_coxph_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_coxph_regression.cpp:819-885`
- Modify: `EDI/R/inference_survival_coxph.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-coxph-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_coxph_rand_bootstrap_parallel_cpp(y0, dead, Xc, i_mat, w_mat, delta, noise_mat, num_cores)`.
- Note: leave the sibling function `compute_coxph_rand_bootstrap_cpp` (same file, lines 776-817) untouched — it is exported but not called from any R private method (verified via grep) and is out of scope.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-coxph-kernel.R`:

```r
library(EDI)

test_that("compute_coxph_rand_bootstrap_parallel_cpp adds noise before the multiplicative sharp-null shift", {
	set.seed(20260727)
	n = 60; B = 4; p_cov = 2
	y0 = rexp(n, rate = 0.2) + 0.5 # strictly positive survival times
	dead = rbinom(n, 1, 0.8)
	Xc = matrix(rnorm(n * p_cov), nrow = n, ncol = p_cov)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	# small noise relative to the time scale so times stay well-behaved
	noise_mat = matrix(rnorm(n * B, sd = 0.1), nrow = n, ncol = B)
	delta = 0.15

	res_noisy = EDI:::compute_coxph_rand_bootstrap_parallel_cpp(
		as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_coxph_rand_bootstrap_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), as.integer(dead), Xc, i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-8)

	res_null = EDI:::compute_coxph_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, NULL, 1L)
	res_zero = EDI:::compute_coxph_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, matrix(0, n, B), 1L)
	expect_equal(res_null, res_zero, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-coxph-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_coxph_regression.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_coxph_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    int num_cores)
{
    const int n      = i_mat.nrow();
    const int nsim   = i_mat.ncol();
    const int n_full = y0.size();
    const int p_cov  = Xc.ncol();
    const int p      = 1 + p_cov;  // treatment + covariates; no intercept in Cox PH

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(p);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 1) = xc_ptr[(size_t)j * n_full + r];
                y_b(i)    = (wt && mult != 1.0) ? y0_ptr[r] * mult : y0_ptr[r];
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            std::vector<CoxData> strata;
            strata.emplace_back(y_b, dead_b, X_b);
            CoxFitResult fit = cox_newton_raphson(strata, R_NilValue, true, fspec,
                                                  true, 20, 1e-6);
            if (fit.converged && !fit.beta.empty() && std::isfinite(fit.beta[0]))
                res_ptr[b] = fit.beta[0];
        }
    }
    return wrap(results);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_coxph_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores)
{
    const int n      = i_mat.nrow();
    const int nsim   = i_mat.ncol();
    const int n_full = y0.size();
    const int p_cov  = Xc.ncol();
    const int p      = 1 + p_cov;  // treatment + covariates; no intercept in Cox PH

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(p);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 1) = xc_ptr[(size_t)j * n_full + r];
                double yv = y0_ptr[r];
                if (has_noise) yv += noise_ptr[(size_t)b * n + i];
                y_b(i)    = (wt && mult != 1.0) ? yv * mult : yv;
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            std::vector<CoxData> strata;
            strata.emplace_back(y_b, dead_b, X_b);
            CoxFitResult fit = cox_newton_raphson(strata, R_NilValue, true, fspec,
                                                  true, 20, 1e-6);
            if (fit.converged && !fit.beta.empty() && std::isfinite(fit.beta[0]))
                res_ptr[b] = fit.beta[0];
        }
    }
    return wrap(results);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_survival_coxph.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			Xc = if (!is.null(private$cox_X_fit_cache) && ncol(private$cox_X_fit_cache) > 1L) {
				as.matrix(private$cox_X_fit_cache[, -1L, drop = FALSE])
			} else {
				X_cov = private$get_X()
				if (!is.null(X_cov) && ncol(as.matrix(X_cov)) > 0L) as.matrix(X_cov) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_coxph_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			Xc = if (!is.null(private$cox_X_fit_cache) && ncol(private$cox_X_fit_cache) > 1L) {
				as.matrix(private$cox_X_fit_cache[, -1L, drop = FALSE])
			} else {
				X_cov = private$get_X()
				if (!is.null(X_cov) && ncol(as.matrix(X_cov)) > 0L) as.matrix(X_cov) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_coxph_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-coxph-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_coxph_regression.cpp EDI/R/inference_survival_coxph.R EDI/tests/testthat/test-brt-smoothed-coxph-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for CoxPH rand-bootstrap"
```

---

### Task 7: Weibull marginal kernel — `compute_weibull_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_weibull_regression.cpp:137-203`
- Modify: `EDI/R/inference_survival_KK_weibull_marginal.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-weibull-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_weibull_rand_bootstrap_parallel_cpp(y0, dead, Xc, i_mat, w_mat, delta, noise_mat, num_cores)`.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-weibull-kernel.R`:

```r
library(EDI)

test_that("compute_weibull_rand_bootstrap_parallel_cpp adds noise before the multiplicative sharp-null shift", {
	set.seed(20260728)
	n = 60; B = 4; p_cov = 2
	y0 = rexp(n, rate = 0.2) + 0.5
	dead = rbinom(n, 1, 0.8)
	Xc = matrix(rnorm(n * p_cov), nrow = n, ncol = p_cov)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.1), nrow = n, ncol = B)
	delta = 0.2

	res_noisy = EDI:::compute_weibull_rand_bootstrap_parallel_cpp(
		as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_weibull_rand_bootstrap_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), as.integer(dead), Xc, i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-8)

	res_null = EDI:::compute_weibull_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, NULL, 1L)
	res_zero = EDI:::compute_weibull_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), Xc, i_mat, w_mat, delta, matrix(0, n, B), 1L)
	expect_equal(res_null, res_zero, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-weibull-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_weibull_regression.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_weibull_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    int num_cores)
{
    const int n       = i_mat.nrow();
    const int nsim    = i_mat.ncol();
    const int n_full  = y0.size();
    const int p_cov   = Xc.ncol();
    const int p       = 2 + p_cov;   // intercept + treatment + covariates
    const int n_params = p + 1;       // +1 for log_sigma

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(n_params);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n), params0(n_params);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = 1.0;
                X_b(i, 1) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 2) = xc_ptr[(size_t)j * n_full + r];
                y_b(i)    = (wt && mult != 1.0) ? y0_ptr[r] * mult : y0_ptr[r];
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            params0.setZero();
            WeibullAFTLikelihood fun(y_b, dead_b, X_b);
            LikelihoodFitResult fit = optimize_fixed_likelihood(
                fun, params0, fspec, 100, 1e-8, "lbfgs", "", 0, nullptr);
            if (fit.converged && fit.params.size() > 1 && std::isfinite(fit.params[1]))
                res_ptr[b] = fit.params[1];
        }
    }
    return wrap(results);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_weibull_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const NumericMatrix& Xc,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores)
{
    const int n       = i_mat.nrow();
    const int nsim    = i_mat.ncol();
    const int n_full  = y0.size();
    const int p_cov   = Xc.ncol();
    const int p       = 2 + p_cov;   // intercept + treatment + covariates
    const int n_params = p + 1;       // +1 for log_sigma

    const double* y0_ptr   = y0.begin();
    const int*    dead_ptr = dead.begin();
    const double* xc_ptr   = (p_cov > 0) ? Xc.begin() : nullptr;
    const int*    i_ptr    = i_mat.begin();
    const int*    w_ptr    = w_mat.begin();
    const double  mult     = (delta != 0.0) ? std::exp(delta) : 1.0;

    FixedParamSpec fspec = make_fixed_param_spec(n_params);

    std::vector<double> results(nsim, NA_REAL);
    double* res_ptr = results.data();

    const bool has_noise = noise_mat.isNotNull();
    NumericMatrix noise_m;
    const double* noise_ptr = nullptr;
    if (has_noise) {
        noise_m = NumericMatrix(noise_mat);
        noise_ptr = noise_m.begin();
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
    {
        Eigen::VectorXd y_b(n), dead_b(n), params0(n_params);
        Eigen::MatrixXd X_b(n, p);

#pragma omp for schedule(dynamic)
        for (int b = 0; b < nsim; ++b) {
            const int* ic = i_ptr + (size_t)b * n;
            const int* wc = w_ptr + (size_t)b * n;

            int n_t = 0, n_c = 0;
            for (int i = 0; i < n; ++i) {
                const int r  = ic[i] - 1;
                const int wt = (wc[i] == 1) ? 1 : 0;
                X_b(i, 0) = 1.0;
                X_b(i, 1) = static_cast<double>(wt);
                for (int j = 0; j < p_cov; ++j)
                    X_b(i, j + 2) = xc_ptr[(size_t)j * n_full + r];
                double yv = y0_ptr[r];
                if (has_noise) yv += noise_ptr[(size_t)b * n + i];
                y_b(i)    = (wt && mult != 1.0) ? yv * mult : yv;
                dead_b(i) = static_cast<double>(dead_ptr[r]);
                n_t += wt;
                n_c += (1 - wt);
            }
            if (n_t < 2 || n_c < 2) continue;

            params0.setZero();
            WeibullAFTLikelihood fun(y_b, dead_b, X_b);
            LikelihoodFitResult fit = optimize_fixed_likelihood(
                fun, params0, fspec, 100, 1e-8, "lbfgs", "", 0, nullptr);
            if (fit.converged && fit.params.size() > 1 && std::isfinite(fit.params[1]))
                res_ptr[b] = fit.params[1];
        }
    }
    return wrap(results);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_survival_KK_weibull_marginal.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_data = private$get_X()
			Xc = if (!is.null(private$best_X_colnames)) {
				cols = intersect(private$best_X_colnames, colnames(X_data))
				if (length(cols) > 0L) as.matrix(X_data[, cols, drop = FALSE]) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			} else if (!is.null(X_data) && ncol(as.matrix(X_data)) > 0L) {
				as.matrix(X_data)
			} else {
				matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_weibull_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			X_data = private$get_X()
			Xc = if (!is.null(private$best_X_colnames)) {
				cols = intersect(private$best_X_colnames, colnames(X_data))
				if (length(cols) > 0L) as.matrix(X_data[, cols, drop = FALSE]) else matrix(numeric(0), nrow = private$n, ncol = 0L)
			} else if (!is.null(X_data) && ncol(as.matrix(X_data)) > 0L) {
				as.matrix(X_data)
			} else {
				matrix(numeric(0), nrow = private$n, ncol = 0L)
			}
			compute_weibull_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), Xc, mats$i_mat, mats$w_mat,
				as.numeric(delta), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		}
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-weibull-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_weibull_regression.cpp EDI/R/inference_survival_KK_weibull_marginal.R EDI/tests/testthat/test-brt-smoothed-weibull-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for Weibull marginal rand-bootstrap"
```

---

### Task 8: Log-rank kernel — `compute_logrank_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_logrank.cpp:148-193`
- Modify: `EDI/R/inference_survival_log_rank.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-logrank-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_logrank_rand_bootstrap_parallel_cpp(y0, dead, i_mat, w_mat, delta, noise_mat, num_cores)`.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-logrank-kernel.R`:

```r
library(EDI)

test_that("compute_logrank_rand_bootstrap_parallel_cpp adds noise before the multiplicative sharp-null shift", {
	set.seed(20260729)
	n = 60; B = 4
	y0 = rexp(n, rate = 0.2) + 0.5
	dead = rbinom(n, 1, 0.8)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.1), nrow = n, ncol = B)
	delta = 0.2

	res_noisy = EDI:::compute_logrank_rand_bootstrap_parallel_cpp(
		as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, noise_mat, 1L
	)
	res_ref = numeric(B)
	for (b in seq_len(B)) {
		res_ref[b] = EDI:::compute_logrank_rand_bootstrap_parallel_cpp(
			as.numeric(y0 + noise_mat[, b]), as.integer(dead), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
			delta, NULL, 1L
		)
	}
	expect_equal(res_noisy, res_ref, tolerance = 1e-8)

	res_null = EDI:::compute_logrank_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, NULL, 1L)
	res_zero = EDI:::compute_logrank_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, matrix(0, n, B), 1L)
	expect_equal(res_null, res_zero, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-logrank-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_logrank.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_logrank_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const double mult = std::exp(delta);
  const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
  if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
  {
    Eigen::VectorXd time_b(n);
    std::vector<int> dead_b(n), w_b(n);

#pragma omp for schedule(static)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      for (int i = 0; i < n; ++i) {
        const int row0 = i_col[i] - 1; // i_mat is 1-based
        const int is_t = (w_col[i] == 1);
        time_b[i] = (is_t && delta != 0.0) ? y0_ptr[row0] * mult : y0_ptr[row0];
        dead_b[i] = dead_ptr[row0];
        w_b[i] = is_t;
      }
      ModelResult r = fast_logrank_internal(time_b, dead_b, w_b);
      res_ptr[b] = (r.b.size() > 0) ? r.b[0] : NA_REAL;
    }
  }

  return wrap(results_vec);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_logrank_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const double mult = std::exp(delta);
  const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

  const bool has_noise = noise_mat.isNotNull();
  NumericMatrix noise_m;
  const double* noise_ptr = nullptr;
  if (has_noise) {
    noise_m = NumericMatrix(noise_mat);
    noise_ptr = noise_m.begin();
  }

#ifdef _OPENMP
  if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
  {
    Eigen::VectorXd time_b(n);
    std::vector<int> dead_b(n), w_b(n);

#pragma omp for schedule(static)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      for (int i = 0; i < n; ++i) {
        const int row0 = i_col[i] - 1; // i_mat is 1-based
        const int is_t = (w_col[i] == 1);
        double yv = y0_ptr[row0];
        if (has_noise) yv += noise_ptr[(size_t)b * n + i];
        time_b[i] = (is_t && delta != 0.0) ? yv * mult : yv;
        dead_b[i] = dead_ptr[row0];
        w_b[i] = is_t;
      }
      ModelResult r = fast_logrank_internal(time_b, dead_b, w_b);
      res_ptr[b] = (r.b.size() > 0) ? r.b[0] : NA_REAL;
    }
  }

  return wrap(results_vec);
}
```

- [ ] **Step 4: Update the R call site**

In `EDI/R/inference_survival_log_rank.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# survival sharp-null shift is multiplicative (delta on the log scale)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_logrank_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			# survival sharp-null shift is multiplicative (delta on the log scale)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_logrank_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-logrank-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_logrank.cpp EDI/R/inference_survival_log_rank.R EDI/tests/testthat/test-brt-smoothed-logrank-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for log-rank rand-bootstrap"
```

---

### Task 9: Survival stat-diff kernel (RMST + KM-diff) — `compute_survival_stat_diff_rand_bootstrap_parallel_cpp`

**Files:**
- Modify: `EDI/src/fast_survival_stats.cpp:343-399`
- Modify: `EDI/R/inference_survival_rmst.R` (`compute_fast_rand_bootstrap_distr`)
- Modify: `EDI/R/inference_survival_km_diff.R` (`compute_fast_rand_bootstrap_distr`)
- Test: `EDI/tests/testthat/test-brt-smoothed-survival-stat-diff-kernel.R`

**Interfaces:**
- Consumes: `mats$noise_mat` from Task 1.
- Produces: `compute_survival_stat_diff_rand_bootstrap_parallel_cpp(y0, dead, i_mat, w_mat, delta, do_rmst, noise_mat, num_cores)`. One kernel serves two R classes (`InferenceSurvivalRestrictedMeanDiff` with `do_rmst = TRUE`, `InferenceSurvivalKMDiff` with `do_rmst = FALSE`); both call sites get the identical `mats$noise_mat` addition.

- [ ] **Step 1: Write the failing test**

Create `EDI/tests/testthat/test-brt-smoothed-survival-stat-diff-kernel.R`:

```r
library(EDI)

test_that("compute_survival_stat_diff_rand_bootstrap_parallel_cpp adds noise before the multiplicative sharp-null shift (both do_rmst values)", {
	set.seed(20260730)
	n = 60; B = 4
	y0 = rexp(n, rate = 0.2) + 0.5
	dead = rbinom(n, 1, 0.8)
	i_mat = matrix(rep(1:n, B), nrow = n, ncol = B)
	w_mat = matrix(rbinom(n * B, 1, 0.5), nrow = n, ncol = B)
	for (b in seq_len(B)) { w_mat[1, b] = 1L; w_mat[2, b] = 0L }
	noise_mat = matrix(rnorm(n * B, sd = 0.1), nrow = n, ncol = B)
	delta = 0.2

	for (do_rmst in c(TRUE, FALSE)) {
		res_noisy = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
			as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, noise_mat, 1L
		)
		res_ref = numeric(B)
		for (b in seq_len(B)) {
			res_ref[b] = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0 + noise_mat[, b]), as.integer(dead), i_mat[, b, drop = FALSE], w_mat[, b, drop = FALSE],
				delta, do_rmst, NULL, 1L
			)
		}
		expect_equal(res_noisy, res_ref, tolerance = 1e-8, info = paste("do_rmst =", do_rmst))

		res_null = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, NULL, 1L)
		res_zero = EDI:::compute_survival_stat_diff_rand_bootstrap_parallel_cpp(as.numeric(y0), as.integer(dead), i_mat, w_mat, delta, do_rmst, matrix(0, n, B), 1L)
		expect_equal(res_null, res_zero, tolerance = 1e-10, info = paste("do_rmst =", do_rmst))
	}
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-survival-stat-diff-kernel.R")'`
Expected: FAIL (wrong number of arguments).

- [ ] **Step 3: Add `noise_mat` to the C++ kernel**

In `EDI/src/fast_survival_stats.cpp`, replace:

```cpp
// [[Rcpp::export]]
NumericVector compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    bool do_rmst,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const double mult = std::exp(delta);

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
  {
    // Per-thread reusable buffers: avoids heap allocation in the hot loop
    std::vector<SurvEntry> y_t(n), y_c(n);
    std::vector<double> utimes_t, utimes_c, sprobs_t, sprobs_c;
    utimes_t.reserve(n); sprobs_t.reserve(n);
    utimes_c.reserve(n); sprobs_c.reserve(n);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      int nt = 0, nc = 0;
      for (int i = 0; i < n; ++i) {
        const int row0 = i_col[i] - 1;  // i_mat is 1-based
        SurvEntry e;
        e.status = dead_ptr[row0];
        if (w_col[i] == 1) {
          e.time = (delta != 0.0) ? y0_ptr[row0] * mult : y0_ptr[row0];
          y_t[nt++] = e;
        } else {
          e.time = y0_ptr[row0];
          y_c[nc++] = e;
        }
      }
      if (nt == 0 || nc == 0) continue;
      double stat_t = km_stat_inline(y_t.data(), nt, do_rmst, utimes_t, sprobs_t);
      double stat_c = km_stat_inline(y_c.data(), nc, do_rmst, utimes_c, sprobs_c);
      if (std::isfinite(stat_t) && std::isfinite(stat_c))
        res_ptr[b] = stat_t - stat_c;
    }
  }
  return wrap(results_vec);
}
```

with:

```cpp
// [[Rcpp::export]]
NumericVector compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
    const NumericVector& y0,
    const IntegerVector& dead,
    const IntegerMatrix& i_mat,
    const IntegerMatrix& w_mat,
    double delta,
    bool do_rmst,
    Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat,
    int num_cores) {

  const int n = i_mat.nrow();
  const int nsim = i_mat.ncol();
  std::vector<double> results_vec(nsim, NA_REAL);
  const double* y0_ptr = y0.begin();
  const int* dead_ptr = dead.begin();
  const int* i_ptr = i_mat.begin();
  const int* w_ptr = w_mat.begin();
  double* res_ptr = results_vec.data();
  const double mult = std::exp(delta);

  const bool has_noise = noise_mat.isNotNull();
  NumericMatrix noise_m;
  const double* noise_ptr = nullptr;
  if (has_noise) {
    noise_m = NumericMatrix(noise_mat);
    noise_ptr = noise_m.begin();
  }

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(num_cores > 1)
  {
    // Per-thread reusable buffers: avoids heap allocation in the hot loop
    std::vector<SurvEntry> y_t(n), y_c(n);
    std::vector<double> utimes_t, utimes_c, sprobs_t, sprobs_c;
    utimes_t.reserve(n); sprobs_t.reserve(n);
    utimes_c.reserve(n); sprobs_c.reserve(n);

#pragma omp for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
      const int* i_col = i_ptr + (size_t)b * n;
      const int* w_col = w_ptr + (size_t)b * n;
      int nt = 0, nc = 0;
      for (int i = 0; i < n; ++i) {
        const int row0 = i_col[i] - 1;  // i_mat is 1-based
        double yv = y0_ptr[row0];
        if (has_noise) yv += noise_ptr[(size_t)b * n + i];
        SurvEntry e;
        e.status = dead_ptr[row0];
        if (w_col[i] == 1) {
          e.time = (delta != 0.0) ? yv * mult : yv;
          y_t[nt++] = e;
        } else {
          e.time = yv;
          y_c[nc++] = e;
        }
      }
      if (nt == 0 || nc == 0) continue;
      double stat_t = km_stat_inline(y_t.data(), nt, do_rmst, utimes_t, sprobs_t);
      double stat_c = km_stat_inline(y_c.data(), nc, do_rmst, utimes_c, sprobs_c);
      if (std::isfinite(stat_t) && std::isfinite(stat_c))
        res_ptr[b] = stat_t - stat_c;
    }
  }
  return wrap(results_vec);
}
```

- [ ] **Step 4: Update both R call sites**

In `EDI/R/inference_survival_rmst.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), TRUE, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), TRUE, mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

In `EDI/R/inference_survival_km_diff.R`, replace:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), FALSE, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

with:

```r
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]]) || !is.null(private[["compiled_cpp_stat_fn"]])) return(NULL)
			if (delta != 0 && !identical(transform_responses, "log")) return(NULL)
			mats = private$rand_bootstrap_draw_matrices(rand_bootstrap_draws)
			if (is.null(mats)) return(NULL)
			compute_survival_stat_diff_rand_bootstrap_parallel_cpp(
				as.numeric(y0_full), as.integer(private$dead), mats$i_mat, mats$w_mat,
				as.numeric(delta), FALSE, mats$noise_mat, private$n_cpp_threads(ncol(mats$w_mat))
			)
		},
```

- [ ] **Step 5: Regenerate bindings, install, and run the test**

```bash
Rscript -e 'Rcpp::compileAttributes("EDI")'
R CMD INSTALL --no-docs EDI
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-survival-stat-diff-kernel.R")'
```
Expected: PASS, 1 test, 0 failures.

- [ ] **Step 6: Commit**

```bash
git add EDI/src/fast_survival_stats.cpp EDI/R/inference_survival_rmst.R EDI/R/inference_survival_km_diff.R EDI/tests/testthat/test-brt-smoothed-survival-stat-diff-kernel.R EDI/src/RcppExports.cpp EDI/R/RcppExports.R
git commit -m "BRT smoothed: noise-aware fast kernel for RMST/KM-diff rand-bootstrap"
```

---

### Task 10: Full regression pass — existing BRT test suite

**Files:**
- None modified — verification only.

**Interfaces:**
- Consumes: the fully-updated package from Tasks 1-9.

- [ ] **Step 1: Run the pre-existing BRT-adjacent test files**

```bash
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-rand-bootstrap.R")'
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-wilcox-regr-bootstrap-fast-path.R")'
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-bootstrap-reused-worker-families.R")'
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-bootstrap-reused-worker-asymp-families.R")'
```
Expected: all PASS, 0 failures — these exercise `compute_fast_rand_bootstrap_distr` and the reused-worker path for `percentile`/`studentized` types, none of which changed behavior (the `has_noise` branches are no-ops when `noise_mat` is `NULL`, which is what every non-smoothed caller still passes).

- [ ] **Step 2: Run the 9 new kernel test files together**

```bash
Rscript -e 'library(EDI); testthat::test_dir("EDI/tests/testthat", filter = "brt-smoothed")'
```
Expected: all PASS (Task 1's plumbing test + the 8 kernel tests from Tasks 2-9 — the RMST/KM-diff kernel is one file covering both `do_rmst` values, so 9 kernel behaviors across 8 files plus the foundation file).

- [ ] **Step 3: If anything fails, stop and diagnose before continuing**

Do not proceed to Task 11 with a red suite. If a specific kernel test fails, re-open that kernel's task and re-check the noise-addition line placement (it must occur before any shift/multiply, for every row regardless of `w`).

- [ ] **Step 4: No commit for this task (verification-only)**

---

### Task 11: End-to-end correctness + performance verification (Wilcox)

**Files:**
- Test: `EDI/tests/testthat/test-brt-smoothed-wilcox-ci-perf.R`

**Interfaces:**
- Consumes: `InferenceAllSimpleWilcox$compute_rand_bootstrap_confidence_interval(type = "smoothed")` (public API, unchanged signature).

This is the test that directly validates the original problem is fixed: it reproduces the profiling setup from the investigation (n=30, B=99, single smoothed p-value took 3.5s before this plan) at the CI level, and asserts both correctness (fast-path output matches the old slow-fallback output for identical seeds) and speed (CI now completes quickly).

- [ ] **Step 1: Write the test**

Create `EDI/tests/testthat/test-brt-smoothed-wilcox-ci-perf.R`:

```r
library(EDI)

SlowInferenceAllSimpleWilcox = R6::R6Class(
	"SlowInferenceAllSimpleWilcox",
	inherit = InferenceAllSimpleWilcox,
	private = list(
		compute_fast_rand_bootstrap_distr = function(y0_full, rand_bootstrap_draws, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			NULL
		}
	)
)

test_that("smoothed CI: fast-kernel result matches the forced-slow-fallback result", {
	set.seed(20260731)
	n = 30
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))

	fast_inf = InferenceAllSimpleWilcox$new(des)
	slow_inf = SlowInferenceAllSimpleWilcox$new(des)

	set.seed(99)
	fast_ci = fast_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)
	set.seed(99)
	slow_ci = slow_inf$compute_rand_bootstrap_confidence_interval(B = 99, type = "smoothed", show_progress = FALSE)

	expect_equal(as.numeric(fast_ci), as.numeric(slow_ci), tolerance = 1e-6)
})

test_that("smoothed p-value: fast kernel is dramatically faster than the forced-slow fallback", {
	set.seed(20260732)
	n = 30
	des = DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
	X = data.frame(x1 = rnorm(n))
	for (i in seq_len(n)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
	des$add_all_subject_responses(rnorm(n))

	fast_inf = InferenceAllSimpleWilcox$new(des)
	slow_inf = SlowInferenceAllSimpleWilcox$new(des)

	set.seed(7)
	t_fast = system.time(
		fast_pv <- fast_inf$compute_rand_bootstrap_two_sided_pval(B = 99, type = "smoothed", show_progress = FALSE)
	)[["elapsed"]]
	set.seed(7)
	t_slow = system.time(
		slow_pv <- slow_inf$compute_rand_bootstrap_two_sided_pval(B = 99, type = "smoothed", show_progress = FALSE)
	)[["elapsed"]]

	expect_equal(fast_pv, slow_pv, tolerance = 1e-6)
	# The investigation measured ~7x for a single call; require at least 3x here to allow for
	# machine variance while still catching a regression back to the fully-slow path.
	expect_true(t_fast < t_slow / 3, info = sprintf("t_fast=%.3fs t_slow=%.3fs", t_fast, t_slow))
})
```

- [ ] **Step 2: Run the test**

```bash
Rscript -e 'library(EDI); testthat::test_file("EDI/tests/testthat/test-brt-smoothed-wilcox-ci-perf.R")'
```
Expected: PASS, 2 tests, 0 failures. If the speed assertion is flaky in CI (shared/loaded machines), that's a signal to loosen the ratio, not to delete the test — do not weaken correctness assertions to fix timing flakiness.

- [ ] **Step 3: Commit**

```bash
git add EDI/tests/testthat/test-brt-smoothed-wilcox-ci-perf.R
git commit -m "BRT smoothed: end-to-end correctness + speedup regression test for Wilcox CI"
```

---

### Task 12: Documentation — roxygen update, perf log, memory

**Files:**
- Modify: `EDI/R/inference_all_abstract_rand_bootstrap.R` (roxygen for `compute_rand_bootstrap_two_sided_pval`, the "Implementation caveat" paragraph)
- Modify: `EDI/R/inference_all_abstract_rand_bootstrap_ci.R` (roxygen pointer paragraph)
- Modify: `package_metadata/perf_experiments_final.md` (new dated entry)
- Modify (memory, outside the repo): `project_smoothed_brt_fast_path_bypass.md` in the auto-memory directory

This task runs **after** Task 11 so the roxygen and perf-log text can cite real measured numbers instead of placeholders.

- [ ] **Step 1: Measure the real speedup to cite**

```bash
Rscript -e '
suppressMessages(library(EDI))
set.seed(1); n <- 30
des <- DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous")
X <- data.frame(x1 = rnorm(n))
for (i in 1:n) des$add_one_subject_to_experiment_and_assign(X[i,,drop=FALSE])
des$add_all_subject_responses(rnorm(n))
inf <- InferenceAllSimpleWilcox$new(des)
t <- system.time(pv <- inf$compute_rand_bootstrap_two_sided_pval(B = 99, type = "smoothed", show_progress = FALSE))
cat("post-fix smoothed pval elapsed:", t[["elapsed"]], "s\n")
'
```
Record the printed elapsed time (pre-fix was 3.512s on this same n=30, B=99 workload, per the original investigation) to use in Step 2 and Step 3 below.

- [ ] **Step 2: Update the roxygen "Implementation caveat" paragraph**

In `EDI/R/inference_all_abstract_rand_bootstrap.R`, replace the paragraph:

```r
		#'   \strong{Implementation caveat (ad hoc, not a certified instance of the above).} The
		#'   bandwidth used here, \eqn{\hat{\sigma}/\sqrt{n}} (the raw SE-of-the-mean scale of the
		#'   response), is a pragmatic engineering choice, not one derived from or validated against
		#'   the bandwidth-selection results in the sources above. Those references smooth the
		#'   \emph{resampling distribution itself} with a bandwidth chosen to trade off bias against
		#'   variance (often shrinking slower than \eqn{n^{-1/2}}, e.g. \eqn{n^{-1/5}}-type
		#'   KDE rates); here, noise is instead added directly to each already-resampled response at
		#'   a fixed \eqn{n^{-1/2}} rate. The bandwidth is not exposed as a parameter, has no
		#'   zero-noise escape hatch (the only way to disable smoothing is to pick a different
		#'   \code{type}), and its coverage behavior has not been validated by simulation in this
		#'   package. Treat it as a discreteness patch that is qualitatively motivated by the
		#'   literature above, not a certified implementation of it.
```

with (fill in the measured elapsed time from Step 1 where marked):

```r
		#'   \strong{Implementation caveat (ad hoc, not a certified instance of the above).} The
		#'   bandwidth used here, \eqn{\hat{\sigma}/\sqrt{n}} (the raw SE-of-the-mean scale of the
		#'   response), is a pragmatic engineering choice, not one derived from or validated against
		#'   the bandwidth-selection results in the sources above. Those references smooth the
		#'   \emph{resampling distribution itself} with a bandwidth chosen to trade off bias against
		#'   variance (often shrinking slower than \eqn{n^{-1/2}}, e.g. \eqn{n^{-1/5}}-type
		#'   KDE rates); here, noise is instead added directly to each already-resampled response at
		#'   a fixed \eqn{n^{-1/2}} rate. The bandwidth is not exposed as a parameter, has no
		#'   zero-noise escape hatch (the only way to disable smoothing is to pick a different
		#'   \code{type}), and its coverage behavior has not been validated by simulation in this
		#'   package. Treat it as a discreteness patch that is qualitatively motivated by the
		#'   literature above, not a certified implementation of it.
		#'
		#'   \strong{Performance.} Every class with a C++ \code{compute_fast_rand_bootstrap_distr}
		#'   kernel that operates on real-valued responses (Wilcox HL, simple mean difference, OLS,
		#'   robust regression, CoxPH, Weibull marginal, log-rank, RMST, KM-diff) accepts the smoothing
		#'   noise directly in its kernel, so \code{"smoothed"} is as fast as \code{"percentile"} for
		#'   those classes (measured: a single \code{compute_rand_bootstrap_two_sided_pval(type =
		#'   "smoothed")} call on \code{n = 30}, \code{B = 99} for \code{InferenceAllSimpleWilcox} went
		#'   from 3.512s to <MEASURED_SECONDS>s). The two ordinal classes (\code{InferenceOrdinalRidit},
		#'   \code{InferenceOrdinalJonckheereTerpstraTest}) still use the slower R-level fallback for
		#'   \code{"smoothed"}, because adding continuous Gaussian noise to integer category codes is
		#'   not statistically meaningful — see the response-type caveat above.
```

- [ ] **Step 3: Add a dated entry to `package_metadata/perf_experiments_final.md`**

Append (matching the file's existing entry style — see the entries already in the file for tone/format):

```
## BRT smoothed fast-kernel noise support — 2026-07-20

`compute_rand_bootstrap_two_sided_pval`/`compute_rand_bootstrap_confidence_interval(type = "smoothed")` unconditionally disabled every class's C++ `compute_fast_rand_bootstrap_distr` kernel whenever a draw carried per-draw `smooth_noise`, because none of the kernels could consume it; every smoothed evaluation fell back to the R-level reused-worker/per-iteration path, which for `InferenceAllSimpleWilcox` calls `stats::wilcox.test()` once per bootstrap draw. Extended the 8 C++ kernels that operate on real-valued responses (`compute_wilcox_hl_rand_bootstrap_parallel_cpp`, `compute_rand_bootstrap_mean_diff_parallel_cpp`, `compute_rand_bootstrap_ols_parallel_cpp`, `compute_robust_rand_bootstrap_parallel_cpp`, `compute_coxph_rand_bootstrap_parallel_cpp`, `compute_weibull_rand_bootstrap_parallel_cpp`, `compute_logrank_rand_bootstrap_parallel_cpp`, `compute_survival_stat_diff_rand_bootstrap_parallel_cpp`) with an optional `Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat` argument, added before the sharp-null shift exactly where the R fallback already adds it (`y_sim = y0[i_b] + smooth_noise`, applied to every row regardless of treatment status). `rand_bootstrap_draw_matrices()` now packs this matrix whenever draws carry `smooth_noise`; the two blanket dispatcher checks that disabled the fast path (`approximate_rand_bootstrap_distribution_beta_hat_T` and `get_brt_distribution_prefix` in `inference_all_abstract_rand_bootstrap.R`) were removed in favor of each class deciding for itself. The two ordinal classes (`InferenceOrdinalRidit`, `InferenceOrdinalJonckheereTerpstraTest`) explicitly decline noise-carrying draws in their own `compute_fast_rand_bootstrap_distr`, since continuous noise on integer category codes isn't meaningful, and are unchanged. Correctness: 8 new per-kernel tests use an identity-resampling equivalence check (`kernel(y0, noise_mat=N, ...)` on identity `i_mat` columns must equal `kernel(y0 + N[,b], noise_mat=NULL, ...)` per column) plus a `noise_mat=NULL` vs. explicit zero-matrix no-op check; a further end-to-end test compares the fast-kernel CI/pval output against a forced-slow-fallback subclass for identical seeds (tolerance 1e-6). Measured: `InferenceAllSimpleWilcox$compute_rand_bootstrap_two_sided_pval(B=99, type="smoothed")` on `n=30` went from 3.512s to <MEASURED_SECONDS>s (single call, matching the investigation that motivated this change). Installed with `R CMD INSTALL --no-docs`.
```

(Replace `<MEASURED_SECONDS>` with the value from Step 1 in both this entry and the roxygen paragraph in Step 2.)

- [ ] **Step 4: Regenerate NAMESPACE/docs**

```bash
Rscript fast_roxygenize.R
```
Expected: completes; the only warnings should be the pre-existing unrelated ones (undocumented `supports_lik_ratio_param_bootstrap`, incidence CI param docs) already present before this plan — do not attempt to fix those here (out of scope).

- [ ] **Step 5: Update the auto-memory file**

The memory file `project_smoothed_brt_fast_path_bypass.md` (in the auto-memory directory referenced by this session, not in the repo) currently ends with: *"As of 2026-07-20 this has not been logged as a TODO in `package_metadata/perf_experiments_final.md`."* Update that closing line to state the fix has shipped, referencing this plan and the measured speedup, and add a one-line index entry update in `MEMORY.md` reflecting the resolved state (not a new file — edit the existing entry's description).

- [ ] **Step 6: Final full-suite sanity check**

```bash
Rscript -e 'library(EDI); testthat::test_dir("EDI/tests/testthat", filter = "brt-smoothed|rand-bootstrap|wilcox-regr-bootstrap")'
```
Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add EDI/R/inference_all_abstract_rand_bootstrap.R EDI/R/inference_all_abstract_rand_bootstrap_ci.R EDI/man EDI/NAMESPACE package_metadata/perf_experiments_final.md
git commit -m "BRT smoothed: document fast-kernel noise support, log perf results"
```

---

## Self-Review

**Spec coverage:** "extend each `compute_fast_rand_bootstrap_distr` C++ kernel to accept an optional per-draw noise vector" — Tasks 2-9 cover all 8 distinct C++ export functions behind the 9 classes with real-valued responses (RMST + KM-diff share one kernel/file, covered together in Task 9). The 2 ordinal kernels are explicitly and separately handled (declined, not silently ignored) in Task 1, with the rationale documented inline and in Task 12's roxygen update. The shared dispatcher plumbing (`rand_bootstrap_draw_matrices`, the two blanket bypass sites) is Task 1. Correctness is TDD'd per kernel (Tasks 2-9) plus a full-suite regression pass (Task 10) plus an end-to-end + perf regression test tied to the original motivating case (Task 11). Documentation/perf-log/memory updates close the loop (Task 12).

**Placeholder scan:** no TBD/"add validation"/"similar to Task N" — every code step shows complete before/after code. The one intentional placeholder, `<MEASURED_SECONDS>`, is explicitly called out as "fill in from Step 1's measured output" in both places it appears, not left as an implementation gap.

**Type consistency:** `noise_mat` is `Rcpp::Nullable<Rcpp::NumericMatrix>` in every kernel, always the parameter immediately before `num_cores`, always passed explicitly (no default) from R as `mats$noise_mat`. `rand_bootstrap_draw_matrices()`'s return list gains exactly one new element, `noise_mat`, consistently named across all 9 consuming call sites.
