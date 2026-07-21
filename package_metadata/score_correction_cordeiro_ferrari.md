# TODO: Cordeiro & Ferrari Modified Score Test

## Scope

This is an implementation plan for adding Cordeiro & Ferrari's (1991, *Biometrika*)
modified score test — the score-test analog of the Bartlett-corrected likelihood-ratio
test already implemented for `lik_ratio_bartlett_approx`/`lik_ratio_bartlett_exact` — to
the package's likelihood-backed inference paths (`InferenceAsympLik` and its
descendants), plus a confidence interval obtained by inverting the corrected test.

It mirrors the existing Bartlett work in `package_metadata/likrat_correction_bartlett.md`
and reuses as much of its plumbing as possible. Read that report first; this document
assumes familiarity with the `approx`/`exact` split, the `InferenceParamBootstrap`
`simulate_under_lik_null()` machinery, and the generic `invert_test_pval_confidence_interval()`
root-finder.

## Background (short version)

Rao's classical score test only requires the null-restricted fit — never the
unrestricted/alternative fit. Like Wald and LR, the classical score statistic `S` is only
χ²_q to *first order* (error `O(n⁻¹)`). Cordeiro & Ferrari's correction is a **polynomial**
adjustment, not a simple rescaling like Bartlett's:

```
S* = S · (1 − (c + b·S + a·S²))
```

where `a`, `b`, `c` are built from third/fourth derivatives of the GLM mean/variance
functions contracted against the leverages of the full and null-restricted models. After
correction, `S*` is χ²_q to `O(n⁻²)`.

**Key finding from the Bartlett effort that changes the risk profile here**: CRAN package
[`mdscore`](https://cran.r-project.org/package=mdscore) (Silva & Silva Jr., GPL, plain
readable R source, `mdscore/R/mdscore.r`) already implements exactly this correction for
single-equation GLMs (gaussian/Gamma/inverse.gaussian/poisson/binomial families;
log/logit/probit/cloglog/cauchit/identity links). I hand-verified its derivative
quantities (`f`, `g`, `b`, `h`) against my own derivation of `b'''(θ)`/`b''''(θ)` for
logistic regression — they reduce algebraically to the same expressions (`g_i = 0`
identically for canonical link is a real identity, not a coincidence; `b_i = h_i =
b''''(θ_i)`). This gives us something the Bartlett-LR effort never had: **a live,
numerically-checkable reference implementation** for the Easy-tier GLM cases, not just a
formula transcribed from a paywalled paper. This meaningfully de-risks the Easy tier
specifically; it does not change the Borderline/Difficult tiers, which fail for
structural reasons (partial likelihood, quadrature, mixtures, multi-threshold links) that
have nothing to do with formula-transcription risk.

## Proposed API surface (mirrors Bartlett exactly)

- Testing types: `"score_modified_approx"`, `"score_modified_exact"` (aliases:
  `score_modified`, `modified_score`, etc. — finalize during Phase 1)
- Public methods on `InferenceAsympLik`:
  - `compute_score_modified_approx_two_sided_pval(delta = 0, B = 99)`
  - `compute_score_modified_approx_confidence_interval(alpha = 0.05, B = 99)`
  - `compute_score_modified_exact_two_sided_pval(delta = 0)`
  - `compute_score_modified_exact_confidence_interval(alpha = 0.05)`
  - `compute_score_modified_two_sided_pval(delta = 0, B = 99)` / `compute_score_modified_confidence_interval(alpha = 0.05, B = 99)`
    — "best available" smart wrapper (exact wins over approx, errors if neither),
    exactly like `compute_lik_ratio_bartlett_two_sided_pval()`. Same `B`-ignored-warning
    behavior when `B` is explicitly supplied but the exact path is used.
- Private hooks:
  - `supports_score_modified_ratio_approx()` / `supports_score_modified_ratio_exact()`
  - `get_score_modified_correction_approx(spec, delta, full_fit, null_fit, B = 99)`
  - `get_score_modified_correction_exact(spec, delta, full_fit, null_fit)`
  - Open question (Phase 1): should these hooks return the corrected **statistic**
    directly, or the **`(a, b, c)` triple** so the base class applies `S*` uniformly?
    Leaning toward `(a, b, c)` for symmetry with Bartlett's single-`factor` return, and
    because it keeps the "how do I contract a,b,c into S*" logic in one place
    (`get_memoized_likelihood_test_pval`) rather than duplicated per family.

## Implementation phases

### Phase 0: Naming/architecture confirmation
- [ ] Confirm testing-type names and aliases with user (mirrors the `AskUserQuestion` that
      settled `lik_ratio_bartlett_approx`/`_exact`)
- [ ] Decide hook return shape: `(a, b, c)` triple vs. corrected statistic directly

### Phase 1: Shared base-class plumbing (`InferenceAsympLik`)
- [ ] Add the two new testing types to `normalize_testing_type`, `set_testing_type` docs,
      `compute_asymp_two_sided_pval`/`compute_asymp_confidence_interval` switches,
      `get_supported_testing_types_with_bartlett()` → generalize/rename to also cover
      score-modified (or add a parallel `get_supported_testing_types_with_corrections()`
      wrapper folding in *both* Bartlett and score-modified appending logic, since the
      "many files hard-override `get_supported_testing_types_impl()`" problem discovered
      during Bartlett work applies identically here)
- [ ] Add `get_memoized_likelihood_test_pval()` branches for `score_modified_approx` /
      `score_modified_exact`: reuse the **existing raw score statistic** computation
      (`score_test_from_score_information_cpp`, already computed for the plain `"score"`
      testing type) as `S`, then apply `S* = S·(1-(c+bS+aS²))` using the family's
      `(a,b,c)` from the appropriate hook, then `pchisq(S*, df=1, lower.tail=FALSE)`
- [ ] Thread `B` through exactly like Bartlett (`bartlett_B` param pattern →
      `score_modified_B` or similar; **watch for the same shadowing bug** found during
      Bartlett work — `InferenceAsympLikStdModCache` and `InferenceContinKKOLSOneLik`
      both override `compute_likelihood_test_two_sided_pval()` directly and will need the
      same parameter threaded through again)
- [ ] CI: extend `invert_test_pval_confidence_interval()`'s eligible-testing-types set,
      same as Bartlett — no new root-finding logic needed, it's fully generic already
- [ ] Base-class default hooks return `FALSE`/`NULL` (unimplemented), matching Bartlett's
      pattern exactly
- [ ] Smart wrapper (`compute_score_modified_two_sided_pval`/`_confidence_interval`):
      copy the exact-wins-over-approx-errors-if-neither logic and `B`-ignored warning
      from the Bartlett smart wrapper

### Phase 2: Approx (Monte-Carlo) path — cheaper than Bartlett's approx
- [ ] Implement `get_score_modified_correction_approx()` in `InferenceParamBootstrap`
      as a **bootstrap-calibrated score test**, not a literal Cordeiro-Ferrari polynomial
      fit: simulate `B` datasets under H0 at `delta` (reuse `simulate_under_lik_null()`),
      compute the **null-refit-only** score statistic `S_b` for each replicate (reuse
      `boot_spec$fit_null(delta)` and the existing score/information machinery), and use
      the empirical tail probability of the observed `S` against `{S_b}`, exactly
      analogous to `compute_lik_ratio_bootstrap_two_sided_pval()`'s empirical LR
      calibration
- [ ] **Efficiency note**: unlike Bartlett-approx (which needs both `fit_null()` *and*
      `full_fit` — an unrestricted refit — per replicate, since it needs the LR
      statistic), the score approach only ever needs the null-restricted refit. Each
      replicate should be roughly half the cost of a Bartlett-approx replicate. Don't
      reuse `compute_param_bootstrap_lr_impl()` as-is (it always computes both fits) —
      write a leaner `compute_param_bootstrap_score_impl()` that skips the unrestricted
      refit entirely
- [ ] `supports_score_modified_ratio_approx()` delegates to
      `supports_lik_ratio_param_bootstrap()` by default, same auto-opt-in-for-everyone
      pattern as Bartlett-approx, with the same carve-out mechanism for families whose
      raw score statistic is independently known to be miscalibrated (check whether the
      Zero-Inflated-Poisson/Hurdle-Poisson carve-out reason — raw LR miscalibration —
      also taints their raw score statistic; if so, carve them out here too)

### Phase 3: Exact (Cordeiro-Ferrari analytic) path — start with `InferenceIncidLogRegr`
- [ ] Adapt (with attribution to Silva & Silva Jr., `mdscore` GPL license) the
      `f`/`g`/`b`/`h` derivative quantities and `Z`/`Z2`/`Z-Z2` leverage-matrix
      construction from `mdscore/R/mdscore.r` for canonical-link logistic regression,
      where several terms simplify (`g_i = 0` identically — verified)
- [ ] Implement `get_score_modified_correction_exact()` returning `(a, b, c)` (or
      corrected statistic, per Phase 0 decision) for `InferenceIncidLogRegr`
- [ ] **Validate three ways** (much stronger than what was possible for Bartlett):
      1. Direct numerical cross-check against `mdscore::mdscore()` itself on matched
         simulated data/design — this is the big win, an actual reference
         implementation to diff against, not just plausibility
      2. Compare against our own `score_modified_approx` (MC) path — should agree within
         Monte Carlo noise, same check used for Bartlett
      3. Large-scale simulation of `E[S | H0]` and the empirical size of the corrected
         test at nominal α, across several `(n, β, X)` configurations
- [ ] Only after all three validations pass: extend to other Easy-tier canonical/simple
      links (`InferenceIncidProbitRegr`, `InferenceCountPoisson`, log-binomial,
      identity-binomial, modified-Poisson variants) — these should be nearly mechanical
      once the logistic case is validated, since `mdscore` already parameterizes over
      link/family
- [ ] `InferenceContinKKOLSOneLik` (Gaussian): flagged as **lowest risk of all** — the
      Gaussian score-modified correction (like Gaussian Bartlett) reduces to standard
      finite-sample F/t-distribution theory, not tensor algebra; consider doing this one
      *first*, even before logistic, as a warm-up with essentially zero derivation risk

### Phase 4: Broader rollout decision
- [ ] Revisit the "opt in all `InferenceParamBootstrap` families to approx" question
      (same one asked and answered "yes" for Bartlett) — likely yes again, same reasoning
- [ ] Decide whether exact rollout follows the difficulty table below strictly
      (Easy → Borderline → Difficult, stop and reassess between tiers) per the original
      Bartlett report's recommendation of a **selective rollout, not simultaneous
      all-path implementation**

### Phase 5: Tests (mirror the three Bartlett test files)
- [ ] `test-score-modified-plumbing.R` — generic dispatch, exact-vs-approx smart wrapper,
      B-ignored warning, default-off behavior (mirrors `test-bartlett-lr-plumbing.R`)
- [ ] `test-score-modified-logit.R` — real correctness tests for `InferenceIncidLogRegr`:
      factor/coefficient value tests, pval exact-formula tests, CI bracket tests, **plus
      a new test class not available for Bartlett: direct numeric comparison against
      `mdscore::mdscore()` output** (skip gracefully if `mdscore` isn't installed, same
      pattern as `skip_if_not_installed("survival")` elsewhere in this test suite)
- [ ] `test-score-modified-approx-smoke-families.R` — cross-family smoke suite (mirrors
      `test-bartlett-lr-approx-smoke-families.R`)
- [ ] Add `compute_score_modified_two_sided_pval`/`compute_score_modified_confidence_interval`
      to `comprehensive_tests.R` (new `should_run_test_family("score_modified")`, add to
      `ALL_TEST_FAMILY_FILTERS`, mirror the `supports_bartlett`/`supports_bartlett_ci`
      variable pattern)
- [ ] Update `package_tests/path_audits_source.R` with a new `score-mod-approx` /
      `score-mod-exact` column pair, same NI/NTS/✓ logic as the Bartlett columns

## Difficulty table (implementing the modified score test per path)

Same 34 likelihood/partial-likelihood paths audited in the Bartlett report. Tiers reflect
whether Cordeiro-Ferrari's GLM mean/variance-derivative + leverage machinery structurally
applies — the same regularity conditions Cordeiro's Bartlett-LR correction needs, since
both are corrections for regular, single-linear-predictor GLM likelihoods. The **Verified
against mdscore?** column is new relative to the Bartlett report and is the real
practical difference: for these specific link/family combinations, `mdscore` already
ships a checkable reference implementation.

### Incidence

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | Easy | Yes — logit is `mdscore`'s native binomial case; `g=0` identity hand-checked |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` (2-cat) | Easy | Yes — `mdscore` has an explicit probit `dmu` case |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | Easy | Yes — underlying likelihood is plain Poisson |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | Easy | Yes, same Poisson likelihood (KK sequential design needs care in leverage bookkeeping, as with Bartlett) |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | Easy | Partially — log link + binomial variance isn't a standard `glm()` family combo `mdscore` wires up automatically, but `mu=exp(η)`, `V=μ(1-μ)` are trivial to supply directly |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | Easy | Partially — identity link gives `d2μ=0` so `f=0` identically (nice simplification), but again not a wired-up `mdscore` family combo out of the box |
| `InferenceIncidKKClogitOneLik` | `fast_logistic_regression_with_var_cpp` | Borderline | No — bespoke combined conditional-logit + reservoir likelihood, not a plain single-equation GLM |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | Difficult | No — quadrature-integrated random effects has no GLM hat-matrix structure |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | Difficult | No — hybrid |

### Count

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | Easy | Yes — `mdscore`'s native log-link Poisson case |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | Easy | Yes, same underlying likelihood (robust SE is a separate concern from the score/LR path) |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | Easy | Yes, same |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | Borderline | No — extra jointly-estimated dispersion/size parameter outside `mdscore`'s fixed/simple-dispersion GLM scope |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | Difficult | No — mixture; likely inherits the same raw-statistic miscalibration that carved it out of Bartlett-approx |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | Difficult | No — same carve-out concern |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | Difficult | No — mixture + dispersion |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | Difficult | No — hurdle + dispersion |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | Difficult | No — quadrature GLMM |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | Difficult | No — hurdle + GLMM |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | Difficult | No — custom combined likelihood |

### Continuous

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy — recommended first target** | Yes — `mdscore` has a native gaussian case; also reduces to textbook finite-sample F/t theory independent of `mdscore` |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | Borderline | No — variance components / quadrature, `mdscore` doesn't extend to LMMs |

### Ordinal

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | Borderline | No — multi-threshold cumulative link; `mdscore` assumes one linear predictor, not `K-1` thresholds |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | Borderline | No — same threshold issue despite `mdscore` supporting the probit link in the single-threshold case |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | Borderline | No — same threshold issue (interesting: `mdscore` does define a cauchit link, but only for single-equation models) |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | Borderline | No — same threshold issue |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | Difficult | No — thresholds + quadrature |

### Proportion

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | Borderline | No — joint mean/precision two-parameter family outside `mdscore`'s single-parameter GLM scope, though dedicated Cordeiro/Ferrari beta-regression correction papers exist and may be independently checkable later |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | Difficult | No — three-component mixture |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | Difficult | No — quadrature GLMM |

### Survival

| Path | Engine / Fitter | Difficulty | Verified against `mdscore`? |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | Borderline | No — partial likelihood isn't a GLM mean/variance structure at all |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | Borderline | No — same, plus strata |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | Borderline | No — same partial-likelihood issue, plus combined-design bespoke-ness |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | Borderline | No — same |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | Borderline | No — AFT-style location/scale structure, closer in spirit to NegBin/Beta than to a plain GLM |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | Difficult | No — highly custom combined event/censoring likelihood |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | Difficult | No — frailty integration |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | Difficult | No — copula dependence |

**Totals**: Easy = 10, Borderline = 12, Difficult = 12 (34 paths) — identical tier counts
to the Bartlett report, since the same structural conditions gate both corrections. The
difference is entirely in the **Verified against `mdscore`** column: 8 of the 10 Easy-tier
paths now have a live reference implementation to check against (vs. zero for Bartlett-LR),
which is the actual reason to prioritize the modified score test's Easy tier over
revisiting Bartlett-exact.

## Open risks / questions to resolve before Phase 1 starts

- Does the ZIP/Hurdle-Poisson raw-LR-miscalibration carve-out (see
  `zero_augmented_model_lrt_bootstrap_disabled()`) also taint the raw *score* statistic
  for those two families? Needs checking before assuming `score_modified_approx`'s
  auto-opt-in is safe for them.
- `mdscore`'s exported `mdscore()` function takes a fitted `glm` object plus the
  nuisance-only design matrix `X1` and computes `a`/`b`/`c` in one pass; our
  architecture needs the equivalent broken into `get_score_modified_correction_exact()`
  called from the shared delta-inversion machinery at *every* candidate `delta` during CI
  search, not just once at the observed estimate. Confirm `Z`/`Z2` (which depend on `W`,
  which depends on `μ`, which depends on `delta` through the null-restricted fit) really
  do need recomputing at each `delta`, or whether some caching is possible (parallel to
  the `entry$bartlett_B` cache-invalidation trick already built for Bartlett).
- `mdscore` is GPL (≥2/3); confirm license compatibility before adapting code into this
  package, and add appropriate attribution in roxygen/comments per its license terms.
