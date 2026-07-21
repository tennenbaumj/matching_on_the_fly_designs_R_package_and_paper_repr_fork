# Bartlett-Correction Report For Likelihood-Ratio Inference Paths

## Scope

This report covers the package’s **likelihood-backed inference paths**, i.e. the
concrete classes that participate in the `get_likelihood_test_spec()` surface
used by `InferenceAsympLik` for:

- likelihood-ratio p-values
- likelihood-ratio confidence intervals obtained by inversion
- score / gradient / LR warm-start caching around constrained null fits

It excludes non-likelihood IVWC / GEE / Wald-only / pass-through paths that do
not expose a likelihood-test spec.

The question here is:

> How hard would it be to add a **Bartlett-type correction** for
> likelihood-ratio tests, and therefore for the confidence intervals obtained by
> inverting those tests, across all likelihood-backed inference paths?

I use three engineering labels:

- **Easy**: the path is regular and smooth enough that a practical Bartlett
  factor is a credible implementation target with the current architecture.
- **Borderline**: the path is smooth enough in principle, but the Bartlett
  factor would require substantial model-specific higher-order work.
- **Difficult**: the path is mixture-, quadrature-, frailty-, copula-, or
  custom-combined-likelihood enough that Bartlett support is not a realistic
  generic extension.

This is an implementation audit, not a statement about abstract asymptotic
existence.

## Short Answer

Bartlett correction is more plausible than Firth for this repo because it does
**not** change the estimator or optimizer. It only changes the calibration of
the **likelihood-ratio statistic** and therefore the LR-based confidence
interval.

That said, the correction factor is usually **model-specific**, often
**null-value-specific**, and may depend on higher-order derivatives or null
expectations. So the current base class can centralize the **plumbing**, but not
the actual mathematics.

## What Can Be Centralized

The existing `InferenceAsympLik` structure is already close to what Bartlett
correction needs. The package already has:

- `get_likelihood_test_spec()`
- null-constrained refit plumbing via `spec$fit_null(delta, start = ...)`
- cached `full_negloglik` and `null_negloglik`
- LR p-value computation through `compute_lik_ratio_two_sided_pval_impl()`
- LR CI inversion through `invert_lik_ratio_ci_newton()`

So the reusable part should live in `InferenceAsympLik`:

1. add a new testing type, e.g. `lik_ratio_bartlett`
2. add public wrappers:
   - `compute_lik_ratio_bartlett_two_sided_pval()`
   - `compute_lik_ratio_bartlett_confidence_interval()`
3. add a shared private method that:
   - computes the raw LR statistic from `full_negloglik` and `null_negloglik`
   - asks the subclass for a Bartlett factor or corrected LR statistic
   - converts the corrected statistic to a `chi-square(1)` p-value
4. reuse the same corrected p-value function when inverting the test for the CI

The subclass-specific part should be exposed through an optional hook, e.g.

```r
get_bartlett_lr_adjustment = function(spec, delta, full_fit, null_fit, full_negloglik, null_negloglik){
  NULL
}
```

or equivalently:

```r
get_bartlett_factor = function(spec, delta, full_fit, null_fit){
  NULL
}
```

with the contract that the base class computes

```text
LR_B(delta) = LR(delta) / c(delta)
```

and then uses `LR_B` for both p-values and CI inversion.

## Why The CI Story Is Straightforward Once The Test Exists

For this repo, LR confidence intervals are already obtained by inversion of the
LR test. So once a Bartlett-corrected LR p-value function exists, the CI change
is conceptually simple:

- current LR CI:
  invert `p_LR(delta)`
- Bartlett LR CI:
  invert `p_Bartlett(delta)`

No new optimizer is needed. The difficult part is only obtaining the correction
factor `c(delta)` in a reliable way for each likelihood family.

One important consequence is that the Bartlett factor may need to be recomputed
at **every null value `delta`** visited during CI inversion. So even where the
method is mathematically feasible, CI inversion can be materially more expensive
than ordinary LR inversion.

## What The Likelihood Spec Would Need

The current spec is almost enough for raw LR inference:

- `full_fit`
- `fit_null`
- `extract_start`
- `score`
- `neg_loglik`
- optionally information matrices

For analytic Bartlett correction, many families would also need either:

- a family-specific closed-form Bartlett factor function, or
- enough derivative/exponential-family structure for a model-specific helper to
  compute that factor outside the generic spec

For difficult families, a bootstrap-calibrated LR test is more realistic than a
true analytic Bartlett correction, but that would no longer be a clean generic
“Bartlett factor” implementation.

## Status Update: Approx Is Implemented, Exact Is Implemented For One Family

This report originally predates any implementation work. Since it was written:

- **Approx (Monte-Carlo) Bartlett is fully implemented** and auto-enabled for
  essentially every `InferenceParamBootstrap`-derived family (`lik_ratio_bartlett_approx`
  testing type; `supports_bartlett_likelihood_ratio_approx()` delegates to
  `supports_lik_ratio_param_bootstrap()` by default, no per-family math needed — it
  reuses the existing `simulate_under_lik_null()` refit machinery and Monte-Carlo
  estimates `c(delta) = mean(simulated LR statistics)`). The only carve-out is
  Zero-Inflated-Poisson / Hurdle-Poisson, whose raw LR statistic is independently known
  to be miscalibrated. See `test-bartlett-lr-logit.R`,
  `test-bartlett-lr-approx-smoke-families.R`.
- **Exact (closed-form analytic) Bartlett is implemented for exactly one family so
  far: `InferenceContinKKOLSOneLik`** (Gaussian OLS). This was *not* a Cordeiro-style
  tensor derivation — because that class's `neg_loglik()` holds `sigma2` fixed at the
  full model's unbiased estimate rather than re-profiling it under each null fit, the
  package's own LR statistic is algebraically identical to the classical partial
  `F(1, n-p)` statistic. The "exact factor" is the value that makes
  `LR(delta)/factor(delta)`, referred to chi-square(1), reproduce the true `F(1,n-p)`
  tail probability exactly — verified directly against base R's `lm()` classical
  t-test/F-test (not merely "plausible," bit-for-bit matched to high tolerance across
  multiple seeds and non-null `delta` values). See `test-bartlett-lr-ols-exact.R`.
  Note this deliberately does **not** match this class's Wald CI numerically — Wald
  uses an HC2 heteroskedasticity-robust sandwich SE, while the exact LR/F pivot requires
  the classical homoskedastic-errors assumption; they answer different robustness
  questions and are expected to differ.
- **Exact is still blocked for every GLM family** (including the "Easy" tier below) —
  see the practical-derivation-risk table immediately below, which reflects the actual
  attempt-and-search history (Cordeiro 1983 paywalled, ar5iv/OCR extraction unreliable,
  CRAN search found `mdscore`/`rqlm` which turned out to implement different
  corrections — see `score_correction_cordeiro_ferrari.md` — not Bartlett-LR itself).

### Practical Derivation Risk (exact path only, current as of this implementation round)

This table answers a different question than the "Audit Table" further down: not
"does closed-form Bartlett theory apply in principle" (the **Theoretical Tier** column,
identical to the original Easy/Borderline/Difficult audit), but "can I *actually*
safely reproduce and ship the exact coefficients *today*, given the sources I could
verify."

| Response Type | Path | Theoretical Tier | Practical Derivation Risk Today |
|---|---|---|---|
| Incidence | `InferenceIncidLogRegr` | Easy | Blocked — attempted; tensor coefficients unverifiable from sources I could access |
| Incidence | `InferenceIncidProbitRegr` | Easy | Blocked — same issue, non-canonical link adds another derivative layer |
| Incidence | `InferenceIncidModifiedPoisson` | Easy | Blocked |
| Incidence | `InferenceIncidKKModifiedPoisson` | Easy | Blocked |
| Incidence | `InferenceIncidLogBinomial` | Easy | Blocked |
| Incidence | `InferenceIncidBinomialIdentityRiskDiff` | Easy | Blocked |
| Incidence | `InferenceIncidKKClogitOneLik` | Borderline | Blocked (compounds with bespoke combined-likelihood risk) |
| Incidence | `InferenceIncidKKGLMM` | Difficult | Not a realistic target regardless of source access |
| Incidence | `InferenceIncidKKClogitPlusGLMMOneLik` | Difficult | Not a realistic target |
| Count | `InferenceCountPoisson` | Easy | Blocked |
| Count | `InferenceCountRobustPoisson` | Easy | Blocked |
| Count | `InferenceCountQuasiPoisson` | Easy | Blocked |
| Count | `InferenceCountNegBin` | Borderline | Blocked (dispersion parameter adds another layer) |
| Count | `InferenceCountZeroInflatedPoisson` | Difficult | Not a realistic target (also carved out of approx — raw LR miscalibrated) |
| Count | `InferenceCountHurdlePoisson` | Difficult | Not a realistic target (same carve-out) |
| Count | `InferenceCountZeroInflatedNegBin` | Difficult | Not a realistic target |
| Count | `InferenceCountHurdleNegBin` | Difficult | Not a realistic target |
| Count | `InferenceCountKKGLMM` | Difficult | Not a realistic target |
| Count | `InferenceCountKKHurdlePoissonOneLik` | Difficult | Not a realistic target |
| Count | `InferenceCountKKCPoissonOneLik` | Difficult | Not a realistic target |
| Continuous | `InferenceContinKKOLSOneLik` | Easy | **✅ Implemented** — exact finite-sample result (F/t-test identity), verified against `lm()`, not a tensor derivation at all |
| Continuous | `InferenceContinKKGLMM` | Borderline | Blocked |
| Ordinal | `InferenceOrdinalPropOddsRegr` | Borderline | Blocked (threshold parameters) |
| Ordinal | `InferenceOrdinalOrderedProbitRegr` | Borderline | Blocked |
| Ordinal | `InferenceOrdinalCauchitRegr` | Borderline | Blocked |
| Ordinal | `InferenceOrdinalCloglogRegr` | Borderline | Blocked |
| Ordinal | `InferenceOrdinalKKGLMM` | Difficult | Not a realistic target |
| Proportion | `InferencePropBetaRegr` | Borderline | Blocked, though promising — CRAN's `betareg`-adjacent literature (Cordeiro & Ferrari beta-regression Bartlett corrections, e.g. arXiv 1501.07551) is unpaywalled; worth checking first if this family becomes a priority |
| Proportion | `InferencePropZeroOneInflatedBetaRegr` | Difficult | Not a realistic target |
| Proportion | `InferencePropKKGLMM` | Difficult | Not a realistic target |
| Survival | `InferenceSurvivalCoxPHRegr` | Borderline | Blocked (partial likelihood) |
| Survival | `InferenceSurvivalStratCoxPHRegr` | Borderline | Blocked |
| Survival | `InferenceSurvivalKKLWACoxOneLik` | Borderline | Blocked |
| Survival | `InferenceSurvivalKKStratCoxOneLik` | Borderline | Blocked |
| Survival | `InferenceSurvivalWeibullRegr` | Borderline | Blocked |
| Survival | `InferenceSurvivalDepCensTransformRegr` | Difficult | Not a realistic target |
| Survival | `InferenceSurvivalKKWeibullFrailtyOneLik` | Difficult | Not a realistic target |
| Survival | `InferenceSurvivalKKClaytonCopulaOneLik` | Difficult | Not a realistic target |

**Key takeaway**: "Blocked" applies to nearly every GLM path regardless of theoretical
tier — the theoretical tier says whether Cordeiro-style theory *exists*, not whether it
can be safely reproduced without either a non-paywalled source or an existing verified
implementation to check against. `InferenceContinKKOLSOneLik` was the one exception,
and *not* because Gaussian is "the easiest GLM" in the tensor-algebra sense — it's
because holding `sigma2` fixed turns the problem into a textbook F/t-test identity with
no tensor algebra involved at all, a structurally different (and easier) situation than
every other row in this table.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Easy** | Standard smooth binary GLM. Bartlett-style LR calibration is classical and the current LR plumbing fits it well. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Easy** | Smooth binary probit likelihood with fixed dispersion. Harder than logit, but still a regular low-dimensional likelihood path. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same Poisson likelihood as ordinary count Poisson. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same companion Poisson likelihood engine. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Easy** | Smooth binomial likelihood with a noncanonical link. Still a regular fixed-effects GLM-type target. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Easy** | Smooth fixed-effects binomial likelihood; less canonical, but still structurally regular. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | The path is smooth, but it is already a bespoke combined-likelihood construction rather than a textbook LR model. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | Quadrature-integrated random-effects likelihood makes analytic higher-order LR correction highly bespoke. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **Difficult** | Hybrid conditional-logit plus GLMM structure is not a realistic package-wide Bartlett target. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Easy** | Canonical smooth GLM and a natural first implementation target. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Easy** | Reported variance is robust, but the LR path is still the ordinary Poisson likelihood. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Easy** | Same comment: the LR path is still Poisson. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Smooth likelihood, but the dispersion parameter makes higher-order LR corrections more bespoke. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Mixture likelihood with inflation block; a stable analytic Bartlett factor is not a clean target. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **Difficult** | Truncation and two-part structure push this into model-specific research territory. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **Difficult** | Mixture plus dispersion parameter is too bespoke for a general Bartlett rollout. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **Difficult** | Same issue with added hurdle structure. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **Difficult** | Quadrature GLMM path; not a practical generic analytic Bartlett implementation. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **Difficult** | Hurdle plus random effects plus quadrature is too structurally complex. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **Difficult** | Custom combined likelihood with nonstandard LR calibration needs. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy — ✅ implemented** | Gaussian likelihood is the cleanest case; turned out to be an exact F/t-test identity (see Status Update above), not a Cordeiro tensor derivation. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Smooth and Gaussian, but variance components make the LR correction more bespoke than OLS. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Regular fixed-effects cumulative-link model, but threshold parameters make the higher-order terms nontrivial. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same threshold issue with a probit link. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Smooth fixed-effects likelihood, but link-specific higher-order work is needed. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same issue as above. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **Difficult** | Ordinal thresholds plus quadrature plus variance parameters are not a realistic first Bartlett tier. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Smooth likelihood, but joint mean/precision structure complicates the Bartlett factor. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **Difficult** | Three-component mixture model. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **Difficult** | Same quadrature GLMM issue as the incidence counterpart. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Partial-likelihood Bartlett corrections can exist, but they are much less plug-in than ordinary GLMs. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same Cox issue with strata-specific risk sets. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Combined design usage of Cox partial likelihood still needs model-specific higher-order treatment. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same with stratification. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival likelihood, but materially more bespoke than the GLM cases. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **Difficult** | Highly custom likelihood with coupled event/censoring structure. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **Difficult** | Frailty integration and extra variance structure make analytic Bartlett correction impractical. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **Difficult** | Copula dependence plus parametric margins plus combined design structure. |

## Easy Tier

These are the paths where Bartlett correction fits the current architecture best:

1. `InferenceIncidLogRegr`
2. `InferenceIncidProbitRegr`
3. `InferenceCountPoisson`
4. `InferenceIncidModifiedPoisson`
5. `InferenceIncidKKModifiedPoisson`
6. `InferenceIncidLogBinomial`
7. `InferenceIncidBinomialIdentityRiskDiff`
8. `InferenceCountRobustPoisson`
9. `InferenceCountQuasiPoisson`
10. `InferenceContinKKOLSOneLik`

These are the best candidates because:

- the raw LR test already fits the existing `InferenceAsympLik` machinery
- the models are smooth and low-dimensional enough that higher-order LR
  calibration is plausible
- the corrected LR p-value can be inverted using the same CI plumbing

## Borderline Tier

These are the families where Bartlett correction is mathematically plausible but
not close to plug-and-play:

- negative binomial
- fixed-effects ordinal cumulative-link models
- beta regression
- Cox and stratified Cox partial-likelihood paths
- Weibull regression
- Gaussian LMM
- the combined conditional-logistic likelihood path

The recurring complications are:

- nuisance dispersion or variance components
- threshold parameters
- partial likelihood rather than full likelihood
- bespoke combined-likelihood constructions

I would not schedule these until the easy tier is working and the base-class API
for correction factors has stabilized.

## Difficult Tier

These are not good targets for a generic package-wide Bartlett rollout:

- zero-inflated and hurdle count models
- zero/one-inflated beta
- all quadrature GLMM likelihood paths
- frailty and copula models
- dependent-censoring transform likelihood
- custom combined-likelihood hybrids like `cPoisson` and `clogit + GLMM`

For these families, “Bartlett correction” quickly becomes one of:

- a family-specific research implementation
- a simulation-based LR calibration
- a bootstrap correction rather than a clean analytic Bartlett factor

That is outside the spirit of a common `InferenceAsympLik` extension.

## Recommended Implementation Plan

### Phase 1: Base-Class Plumbing

Add optional Bartlett support to `InferenceAsympLik`:

1. new testing type: `lik_ratio_bartlett`
2. new dispatcher methods for p-values and CIs
3. new subclass hook:
   - `supports_bartlett_likelihood_ratio()`
   - `get_bartlett_factor(...)` or `get_bartlett_lr_adjustment(...)`
4. reuse the existing LR memoization and null-fit warm-start machinery

This phase is a moderate refactor but does not require any family-specific math
yet.

### Phase 2: Easy Families First

Implement and test Bartlett factors for:

1. ordinary incidence logit
2. ordinary incidence probit
3. Poisson
4. modified-Poisson incidence wrappers
5. log-binomial
6. identity-binomial
7. Gaussian OLS one-likelihood path

This is the tier most likely to produce value without destabilizing the package.

### Phase 3: Selected Borderline Families

Only after the easy tier is stable, consider:

1. ordinal fixed-effects models
2. Cox / stratified Cox
3. negative binomial
4. Weibull

I would leave beta regression, LMMs, and combined clogit likelihoods for later,
if at all.

## Bottom Line

Bartlett correction is a realistic extension for the package’s **simplest smooth
likelihood paths**, and it integrates naturally with the current LR p-value and
CI inversion architecture.

The correct place for the shared logic is `InferenceAsympLik`, but only as
**plumbing**. The correction factor itself remains family-specific.

Pragmatically:

- **Easy**: GLM-like fixed-effects likelihoods and Gaussian OLS
- **Borderline**: ordinal, Cox, Weibull, negative binomial, beta, LMM, and
  bespoke combined fixed-effects paths
- **Difficult**: mixtures, quadrature GLMMs, frailty, copula, and other custom
  latent-structure likelihoods

So if the package wants Bartlett-corrected LR tests and inverted CIs, the right
strategy is a **selective rollout**, not a simultaneous all-path implementation.
