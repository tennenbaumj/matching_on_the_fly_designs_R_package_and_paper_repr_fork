# TODO: Vargas-Ferrari-Lemonte Bartlett-Type Gradient-Test Correction

## Scope

This is an implementation plan for adding a higher-order-corrected version of the
**gradient test** (Terrell, 2002) — the fourth member of the "trinity of tests"
(Wald, score, LR, gradient) already computed by this package as `testing_type = "gradient"`
— alongside the Bartlett-corrected LR (`lik_ratio_bartlett_approx`/`_exact`, already
implemented) and the planned Cordeiro-Ferrari modified score test
(`package_metadata/score_correction_cordeiro_ferrari.md`).

Read those two documents first; this plan reuses their architecture almost verbatim.

## Background

The gradient statistic is `S = U(θ̃)'(θ̂ - θ̃)`, where `θ̃` is the null-restricted MLE and
`θ̂` is the unrestricted MLE — a score-times-estimate-gap product, needing no information
matrix at all. Like Wald, score, and LR, it is only χ²_q to first order (`O(n⁻¹)` error).

**Source**: Vargas, T.M., Ferrari, S.L.P., Lemonte, A.J. (2013), *"Gradient statistic:
higher-order asymptotics and Bartlett-type correction"*, arXiv:1206.2206 (published,
*Electronic Journal of Statistics*, open access). Unlike the Cordeiro (1983) logistic-LR
search, **I obtained the arXiv LaTeX source directly** (`arxiv.org/e-print/1206.2206`),
not a lossy AI-summarized PDF extraction — every equation below is transcribed verbatim
from that source, not recalled from memory or OCR'd.

### The correction, verbatim

For testing `H0: θ₁ = θ₁₀` (`q` of `p` parameters, `p-q` nuisance), with `Uⱼ = ∂ℓ/∂θⱼ`,
`Uⱼᵣ = ∂²ℓ/∂θⱼ∂θᵣ`, etc., and cumulants `κⱼᵣ = E(Uⱼᵣ)`, `κⱼᵣₛ = E(Uⱼᵣₛ)`,
`κⱼᵣₛᵤ = E(Uⱼᵣₛᵤ)` (plus several mixed forms — see the paper's Theorem 1 for the full
list), the CDF of `S` expands as:

```
Pr(S ≤ x) = G_q(x) + (1/24n) Σ_{i=0}^{3} R_i · G_{q+2i}(x) + o(n⁻¹)
```

with `R₁=3A₃-2A₂+A₁`, `R₂=A₂-3A₃`, `R₃=A₃`, `R₀=-(R₁+R₂+R₃)`, and `A₁,A₂,A₃` given by full
tensor (Einstein-summation) expressions over the cumulants above (Theorem 1 in the
source; reproduce verbatim when implementing — do not re-derive from memory).

**The corrected statistic** (Corollary, directly citable):

```
S* = S · {1 - (c + b·S + a·S²)}
a = A₃ / (12·n·q·(q+2)·(q+4))
b = (A₂ - 2·A₃) / (12·n·q·(q+2))
c = (A₁ - A₂ + A₃) / (12·n·q)
```

`S*` is χ²_q to `o(n⁻¹)`.

**Notice this is the exact same functional form as Cordeiro-Ferrari's modified score
correction** (`score_correction_cordeiro_ferrari.md`) — same `a,b,c` formulas in terms of `A₁,A₂,A₃`. This
is not a coincidence: the paper states explicitly that whenever a statistic's CDF admits
this `G_q(x) + (1/24n)ΣR_i G_{q+2i}(x)` expansion — true for score, gradient, and (with a
degenerate simplification) LR — Cordeiro & Ferrari's (1991) general polynomial-correction
result applies uniformly. **This means the final "turn A₁,A₂,A₃ into a corrected
statistic" step is shared, reusable plumbing across score and gradient** — only the
per-statistic computation of `A₁,A₂,A₃` from the model's cumulants differs. Bartlett's
1937 LR correction is the odd one out: LR's own bias term happens to be `S`-independent
(a plain constant), which is why it collapses to a simple rescaling `LR/(1+B/n)` instead
of a quadratic-in-`S` polynomial — a genuine structural fact about LR, not a
simplification we chose.

**Practical implication for our codebase**: implement one shared base-class helper
```r
apply_bartlett_type_polynomial_correction(S, A1, A2, A3, q, n)
```
usable by *both* `score_modified_exact` and `gradient_modified_exact` (see Phase 1).

### What's directly usable vs. what needs adaptation

- The paper gives **fully worked, self-verified one-parameter examples** (exponential,
  normal, inverse-normal, gamma, Pareto, Laplace, truncated extreme value) where the
  approximate moments derived from the general theorem are checked against **exact**
  moments computed from the known distribution (e.g. the exponential case: `n·x̄` is
  exactly Gamma-distributed, so exact vs. approximate moments are compared directly).
  These are reproducible unit tests for the *general machinery*, not GLM-specific, but
  they're a strong transcription-correctness check before touching any regression model.
- A "two orthogonal parameters" section gives a nuisance-parameter case with its own
  Monte Carlo validation (Birnbaum-Saunders example) — **does not directly apply** to our
  regression coefficients (real covariates aren't orthogonal in general), but confirms the
  paper's own validation methodology, worth mirroring.
- **No GLM-specific worked example exists in this paper** (checked — despite the source
  filename being `BartlettGrad_Bernoulli.tex`, there is no Bernoulli/binomial/logistic
  content in the final text; that filename is presumably a stale internal working title).
  Unlike Cordeiro-Ferrari's score correction, **there is no `mdscore`-equivalent CRAN
  package to numerically diff against for this one** — the confidence gain here is
  entirely "verbatim equations, not memory/OCR," not "an existing reference
  implementation." Adapting Theorem 1's general `p`-parameter tensor formula to a
  specific GLM (logistic regression) is still original derivation work, just with a
  citable, exactly-transcribed starting point instead of a remembered one.

## Proposed API surface (mirrors Bartlett and Cordeiro-Ferrari exactly)

- Testing types: `"gradient_modified_approx"`, `"gradient_modified_exact"`
- Public methods on `InferenceAsympLik`:
  - `compute_gradient_modified_approx_two_sided_pval(delta = 0, B = 99)`
  - `compute_gradient_modified_approx_confidence_interval(alpha = 0.05, B = 99)`
  - `compute_gradient_modified_exact_two_sided_pval(delta = 0)`
  - `compute_gradient_modified_exact_confidence_interval(alpha = 0.05)`
  - `compute_gradient_modified_two_sided_pval(delta = 0, B = 99)` /
    `compute_gradient_modified_confidence_interval(alpha = 0.05, B = 99)` — "best
    available" smart wrapper, identical exact-wins/errors-if-neither/`B`-ignored-warning
    logic as the other two corrections
- Private hooks:
  - `supports_gradient_modified_ratio_approx()` / `supports_gradient_modified_ratio_exact()`
  - `get_gradient_modified_correction_approx(spec, delta, full_fit, null_fit, B = 99)`
  - `get_gradient_modified_correction_exact(spec, delta, full_fit, null_fit)`
  - Same Phase-0 open question as the other two docs: return `(a,b,c)` or the corrected
    statistic directly? Given the shared-polynomial-helper insight above, leaning even
    more strongly toward `(A1, A2, A3, q)` (or `(a,b,c)` after the shared helper applies
    the `12n·q·(q+2)·(q+4)`-type denominators) so `apply_bartlett_type_polynomial_correction()`
    is the *only* place the final combination logic lives, shared with score-modified

## Implementation phases

### Phase 0: Naming/architecture confirmation
- [ ] Confirm testing-type names/aliases
- [ ] Confirm hook return shape, and confirm the shared
      `apply_bartlett_type_polynomial_correction(S, A1, A2, A3, q, n)` helper design with
      whoever implements Cordeiro-Ferrari first (whichever of the two correction types
      lands first should build this helper; the second should reuse it, not duplicate it)

### Phase 1: Shared base-class plumbing (`InferenceAsympLik`)
- [ ] Add `gradient_modified_approx`/`gradient_modified_exact` to `normalize_testing_type`,
      `set_testing_type` docs, the two dispatch switches, and the
      `get_supported_testing_types_with_bartlett()`-style wrapper (by this point that
      wrapper should probably be renamed to something like
      `get_supported_testing_types_with_corrections()` since it now needs to fold in
      three independent correction families — Bartlett, score-modified,
      gradient-modified — each with their own approx/exact opt-in flags)
- [ ] Add `get_memoized_likelihood_test_pval()` branches: reuse the **existing raw
      gradient statistic** computation (`gradient_test_from_restricted_score_cpp`,
      already computed for the plain `"gradient"` testing type) as `S`, then apply the
      shared polynomial correction
- [ ] Thread `B` through exactly like the other two corrections, **including the same
      shadowed-override fix** already needed twice now
      (`InferenceAsympLikStdModCache`, `InferenceContinKKOLSOneLik` both override
      `compute_likelihood_test_two_sided_pval()` directly)
- [ ] CI: extend `invert_test_pval_confidence_interval()`'s eligible-testing-types set —
      no new root-finding logic, fully generic already
- [ ] Smart wrapper: copy the exact-wins-over-approx pattern verbatim

### Phase 2: Approx (Monte-Carlo) path — same cost tier as Bartlett-approx, *not* score-approx
- [ ] Implement `get_gradient_modified_correction_approx()` in `InferenceParamBootstrap`
      as a bootstrap-calibrated gradient test: simulate `B` datasets under H0 at `delta`
      (reuse `simulate_under_lik_null()`), compute the raw gradient statistic `S_b` for
      each replicate, empirical tail probability against the observed `S`
- [ ] **Cost note, important distinction from score-modified**: the gradient statistic
      needs *both* the null-restricted fit (for `U(θ̃)`) and the unrestricted fit (for
      `θ̂`, to form the estimate gap) — same two-refits-per-replicate cost as
      Bartlett-approx, unlike score-modified-approx which only needs the null refit. Do
      not assume gradient-approx is as cheap as score-approx.
- [ ] `supports_gradient_modified_ratio_approx()` delegates to
      `supports_lik_ratio_param_bootstrap()` by default, same auto-opt-in pattern, same
      need to check whether the Zero-Inflated-Poisson/Hurdle-Poisson raw-statistic
      miscalibration carve-out applies here too (check the raw gradient statistic
      specifically, don't assume it transfers from the LR/score carve-out reasoning
      without checking)

### Phase 3: Exact (analytic) path
- [ ] **Recommended first target: re-derive and validate the paper's own one-parameter
      worked examples** (exponential, normal, gamma, etc.) as a pure transcription-check,
      *before* touching any GLM — this doesn't require adapting anything to our
      package's model classes, it's a direct test that Theorem 1 was transcribed
      correctly, using the paper's own self-verified exact-vs-approximate moment
      comparisons
- [ ] Only after that passes: adapt Theorem 1's general tensor formula (not the
      orthogonal-parameters special case, which doesn't apply to correlated regression
      covariates) to canonical-link logistic regression — this is original derivation
      work (no GLM-specific worked example or reference implementation exists for this
      one, unlike score-modified's `mdscore`)
- [ ] Validate the logistic-regression adaptation three ways (same as planned for
      score-modified, minus the reference-implementation cross-check which isn't
      available here):
      1. Compare against `gradient_modified_approx` (MC) — should agree within MC noise
      2. Large-scale simulation of `E[S | H0]` and empirical size at nominal α across
         several `(n, β, X)` configurations
      3. Cross-check against `score_modified_exact` if that lands first — both should be
         computing overlapping cumulant quantities (`κ_jrs`, `κ_jrsu`) for the *same*
         logistic-regression likelihood, so the shared intermediate quantities should
         match even though the final polynomials differ
- [ ] `InferenceContinKKOLSOneLik` (Gaussian): same "do this one first" recommendation as
      the other two docs — lowest derivation risk, standard finite-sample theory
      independent of tensor algebra

### Phase 4: Broader rollout & tests
- [ ] Same broad-opt-in-for-approx decision as Bartlett/score-modified
- [ ] `test-gradient-modified-plumbing.R`, `test-gradient-modified-logit.R` (including a
      dedicated test replicating the paper's exponential-distribution worked example as a
      pure transcription check, independent of any GLM machinery),
      `test-gradient-modified-approx-smoke-families.R` — mirror the Bartlett test file
      structure exactly
- [ ] Add to `comprehensive_tests.R` (`"gradient_modified"` test family,
      `ALL_TEST_FAMILY_FILTERS` entry, `supports_gradient_modified`/`_ci` variables)
- [ ] Add `grad-mod-approx`/`grad-mod-exact` column pair to `path_audits_source.R`

## Difficulty table (implementing the modified gradient test per path)

Same 34 paths as the Bartlett and Cordeiro-Ferrari reports. Tiers are identical in
structure to `score_correction_cordeiro_ferrari.md`'s table for the same reason given there: both
corrections need the same regularity conditions (regular, single-linear-predictor,
third/fourth-derivative-differentiable likelihood). Unlike that table, there is no
"verified against an existing package" column here — see the honesty note above.

### Incidence

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | Easy | Recommended second GLM target (after Gaussian OLS); canonical link, `g=0`-style simplifications should carry over from the score-modified derivation once that's validated |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` (2-cat) | Easy | Non-canonical link adds a derivative layer, same as noted for Bartlett-exact |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | Easy | Plain Poisson likelihood underneath |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | Easy | Same, KK sequential design needs leverage-bookkeeping care |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | Easy | Log link + binomial variance; mean/variance derivatives trivial to supply directly |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | Easy | Identity link simplifies several derivative terms (`d²μ=0`) |
| `InferenceIncidKKClogitOneLik` | `fast_logistic_regression_with_var_cpp` | Borderline | Bespoke combined conditional-logit + reservoir likelihood |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | Difficult | Quadrature-integrated random effects |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | Difficult | Hybrid |

### Count

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | Easy | Canonical log-link Poisson |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | Easy | Same underlying likelihood |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | Easy | Same |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | Borderline | Extra jointly-estimated dispersion parameter |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | Difficult | Mixture; check whether the raw-gradient-statistic carve-out applies (Phase 2) |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | Difficult | Same carve-out check |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | Difficult | Mixture + dispersion |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | Difficult | Hurdle + dispersion |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | Difficult | Quadrature GLMM |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | Difficult | Hurdle + GLMM |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | Difficult | Custom combined likelihood |

### Continuous

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Easy — recommended first target overall** | Gaussian; reduces to textbook finite-sample theory, lowest derivation risk in the whole table |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | Borderline | Variance components / quadrature |

### Ordinal

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | Borderline | Multi-threshold cumulative link — Theorem 1's `p`-parameter framework technically covers this (thresholds are just extra parameters), but the resulting tensor algebra grows substantially with `K-1` threshold parameters |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | Borderline | Same |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | Borderline | Same |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | Borderline | Same |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | Difficult | Thresholds + quadrature |

### Proportion

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | Borderline | Joint mean/precision two-parameter family |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | Difficult | Three-component mixture |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | Difficult | Quadrature GLMM |

### Survival

| Path | Engine / Fitter | Difficulty | Notes |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | Borderline | Partial likelihood, not a full-likelihood GLM structure Theorem 1 assumes |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | Borderline | Same, plus strata |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | Borderline | Same, plus combined-design bespoke-ness |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | Borderline | Same |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | Borderline | AFT-style location/scale, closer to NegBin/Beta in spirit |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | Difficult | Highly custom combined event/censoring likelihood |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | Difficult | Frailty integration |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | Difficult | Copula dependence |

**Totals**: Easy = 10, Borderline = 12, Difficult = 12 (34 paths) — same tier counts as
both other reports, for the same structural reason.

## Open risks / questions to resolve before Phase 1 starts

- Same ZA-carve-out question as `score_correction_cordeiro_ferrari.md`: does the
  Zero-Inflated-Poisson/Hurdle-Poisson raw-statistic miscalibration extend to the
  *gradient* statistic specifically? Must be checked independently — don't assume it
  transfers from the LR/score reasoning.
- Confirm the shared `apply_bartlett_type_polynomial_correction()` helper design with
  whichever of {score-modified, gradient-modified} gets implemented first; avoid
  duplicating the `a,b,c` combination logic in two places.
- The paper's general Theorem 1 is for **any** regular parametric model, not
  GLM-specific — unlike Cordeiro-Ferrari's `mdscore`, there's no existing GLM
  implementation to adapt from directly. The one-parameter worked-example
  transcription check (Phase 3, first bullet) is the substitute safety net; treat it as a
  hard gate, not optional, before adapting anything to logistic regression.
- Both `score_modified_exact` and `gradient_modified_exact` will independently need the
  same third/fourth-order cumulant machinery (`κ_jrs`, `κ_jrsu`) for whichever GLM family
  is targeted. Consider extracting a shared "GLM cumulant helper" (leverages, `b'''(θ)`,
  `b''''(θ)`, etc.) once *both* are being implemented, rather than duplicating derivative
  bookkeeping across the two correction types.
