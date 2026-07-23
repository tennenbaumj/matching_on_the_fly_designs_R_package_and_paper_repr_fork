# Cox-Snell Bias-Correction Audit of EDI Likelihood and Partial-Likelihood Paths

Date: 2026-07-23

## Relationship to the existing Cordeiro-McCullagh report

`package_metadata/cordeiro_mccullagh_bias_correction_report.md` already contains
a thorough per-concrete-class audit of explicit first-order bias correction
across every path on the `get_likelihood_test_spec()` surface, using
**exactly** this theory: Cordeiro & McCullagh (1991) is the closed-form,
canonical-link-GLM specialization of the general Cox & Snell (1968) cumulant
expansion. The Easy/Borderline/Difficult tiering in that report is the
authoritative per-class verdict and this document does not re-derive it or
contradict it.

This document adds four things that report does not cover:

1. **File-level grounding** — which `src/*.cpp` kernel and `R/*.R` class each
   path lives in, and what Fisher-information/Hessian infrastructure already
   exists there to build on.
2. **An explicit theoretical account of *why* partial-likelihood (Cox PH,
   conditional logit) and GLMM/frailty (random-effects, Laplace-approximated)
   paths resist the standard Cox-Snell derivation** — the existing report
   labels these Borderline/Difficult but doesn't explain the mechanism.
3. **Engineering effort estimates** (days/weeks) and **where the correction
   would physically live** (R-side post-fit vs. new C++ third-derivative
   kernel).
4. **An explicit ruling on GEE**, which sits outside the
   `get_likelihood_test_spec()` surface the other report scopes to.

Where the two documents overlap, tiers below are translated as: Easy → Easy,
Borderline → Moderate, Difficult → Hard or Very Hard depending on whether a
concrete (if laborious) derivation path exists (Hard) or the family requires
a different theoretical framework altogether (Very Hard / research problem).

## Core formula (recap)

General Cox-Snell bias of MLE component θ̂_s, using expected 2nd/3rd
log-likelihood cumulants and inverse Fisher information i^rs:

  bias(θ̂_s) = Σ_{r,u,v} i^sr i^uv × [(1/2)κ_ruv + κ_ru,v]

For a canonical-link exponential-family GLM, this collapses to the
Cordeiro-McCullagh closed form:

  bias(β̂) = (XᵀWX)⁻¹ Xᵀ W ξ

— a couple of matrix products reusing the working-weight matrix W and the
`XᵀWX` (= observed/Fisher information) that these kernels already compute for
`vcov()`. This is why the canonical GLM family is the cheap case: no
third-derivative cumulants need to be hand-derived or coded.

**Why partial likelihood breaks the derivation.** The Cox-Snell expansion
relies on log-likelihood cumulants being additive over n i.i.d. observations
(κ = Σᵢ κᵢ), which lets you organize the bias as an O(1/n) series. Cox partial
likelihood is a product over *risk sets*, not over independent observations —
each risk-set contribution is correlated with every later one through the
shared at-risk membership, so the ordinary Cox-Snell cumulant-additivity step
does not go through unmodified. A partial-likelihood analog exists in the
survival literature (score/information cumulants computed against the
counting-process martingale representation instead of i.i.d. observations),
but it is a distinct, family-specific derivation — not a drop-in application
of the GLM formula above.

**Why GLMMs break the derivation differently.** The fitted quantity in this
codebase's GLMMs is not a genuine log-likelihood but a Laplace/adaptive-
quadrature *approximation* to the marginal likelihood ∫ L(β, b) φ(b; Σ) db
integrated over random effects b. Cox-Snell bias correction assumes the exact
log-likelihood's cumulants; applying it naively to the Laplace-approximated
surface would correct for the wrong object (finite-sample bias of the
approximate objective, conflated with quadrature/Laplace approximation
error itself, which is typically the dominant error term for small cluster
sizes). A correct treatment requires either bias-correcting the fixed effects
conditional on the estimated variance components under the Laplace
approximation's own error structure, or abandoning Cox-Snell in favor of a
different higher-order correction designed for this setting. This is a
research question, not routine implementation.

**GEE is excluded entirely.** `fast_gee.cpp` optimizes a *quasi*-log-likelihood
(`quasi_loglik`, confirmed in source) via an estimating equation with a
working-correlation sandwich variance — there is no true likelihood whose
cumulants Cox-Snell could act on. Bias correction for GEE estimating
equations is a different (bias-corrected sandwich / small-sample GEE)
literature, out of scope here.

## Summary table

Sorted easiest/highest-value first, then descending difficulty.

| Model family | R class | Kernel file(s) | Type | Difficulty | Where correction lives | Effort |
|---|---|---|---|---|---|---|
| Logistic regression | `InferenceIncidLogRegr` | `fast_logistic_regression.cpp` | Canonical GLM | **Easy** | R-side post-fit; reuses `hess_fisher_info_matrix` / `get_information_matrix()` | ~1 week (build shared helper) |
| Binary probit | `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression.cpp` (2-cat case) | Non-canonical GLM (closed form still exists) | **Easy** | R-side post-fit, shared helper | ~1-2 days once helper exists |
| Poisson (+ modified/robust/quasi variants) | `InferenceCountPoisson`, `InferenceCountRobustPoisson`, `InferenceCountQuasiPoisson`, `InferenceIncidModifiedPoisson`, `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression.cpp` | Canonical GLM | **Easy** | R-side post-fit, shared helper | ~1 day each (5 paths share one fit) |
| OLS / Gaussian | `InferenceContinKKOLSOneLik` | `fast_ols.cpp` | Canonical GLM (Gaussian) | **Easy** | R-side post-fit, shared helper; note bias is typically negligible/zero for the mean-model coefficients in OLS — main value is the σ² bias term | ~1-2 days |
| Log-binomial | `InferenceIncidLogBinomial` | `fast_log_binomial_regression.cpp` | Non-canonical GLM, boundary-sensitive (μ ≤ 1 constraint) | **Moderate** | R-side, but needs boundary guard before applying the correction | ~3-5 days |
| Identity-binomial (risk difference) | `InferenceIncidBinomialIdentityRiskDiff` | `fast_log_binomial_regression.cpp` (identity-link path) | Non-canonical GLM, boundary-sensitive | **Moderate** | R-side, boundary guard required | ~3-5 days |
| Negative binomial | `InferenceCountNegBin` | `fast_negbin_regression.cpp` | Custom 2-block likelihood (β, log-dispersion) | **Moderate** | R-side reuses existing `expected_hessian()` (analytic Fisher already implemented); needs new 3rd-derivative terms for the dispersion block and its cross terms with β | ~1-2 weeks |
| Beta regression | `InferencePropBetaRegr` | `fast_beta_regression.cpp` | Custom 2-block likelihood (mean β, precision φ) | **Moderate** | Needs new analytic 3rd-derivative kernel for mean/precision cross terms | ~1-2 weeks |
| Ordinal cumulative-link (logit/probit/cauchit/cloglog) | `InferenceOrdinalPropOddsRegr`, `InferenceOrdinalOrderedProbitRegr`, `InferenceOrdinalCauchitRegr`, `InferenceOrdinalCloglogRegr` | `fast_ordinal_regression.cpp`, `fast_ordinal_probit_regression.cpp`, `fast_ordinal_cauchit_regression.cpp`, `fast_ordinal_cloglog_regression.cpp` | Custom likelihood, K-1 threshold nuisance params | **Moderate** | New 3rd-derivative kernel per link; threshold-block cross terms multiply with #categories | ~1-2 weeks per link (4 links, some shared structure) |
| Weibull AFT | `InferenceSurvivalWeibullRegr` | `fast_weibull_regression.cpp` | Custom full likelihood (shape + scale/β) | **Moderate** | New 3rd-derivative kernel, 2-block (shape, β) | ~1-2 weeks |
| Cox PH (+ stratified) | `InferenceSurvivalCoxPHRegr`, `InferenceSurvivalStratCoxPHRegr`, `InferenceSurvivalKKLWACoxOneLik`, `InferenceSurvivalKKStratCoxOneLik` | `fast_coxph_regression.cpp` | Partial likelihood | **Moderate–Hard** | Needs partial-likelihood-specific cumulant derivation (see above); existing analytic `hessian()` in the kernel is reusable for the information-matrix half only | ~2-3 weeks; literature search + derivation dominates the effort, not coding |
| KK combined conditional-logit path | `InferenceIncidKKClogitOneLik` | `fast_logistic_regression.cpp` (stacked design) | Custom combined likelihood, no off-the-shelf formula | **Hard** | Would need a bespoke derivation for the stacked conditional+reservoir likelihood | ~3-4 weeks, high risk of dead end |
| Zero-augmented/hurdle Poisson | `InferenceCountZeroInflatedPoisson`, `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson.cpp` | Mixture/truncation, 2 likelihood components | **Hard** | New kernel per component; cross terms between hurdle-prob and count-model blocks are the hard part | ~3-4 weeks |
| ZINB / hurdle negbin | `InferenceCountZeroInflatedNegBin`, `InferenceCountHurdleNegBin` | `fast_zinb.cpp`, `fast_hurdle_negbin.cpp` | Mixture + dispersion, 3 blocks | **Hard** | Same as above plus dispersion cross terms | ~4-5 weeks |
| Zero-one-inflated beta | `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta.cpp` | 3-part mixture (point masses at 0/1 + beta body) | **Hard** | 3+ parameter blocks, dense cross-cumulant table | ~4-5 weeks |
| Custom combined Poisson (KK cPoisson) | `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined.cpp` | Custom combined likelihood, nonstandard nuisance geometry | **Hard** | Bespoke derivation | ~3-4 weeks |
| Dependent-censoring transform | `InferenceSurvivalDepCensTransformRegr` | `fast_survival_models_optim.cpp` | Coupled event/censoring custom likelihood | **Hard** | Bespoke derivation, coupling makes cumulants nonstandard | ~4+ weeks |
| Clayton copula Weibull AFT | `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_survival_models_optim.cpp` | Copula-coupled survival margins | **Very Hard** | Copula dependence parameter adds another nuisance block with no standard closed form | research-scale, no confident estimate |
| All GLMM paths (logistic, Poisson, ordinal, hurdle-Poisson-GLMM, Gaussian LMM) | `InferenceIncidKKGLMM`, `InferenceCountKKGLMM`, `InferencePropKKGLMM`, `InferenceOrdinalKKGLMM`, `InferenceCountKKHurdlePoissonOneLik`, `InferenceContinKKGLMM` | `fast_logistic_glmm.cpp`, `fast_poisson_glmm.cpp`, `fast_ordinal_glmm.cpp`, `fast_hurdle_poisson_glmm.cpp`, `fast_gaussian_lmm.cpp` | Marginal (Laplace-approx) likelihood, random effects | **Very Hard** | Needs a different theoretical framework — see rationale above; not a coding task | research problem |
| Clogit+GLMM hybrid | `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm.cpp` | Hybrid partial-likelihood + GLMM | **Very Hard** | Compounds both the partial-likelihood and GLMM issues above | research problem |
| Weibull frailty | `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty.cpp` | Frailty (random-effect survival), integrated likelihood | **Very Hard** | Same marginal-likelihood issue as GLMM | research problem |
| GEE (Poisson) | `InferenceCountPoissonKKGEE` | `fast_gee.cpp` | Quasi-likelihood, estimating equation | **N/A** | Cox-Snell does not apply; a bias-corrected-sandwich approach is a different project | N/A |

## Excluded, not likelihood paths

Permutation/randomization generators, matching/blocking/rerandomization
design machinery (Pocock-Simon, Atkinson, D-optimal, greedy pair-switching),
bootstrap resampling loops, and all post-fit-only statistics (CMH, Newcombe,
Miettinen-Nurminen, exact tests, rank-based tests without an underlying
parametric likelihood) are excluded — none of these fit a model by MLE or
partial likelihood, so Cox-Snell bias correction has no target parameter to
apply to in this codebase.

## Existing infrastructure worth reusing

- `InferenceMixinInformationMatrix$get_information_matrix()`
  (`R/inference_mixin_information_matrix.R`) already exposes both Fisher and
  observed information for any class implementing
  `get_likelihood_test_spec()`, with an `auto`-preference fallback. Any R-side
  bias-correction helper should sit on top of this rather than re-deriving
  the information matrix.
- `fast_negbin_regression.cpp` and `fast_coxph_regression.cpp` both already
  expose analytic (not numeric) Hessians / expected Hessians
  (`hessian()`, `expected_hessian()`, `hess_fisher_info_matrix`), so the
  information-matrix half of the Cox-Snell formula is free for those two
  families — only the third-derivative cumulant term is new work.
- The `InferenceMixinBartlettApprox`, `InferenceMixinCordeiroFerrariApprox`,
  and `InferenceMixinLemonteGradientApprox` stubs (`R/inference_mixin_*.R`)
  are a *different* higher-order-asymptotics project (score/gradient/LR test
  corrections, not point-estimate bias) but establish the package's existing
  Pattern-1 mixin convention (`list(public=..., private=...)`, spliced via
  `utils::modifyList`) that a new `InferenceMixinCoxSnellBiasCorrection` mixin
  should follow for consistency.

## Recommended implementation order

1. **Canonical-link GLM family first**: logistic, probit (binary), Poisson
   (+ its 4 thin wrappers), OLS. Build one shared `XᵀWX`-based helper against
   `get_information_matrix()`; every subsequent canonical-GLM path is then a
   ~1-day integration, not a new derivation. This matches the existing
   Cordeiro-McCullagh report's Easy tier exactly.
2. **Log-binomial / identity-binomial** next, once the boundary-guard logic
   (μ near 0/1) is worked out once and shared between the two link variants
   that already live in the same file.
3. **Negative binomial and Cox PH** as the first non-GLM targets: both
   already have analytic Fisher-information machinery in their C++ kernels,
   so only the third-derivative term is genuinely new work, and both are
   high-value families (heavily used in this package's count and survival
   inference paths).

Explicitly out of scope for a near-term rollout: all GLMM/frailty/copula
paths (research problem, not an implementation task) and GEE (wrong
statistical framework). This matches the existing Cordeiro-McCullagh report's
recommendation to treat those as a later-stage or separate research project
rather than part of a package-wide bias-correction rollout.
