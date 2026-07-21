# Closed-Form Firth-Gradient Audit For Likelihood Paths

## Scope And Decision Rule

This audit covers the package’s likelihood-backed inference paths, grouped by the
**actual likelihood engine** they use rather than by every wrapper class that
calls the same score/Hessian code.

Operationally, I audited the paths that participate in the package’s
likelihood-backed testing/inversion framework, i.e. the classes/files that
implement `get_likelihood_test_spec()`. That is the relevant surface for the
question “can LBFGS continue uninterrupted under a Firth penalty gradient?”.

This means the audit is complete for the current **likelihood-test paths**. It
does **not** claim to cover every asymptotic model in the package, because some
classes use Wald/GEE/pass-through logic without a likelihood-test spec.

The expanded table below lists **every concrete class** on that
likelihood-backed surface. It excludes:

- abstract base classes
- alias names that point to the same concrete class object
- non-likelihood IVWC/GEE/Wald/pass-through classes without `get_likelihood_test_spec()`

The question audited is:

> Does this likelihood path admit a **closed-form Firth / Jeffreys penalty
> gradient** that is realistic enough that the current **L-BFGS** workflow could
> continue without switching to a derivative-free optimizer?

I use three labels:

- **Yes**: a closed-form gradient is standard or straightforward enough that an
  analytic adjusted score / Jeffreys-penalty gradient is realistic.
- **Borderline**: a closed-form gradient exists in principle, but would require
  a bespoke derivation that is algebraically heavy enough that I would not count
  it as a clean “LBFGS continues uninterrupted” path without dedicated work.
- **No**: the path is mixture-, latent-, quadrature-, copula-, or custom
  combined-likelihood enough that a practical closed-form Firth gradient is not
  a credible generic implementation target.

This is a **practical engineering audit**, not a statement about abstract
mathematical existence for arbitrary symbolic differentiation.

## Audit Table

### Incidence

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceIncidLogRegr` | `fast_logistic_regression_cpp` | **Yes** | Standard Bernoulli-logit GLM; classical Firth setting with realistic analytic adjusted score. |
| `InferenceIncidProbitRegr` | `fast_ordinal_probit_regression_cpp` in the 2-category case | **Yes** | Smooth fixed-effects binary probit likelihood; closed-form Jeffreys/Firth gradient is realistic in the same sense as other simple binary GLMs. |
| `InferenceIncidModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson likelihood engine as ordinary Poisson regression. |
| `InferenceIncidKKModifiedPoisson` | `fast_poisson_regression_cpp` | **Yes** | Same Poisson companion likelihood engine as above. |
| `InferenceIncidLogBinomial` | `fast_log_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with noncanonical link; analytic Firth gradient is still realistic. |
| `InferenceIncidBinomialIdentityRiskDiff` | `fast_identity_binomial_regression_cpp` | **Yes** | Smooth binomial likelihood with fixed dispersion; less pleasant algebra than logit, but still tractable. |
| `InferenceIncidKKClogitOneLik` | stacked conditional-logistic + reservoir-logistic path via `fast_logistic_regression_with_var_cpp` | **Borderline** | Smooth logistic-based combined likelihood, but not the ordinary unconditional logit Firth case. |
| `InferenceIncidKKGLMM` | `fast_logistic_glmm_cpp` | **No** | GH-quadrature integrated random-effects likelihood with variance parameter. |
| `InferenceIncidKKClogitPlusGLMMOneLik` | `fast_clogit_plus_glmm_cpp` | **No** | Hybrid of conditional logit and quadrature GLMM pieces with shared coefficients. |

### Count

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceCountPoisson` | `fast_poisson_regression_cpp` | **Yes** | Canonical Poisson GLM with explicit score/Hessian. |
| `InferenceCountRobustPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is robust/sandwich, but the likelihood-test path is still plain Poisson. |
| `InferenceCountQuasiPoisson` | `fast_poisson_regression_cpp` | **Yes** | Reported estimator is quasi, but the likelihood-test path is still plain Poisson. |
| `InferenceCountNegBin` | `fast_neg_bin_cpp`, `fast_neg_bin_with_var_cpp` | **Borderline** | Smooth full likelihood, but dispersion-parameter derivatives make the Firth gradient bespoke and polygamma-heavy. |
| `InferenceCountZeroInflatedPoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Mixture likelihood with count and inflation blocks; practical analytic Firth gradient is not a clean target. |
| `InferenceCountHurdlePoisson` | `fast_zero_augmented_poisson_cpp` | **No** | Truncation/mixture structure makes the Jeffreys penalty gradient highly bespoke. |
| `InferenceCountZeroInflatedNegBin` | `fast_zinb_cpp` | **No** | Mixture plus dispersion parameter makes closed-form Firth support impractical. |
| `InferenceCountHurdleNegBin` | `fast_hurdle_negbin_cpp` | **No** | Same issue as above with additional hurdle/truncation structure. |
| `InferenceCountKKGLMM` | `fast_poisson_glmm_cpp` | **No** | Random-effects quadrature likelihood; not a practical closed-form Firth path. |
| `InferenceCountKKHurdlePoissonOneLik` | `fast_hurdle_poisson_glmm_cpp` | **No** | Truncated count + random effects + quadrature is too structurally complex. |
| `InferenceCountKKCPoissonOneLik` | `fast_cpoisson_combined_with_var_cpp` | **No** | Hybrid conditional-plus-marginal likelihood; combined information penalty is bespoke. |

### Continuous

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceContinKKOLSOneLik` | `fast_ols_with_var_cpp` | **Yes** | Gaussian likelihood has explicit information and simple matrix derivatives. |
| `InferenceContinKKGLMM` | `fast_gaussian_lmm_cpp` | **Borderline** | Gaussian structure helps, but variance-component derivatives mean this is no longer a simple drop-in Firth path. |

### Ordinal

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceOrdinalPropOddsRegr` | `fast_ordinal_regression_cpp` | **Borderline** | Fixed-effects ordinal likelihood with thresholds; analytic path exists but needs dedicated derivation. |
| `InferenceOrdinalOrderedProbitRegr` | `fast_ordinal_probit_regression_cpp` | **Borderline** | Same threshold issue as above, with link-specific derivation for probit. |
| `InferenceOrdinalCauchitRegr` | `fast_ordinal_cauchit_regression_cpp` | **Borderline** | Smooth ordinal likelihood, but link-specific adjusted-score derivation is needed. |
| `InferenceOrdinalCloglogRegr` | `fast_ordinal_cloglog_regression_cpp` | **Borderline** | Same as above for the cloglog link. |
| `InferenceOrdinalKKGLMM` | `fast_ordinal_glmm_cpp` | **No** | Cutpoints plus quadrature plus variance parameter make this too bespoke. |

### Proportion

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferencePropBetaRegr` | `fast_beta_regression_cpp` | **Borderline** | Smooth likelihood, but mean-plus-precision structure makes the Jeffreys adjustment nontrivial. |
| `InferencePropZeroOneInflatedBetaRegr` | `fast_zero_one_inflated_beta_cpp` | **No** | Three-component mixture; not a realistic clean analytic Firth target. |
| `InferencePropKKGLMM` | `fast_logistic_glmm_cpp` | **No** | Same engine and same quadrature/variance-parameter issue. |

### Survival

| Concrete likelihood-based inference path | Engine / fitter | Audit result | Why |
|---|---|---|---|
| `InferenceSurvivalCoxPHRegr` | `fast_coxph_regression_cpp` | **Borderline** | Cox bias reduction is plausible, but risk-set derivatives make this a bespoke implementation rather than a plug-in GLM case. |
| `InferenceSurvivalStratCoxPHRegr` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same Cox issue, with extra stratification structure. |
| `InferenceSurvivalKKLWACoxOneLik` | `fast_coxph_regression_cpp` | **Borderline** | Uses Cox partial likelihood over combined data; analytic bias reduction is plausible but still bespoke. |
| `InferenceSurvivalKKStratCoxOneLik` | `fast_stratified_coxph_regression_cpp` | **Borderline** | Same as above with strata-specific risk sets. |
| `InferenceSurvivalWeibullRegr` | `fast_weibull_regression_cpp` | **Borderline** | Smooth parametric survival likelihood, but materially more bespoke than the GLM cases. |
| `InferenceSurvivalDepCensTransformRegr` | `fast_dep_cens_transform_optim_cpp` | **No** | Bespoke transformation likelihood with coupled event/censoring parameter blocks. |
| `InferenceSurvivalKKWeibullFrailtyOneLik` | `fast_weibull_frailty_cpp` | **No** | Frailty integration and variance parameter make analytic Firth support impractical. |
| `InferenceSurvivalKKClaytonCopulaOneLik` | `fast_clayton_weibull_aft_optim_cpp` | **No** | Copula dependence parameter plus Weibull margins plus combined design structure. |

## Concrete Conclusions

### Clean “yes” paths

These are the paths where I would expect a closed-form Firth gradient to be a
realistic extension that preserves the L-BFGS workflow:

- Bernoulli logit GLM
- Poisson GLM
- log-binomial GLM
- binomial identity-link GLM
- Gaussian linear model

These are the best first targets if the goal is to keep the present optimizer
stack and avoid numerical differentiation of the Jeffreys penalty.

### Borderline paths

These paths are analytically smooth enough that a closed-form Firth gradient is
not impossible, but I would not treat them as “LBFGS just keeps going” without
substantial model-specific derivation work:

- negative-binomial regression
- beta regression
- fixed-effects ordinal cumulative-link models, including ordered probit
- Cox / stratified Cox / LWA Cox
- Weibull regression
- Gaussian LMM
- conditional logistic matched-pair combined likelihood

In practice I would split these into two subgroups:

1. **likely worth it**:
   logit/Poisson-adjacent models, some ordinal models, maybe Cox
2. **probably not worth it early**:
   beta, negative binomial with dispersion, Gaussian LMM

### Clear “no” paths

These paths are too mixture- or latent-structure-heavy for a practical generic
closed-form Firth gradient:

- zero-inflated / hurdle count models
- zero/one-inflated beta
- all quadrature GLMM paths
- hurdle Poisson GLMM combined likelihood
- cPoisson combined likelihood
- clogit + GLMM hybrid
- frailty models
- copula models
- dependent-censoring transformation model

For these, a Firth implementation would either become:

- a bespoke research project per family, or
- a numerical penalty-gradient approximation, which defeats the goal of letting
  L-BFGS continue cleanly.

## Recommendation

If the package wants Firth support while preserving the current L-BFGS-based
architecture, I would limit the first implementation set to:

1. `InferenceIncidLogRegr`
2. `InferenceCountPoisson`
3. `InferenceIncidModifiedPoisson` and the Poisson companion-likelihood paths
4. `InferenceIncidLogBinomial`
5. `InferenceIncidBinomialIdentityRiskDiff`
6. `InferenceContinKKOLSOneLik`

After that, the next tier worth evaluating would be:

7. fixed-effects ordinal cumulative-link models, including `InferenceOrdinalOrderedProbitRegr`
8. Cox / stratified Cox
9. `InferenceIncidKKClogitOneLik`
10. maybe negative binomial

I would not plan on package-wide Firth support across the quadrature, mixture,
copula, frailty, and custom combined-likelihood engines if the requirement is
“closed-form gradient so L-BFGS continues uninterrupted.”

## Appendix: Closed-Form S3 Risk-Set Moment Derivation For Cox PH

This appendix works out the specific "bespoke" algebra referenced in the
`InferenceSurvivalCoxPHRegr` / `InferenceSurvivalStratCoxPHRegr` rows above,
so that a future implementation does not have to re-derive it. It follows the
Breslow tie-handling and risk-set-sum structure already implemented in
[EDI/src/fast_coxph_regression.cpp](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_coxph_regression.cpp)
(the `S0`/`S1`/`S2` accumulation loop around lines 154–214), and shows exactly
what a third moment tensor `S3` would need to add.

### Existing risk-set moments

For each distinct event time `t_k` with risk set `R_k` and Breslow event count
`d_k`, the code already accumulates, per covariate vector `x_i` and Cox weight
`w_i = exp(x_i'β)`:

```
S0(β,t_k) = Σ_{i∈R_k} w_i                       (scalar)
S1(β,t_k) = Σ_{i∈R_k} w_i x_i                    (p-vector)
S2(β,t_k) = Σ_{i∈R_k} w_i x_i x_i'               (p×p matrix)
```

which is exactly `workspace.S1` / `workspace.S2` in the code. From these the
code forms the risk-set mean `x̄_k = S1/S0` (`workspace.e_z`) and the risk-set
weighted covariance `V_k = S2/S0 - x̄_k x̄_k'` (the `hess` accumulation at line
~206–214). Summed with Breslow weights `d_k`, these give the ordinary score
and observed information already returned as `fisher_information`:

```
U(β)  = Σ_k [ Σ_{j∈D_k} x_j - d_k x̄_k ]
I(β)  = Σ_k d_k V_k(β)
```

### Why this is a cumulant-generating-function recursion

`S0(β,t_k) = Σ_{i∈R_k} exp(x_i'β)` is literally the cumulant generating
function (in `β`) of the discrete distribution that puts weight
`exp(x_i'β)/S0` on each risk-set member's covariate vector `x_i` (a per-risk-set
softmax). That means successive `β`-derivatives of `log S0(β,t_k)` are exactly
the successive cumulants of that distribution:

```
∂  log S0 / ∂β        =  x̄_k          (1st cumulant = mean)
∂² log S0 / ∂β∂β'     =  V_k          (2nd cumulant = covariance)
∂³ log S0 / ∂β∂β∂β    =  M3_k         (3rd cumulant = third central moment)
```

`U(β)` and `I(β)` are Breslow-weighted sums of the 1st and 2nd of these. The
Firth penalty gradient needs the derivative of `I(β)`, so it needs the natural
next term in the same recursion: the 3rd cumulant, `M3_k`.

### The S3 moment and the resulting M3 tensor

Define the third weighted moment tensor over the risk set, the direct
extension of `S1`/`S2`:

```
S3(β,t_k) = Σ_{i∈R_k} w_i · (x_i ⊗ x_i ⊗ x_i)     (p×p×p tensor)
```

i.e. `S3_{jlm} = Σ_{i∈R_k} w_i x_{i,j} x_{i,l} x_{i,m}`, computed with the same
incremental risk-set scan already used for `S1`/`S2` (same loop, one more
nested index, same `O(n·p³)` pass over the sorted event/censoring times).

The 3rd cumulant (= 3rd central moment, since this is a natural exponential
family) works out to:

```
M3_{jlm}(β,t_k) = S3_{jlm}/S0 − x̄_j V_{lm} − x̄_l V_{jm} − x̄_m V_{jl} − x̄_j x̄_l x̄_m
```

with `V_{jl} = V_k(β)_{jl}` and `x̄ = x̄_k(β)` as already computed. Equivalently,
`M3_k = ∂V_k/∂β` component-wise: differentiating the risk-set covariance with
respect to each `β_m` reproduces this same tensor — a consistency check on the
formula above.

### Assembling the Firth-adjusted score for Cox

Summing with the same Breslow event-count weights `d_k` used for `I(β)`:

```
∂I(β)/∂β_m = Σ_k d_k · M3_k(β)[:,:,m]              (p×p matrix, one per m)
```

and the Firth penalty gradient follows the generic trace identity:

```
∂/∂β_m log|I(β)| = tr( I(β)^{-1} · ∂I(β)/∂β_m )
```

giving the adjusted score used for the penalized fit / profile score test:

```
U*(β) = U(β) + 0.5 · [ tr(I⁻¹ Σ_k d_k M3_k[:,:,1]), ... , tr(I⁻¹ Σ_k d_k M3_k[:,:,p]) ]
```

This matches the closed-form correction in Heinze & Schemper (2001,
*Biometrics*, "A solution to the problem of monotone likelihood in Cox
regression"), the basis of the `coxphf` R package, confirming the algebra
above is the established Firth-Cox result rather than a novel derivation.

### Why this stays "Borderline" and not "Yes"

- It is a genuine closed form — no numerical third-derivative differentiation
  of `I(β)` is required, so L-BFGS can in principle continue uninterrupted.
- But it requires a new `S3` accumulator (a `p×p×p` tensor) alongside the
  existing `S0`/`S1`/`S2` accumulators, at `O(p³)` storage/compute per risk
  set versus `O(p²)` for the existing Hessian pass — a real cost and code
  change, not a drop-in reuse of the GLM exponential-family shortcut that
  makes logit/Poisson/OLS "Yes".
- `InferenceSurvivalStratCoxPHRegr` / `InferenceSurvivalKKStratCoxOneLik`
  need the same derivation applied independently within each stratum's risk
  sets, and `InferenceSurvivalKKLWACoxOneLik` needs it applied over the
  combined-likelihood risk sets — mechanically the same recursion, but each
  is its own bespoke wiring into the corresponding fitter rather than a
  single shared implementation.
- Efron tie-handling (as opposed to the Breslow approximation assumed above,
  which matches the current `fast_coxph_regression.cpp` implementation) would
  add further per-tie correction terms to `S1`/`S2`/`S3` that are not derived
  here.
