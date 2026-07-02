# EDI Native Kernel Performance Report

Dates: 2026-05-14 (initial pass), 2026-05-17–18 (v2 follow-up)

---

## Objective

Profile the native C++ kernels in `EDI/src`, rank them by optimization payoff, and document all experiments — retained and reverted. This document is the merged and completed record of `perf_experiments.md` and `perf_experiments_v2.md`.

---

## Method

### Build

```bash
R CMD INSTALL --no-test-load --library=/tmp/edi_profile_lib EDI
```

### Timing

In-process 60-repetition median loop using `proc.time()["elapsed"]`. Reps scaled up for profiling passes (800–3000 where noted).

### Profiling

```bash
perf record -F 199 --output=/tmp/perf_<kernel>.data -- \
  Rscript /tmp/profile_one_kernel.R --lib=/tmp/edi_profile_lib --kernel=<kernel> --reps=<N>
perf report --stdio --no-children --sort symbol -i /tmp/perf_<kernel>.data
```

### Important caveat

The installed `EDI.so` did not provide clean line-level attribution through `perf report`, so the loop ranking is an inference from: (1) per-kernel elapsed time, (2) top sampled symbols, and (3) the surrounding source structure. This is sufficient to rank optimization targets but is loop-level reasoning rather than exact per-line accounting.

---

### July 2026 Comprehensive Re-profiling

A full sweep across all 36 EDI kernels × 2 paths (estimate-only `_est` and full/variance `_var`) was conducted in July 2026 using the unified profiling script `/tmp/edi_kernel_profiler.R`. For each kernel, three file types were generated in `/tmp/`:

- **`perf_<kernel>.data`** — raw sample data from `perf record -F 199`. 67 files total (some kernels have only one path; some older files from the v2 re-profiling run are included).
- **`perf_stat_<kernel>.txt`** — hardware performance counters from `perf stat`, capturing instructions, cycles, branches, and branch-misses. 59 files total. Commands:
  ```bash
  perf stat -e instructions:u,cycles:u,branches:u,branch-misses:u \
    Rscript /tmp/edi_kernel_profiler.R <kernel> 2> /tmp/perf_stat_<kernel>.txt
  ```
- **`perf_annotate_<kernel>.txt`** — per-instruction hotspot attribution from `perf annotate --stdio`. 57 files total (some are 0 bytes where perf could not resolve symbols — typically kernels whose dominant cost is R bytecode, GC, or Eigen expression templates inlined into `EDI.so` without debug info). Commands:
  ```bash
  perf record -F 199 --output=/tmp/perf_<kernel>.data -- \
    Rscript /tmp/edi_kernel_profiler.R <kernel>
  perf annotate --stdio -i /tmp/perf_<kernel>.data > /tmp/perf_annotate_<kernel>.txt
  ```

**Kernel naming convention:** `<type>_<path>` where path is `est` (estimate-only / optimizer loop) or `var` (full path including Hessian/variance). Examples: `logistic_est`, `logistic_var`, `hurdle_p_est`, `hurdle_p_var`.

**IPC interpretation** (from `perf_stat` instructions ÷ cycles):
- IPC < 1.5 — memory-bandwidth-limited; loads/stores stall the pipeline
- IPC 1.5–2.5 — mixed arithmetic + memory
- IPC > 2.5 — arithmetic-bound; FPU/FMA throughput is the bottleneck

**Symbol resolution limits:** `libR.so`, `libm.so`, and Eigen templates inlined into `EDI.so` are visible in annotate output. R overhead appears as `bcEval_loop`, `R_gc_internal`; allocator churn as `_int_malloc`, `cfree`, `__libc_calloc`, `rep stosb` (memset from `operator new`).

---

## Initial Kernel Weights (Estimate-Only, Pre-Optimization Baseline)

| Kernel | Median elapsed (s) | Share of total |
|---|---:|---:|
| `weibull_frailty` | 0.0480 | 29.1% |
| `ordinal_clmm` | 0.0470 | 28.5% |
| `logistic_glmm` | 0.0370 | 22.4% |
| `clogit_plus_glmm` | 0.0185 | 11.2% |
| `poisson_glmm` | 0.0115 | 7.0% |
| `adjacent_category_logit` | 0.0030 | 1.8% |

`weibull_frailty`, `ordinal_clmm`, and `logistic_glmm` together account for 80% of native estimate-only time.

### Non-estimate_only (full-path) weights

| Kernel | Median elapsed (s) |
|---|---:|
| `ordinal_clmm_full` | 0.147 |
| `weibull_frailty_full` | 0.078 |
| `logistic_glmm_full` | 0.048 |
| `poisson_glmm_full` | 0.024 |
| `clogit_plus_glmm_full` | 0.013 |
| `adjacent_category_logit_var` | 0.002 |

`ordinal_clmm_full` is substantially more expensive than estimate-only due to the Hessian/covariance path.

### Broader non-GLMM sweep weights

| Kernel | Median elapsed (s) | Share within sweep |
|---|---:|---:|
| `d_optimal_search` | 2.0345 | 37.6% |
| `zinb_full` | 1.3920 | 25.7% |
| `zap_full` | 1.3350 | 24.7% |
| `negbin_var` | 0.2900 | 5.4% |
| `pair_distance_matrix` | 0.2105 | 3.9% |
| `gaussian_lmm_full` | 0.1420 | 2.6% |
| `coxph_full` | 0.0055 | 0.1% |

---

## Initial Perf Highlights

### Ordinal CLMM

Top symbols: `GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()`, `__ieee754_exp_fma`, `__ieee754_log_fma`, `malloc`, `cfree`, Eigen GEMV.
Allocation overhead (`malloc`/`cfree`) was prominent — a sign of repeated tiny allocations inside the derivative path.

### Logistic GLMM

Top symbols: `__log1p_fma`, `plogis_array_safe`, Eigen GEMV, `log1pexp_array_safe`, `LogisticGLMMObjective::operator()`.
Dominated by logistic math and matrix-vector products. SIMD already active; the headroom is in reducing passes, temporaries, and redundant reductions.

### Weibull frailty

Top symbols: `WeibullFrailtyLikelihood::operator()`, vectorized Eigen packet math, Eigen GEMV, vectorized exponentials.
High-throughput numeric kernel; remaining payoff is from dataflow simplification and fewer materialized temporaries.

### Clogit-plus-GLMM

Top symbols: `__log1p_fma`, `plogis_array_safe`, `log1pexp_array_safe`, Eigen GEMV, `ClogitPlusGLMMObjective::operator()`.
Same structural cost pattern as logistic GLMM.

### Poisson GLMM

Top symbols: `PoissonGLMMObjective::operator()`, Eigen GEMV, Eigen packet math, `log_sum_exp_p`.
Dominated by the repeated node/group sweep and weighted matrix-vector work.

### Adjacent-category logit

Top symbols: `AdjacentCategoryLogitNegLogLik::operator()`, `hessian()`, Eigen dense assignment and outer-product routines.
Locally hot but globally tiny.

### D-optimal search

Top symbols: `d_optimal_search_cpp`, repeated `P(i,j)` and vector coefficient access, negligible `std::shuffle` time.
Not an assembly problem; a data-access and algorithmic-structure problem dominated by the nested swap-search loop.

### Pair-distance matrix

Top symbols: `compute_pair_distance_matrix_cpp`, `Rcpp::Matrix::operator()`, `Rcpp::Matrix::offset`, `Rcpp::Vector::operator[]`.
Dominated by accessor overhead, not floating-point arithmetic. A clear low-risk rewrite target.

### ZINB (initial sweep)

Top symbols: `ZeroInflatedNegBin::hessian`, `Rf_dpsifn`, `__ieee754_pow_fma`, `__ieee754_log_fma`.
Unlike GLMM-family kernels, ZINB is visibly Hessian-heavy. An important exception to the "optimize objective loops first" rule.

### ZAP

Top symbols: `ZeroAugmentedPoisson::hessian`, exponential/logistic scalar operations, Eigen expression overhead.
Similar structure to ZINB but without negative-binomial special functions.

---

## Phase 1: Ordinal CLMM Core Path

### Planned

1. Replace generic ordinal GLMM derivative path in `_glmm_engine.h` for the CLMM logit case.
2. Cache threshold transforms and eliminate repeated per-row `dp` allocation in `fast_ordinal_clmm.cpp`.
3. Preallocate node/derivative workspaces to suppress `malloc`/`cfree` churn.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `ordinal_clmm` | estimate-only fit | 0.0470 | 0.0250 | **keep** threshold/allocation fixes |
| `ordinal_clmm` | dedicated specialization attempt | 0.0785 | 0.0885 | **revert** specialization |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | **keep** direct accumulation |
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0360 | **keep** direct accumulation |

**Kept:**
- `OrdinalLikelihood::log_prob_derivs()` tightened: derivative path no longer reconstructs the full threshold vector on every row/node evaluation.
- Generic ordinal path stopped allocating tiny `dp` vectors inside the innermost loop; now reuses work buffers.
- Logistic and Weibull objective/gradient loops stopped materializing `post_k_expanded`; now accumulate directly at the group level.

**Rejected:**
- Dedicated ordinal CLMM objective bypassing `GLMMObjective` was estimate-stable but benchmarked slower than the optimized generic path.

---

## Phase 2: Shared GLMM Objective-Loop Refactor

### Planned

Refactor logistic, Poisson, Weibull, and clogit-concordant kernels around a shared pattern: node sweep, group posterior weights, direct accumulation without full expanded posterior vectors.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0360 | **keep** |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | **keep** |
| `clogit_plus_glmm` | fit fingerprint preserved | same estimates | same estimates | **keep** |
| `poisson_glmm` | larger KK-shaped fit | 0.0520 | 0.0600 | **revert** |

**Kept:**
- Logistic, Weibull, and clogit-concordant: node sweep + group posterior weights + direct accumulation without dense expanded posterior vectors.
- Adjacent passes fused so group log terms, residuals, and weighted reductions are reused.
- Persistent scratch buffers for group layout metadata, node work matrices, and reusable residual/work vectors.
- A shared contiguous-group layout helper introduced in `EDI/src/_helper_functions.h`.

**Rejected:**
- Poisson GLMM version of the same rewrite preserved the math but benchmarked slower — reverted.
- A more aggressive single abstraction across all four kernels was not adopted.

---

## Phase 3: Math and Reduction Tightening

### Planned

Audit repeated `exp`, `log1p`, and `log_sum_exp` calls; cache node factors; reduce redundant segment extraction and temporary `Eigen::VectorXd` creation.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | 0.0430 | 0.0390 | **keep** |
| `logistic_glmm` | Hessian | 0.0070 | 0.0030 | **keep** |
| `weibull_frailty` | neg log-likelihood eval | 0.000680 | 0.001985 | **revert** |
| `weibull_frailty` | Hessian | 0.00330 | 0.00744 | **revert** |
| `clogit_plus_glmm` | estimate-only fit | 0.0500 | 0.0360 | **keep** |
| `clogit_plus_glmm` | Hessian | 0.01015 | 0.00515 | **keep** |

**Estimate consistency verified:**
- logistic GLMM: parameters and negative log-likelihood matched
- Weibull frailty: `nll = 5407.2669016122754`, `sum(diag(H)) = 316963.61731134733`
- clogit-plus-GLMM: `nll = 398.6283207650863`, `param_sum = -1.7298381109594083`, `sum(diag(H)) = -566.51882986578255`

**Kept:** logistic GLMM and clogit-plus-GLMM Phase 3 math tightening.
**Rejected:** Weibull frailty Phase 3 (preserved math, slower on both objective and Hessian paths).

---

## Phase 4: Hessian-Path Optimization and Estimate-Only Early Returns

### Planned

Optimize Hessian accumulation in logistic, Poisson, Weibull, and clogit-plus-GLMM; add internal early-return paths in native fitters for `estimate_only = TRUE` callers.

### Caller audit finding

Many native fit calls already know they only need point estimates. The highest-value partial-information mode was an internal early-return path rather than a new R-facing API.

### Outcome

| Kernel | Path | Before (s) | After (s) | Estimate parity | Decision |
|---|---:|---:|---:|---|---|
| `logistic_glmm` | estimate-only fit | 0.0090 | 0.0060 | unchanged | **keep** |
| `logistic_glmm` | Hessian | 0.0010 | 0.0010 | unchanged | **keep** |
| `clogit_plus_glmm` | estimate-only fit | 0.0230 | 0.0200 | unchanged | **keep** |
| `clogit_plus_glmm` | Hessian | 0.00175 | 0.00170 | unchanged | **keep** |
| `poisson_glmm` | estimate-only fit | 0.0020 | 0.0010 | unchanged | **keep** |
| `poisson_glmm` | Hessian cleanup attempt | 0.0010 | 0.0010 | unchanged | **reject rewrite; keep early return only** |
| `weibull_frailty` | Hessian tightening attempt | 0.00330 | 0.00744 | unchanged | **keep reverted baseline** |

**Kept:**
- `fast_logistic_glmm_cpp`: early return after optimization when `estimate_only = TRUE`; Hessian blockwise cleanup.
- `fast_clogit_plus_glmm_cpp`: same early return + blockwise cleanup.
- `fast_poisson_glmm_cpp`: early return only; Hessian rewrite reverted.

**Rejected:**
- Grouped Poisson Hessian rewrite: estimate-stable but not faster.
- Additional Weibull Hessian tightening: same pattern as Phase 3; preserved math but slower.

---

## Phase 5: Adjacent-Category Logit

### Planned

Cache `LinSpaced` category grids, remove per-row tiny temporaries, simplify Hessian cross-block updates.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `adjacent_category_logit` | direct objective path | 0.00010 | 0.00011 | **revert** |
| `adjacent_category_logit` | direct score path | 0.00010 | 0.00010 | **revert** |
| `adjacent_category_logit` | direct Hessian path | 0.00015 | 0.00020 | **revert** |
| `adjacent_category_logit` | fit path | 0.0075 | 0.0080 | **revert** |

**Estimate consistency verified:** `nll = 1392.8592205295745`, `score_sum = 45.46632076388376`, `h_trace = -7778.655883166635`, `fit_param_sum = -0.0542012540459799`. VGAM parity: `max_abs_beta_diff = 1e-16`, `max_abs_vcov_diag_diff = 1.57e-08`.

**All changes reverted.** The extra bookkeeping did not beat the compiler's handling of the simpler baseline. This kernel is small enough globally that it is not worth further effort unless its profiled share grows.

---

## Additional Non-GLMM Sweep

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `pair_distance_matrix` | direct distance build | 0.0080 | 0.0020 | **keep** |
| `d_optimal_search` | direct search | 0.0065 | 0.0055 | **keep** |
| `a_optimal_search` | direct search | 0.0130 | 0.0120 | **keep** |
| `zinb` | Hessian | 0.0010 | 0.0010 | **revert** |
| `zinb` | fit | 0.0835 | 0.1470 | **revert** |
| `zero_augmented_poisson` | Hessian | ~0 | ~0 | **revert** |
| `zero_augmented_poisson` | fit | 0.0010 | 0.0035 | **revert** |
| `negbin_with_var` | score/objective | 0.0010 | 0.0010 | **revert** |
| `negbin_with_var` | fit | 0.0070 | 0.0085 | **revert** |

**Pair-distance matrix** (`pair_dist_helpers.cpp`): replaced Rcpp element access inside the innermost loop with raw column-major pointers — 4x speedup. Numerical result unchanged: upper-triangle sum = `2087702.3513132515` before and after.

**D-optimal / A-optimal search** (`optimal_design_search.cpp`): cached matrix diagonals; switched to raw storage access for `P(i,j)` and `H(i,j)` in the nested scan. Modest but real wins; stochastic output shape and treatment-count constraint preserved.

**ZINB, ZAP, negbin**: Hessian/objective restructuring experiments were all estimate-stable but slower. Reverted in full. Root cause identified: the actual bottleneck is special-function evaluation (see v2 section below), not control flow or memory access patterns.

---

## v2 Re-Profiling Pass (2026-05-17–18)

Three open items from the initial report's "Next work" section.

### v2.1 — Ordinal CLMM: Re-Profile on Retained Tree

#### Timing (retained tree)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `ordinal_clmm` | estimate-only fit | 0.0180 |
| `ordinal_clmm` | full (with Hessian) | 0.0245 |

Phase 1 baseline was 0.0470s; retained tree at 0.0180s is consistent with the kept fixes.

#### Perf highlights (800 reps, full path)

| % samples | Symbol |
|---:|---|
| 35.03% | `GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()` |
| 21.89% | `__ieee754_exp_fma` |
| 12.06% | `__ieee754_log_fma` |
| 10.30% | `exp@@GLIBC_2.29` |
| 3.77% | Eigen GEMV |
| 2.40% | `glmm::log_sum_exp` |
| 1.21% | `malloc` |
| 1.15% | `cfree` |

Allocation pressure reduced (2.4% combined, down from previously dominant). Threshold-loop work no longer appears as a separate symbol. Dominant cost is now the generic engine loop (35%) and exp/log evaluation (44% combined).

#### Exp-reduction attempt — REVERTED

**Attempted changes:**

1. Threshold precomputation via `mutable` fields: added `mutable std::vector<double> alpha_` and `mutable std::vector<double> exp_diffs_` to `OrdinalLikelihood`, populated once per optimizer step in a `precompute(par)` const method called from `GLMMObjective::value()` and `operator()()`.
2. `pdf_from_cdf` optimization: added `static inline double pdf_from_cdf(double x, double F)` to all four link structs in `_glmm_links.h`. For `LogitLink`, returns `F*(1-F)` directly, avoiding the redundant `cdf(x)` call inside the old `pdf(x)` method.

Theoretical exp call reduction (K=4, 100 groups, 4 rows/group, 10 GH nodes): ~34,000 per `operator()()` → ~6,000. Parity vs `ordinal::clmm(nAGQ=20)`: max diff < 3e-6. Math was correct.

| Version | Median (s) |
|---|---:|
| Baseline (committed) | 0.021 |
| Full precompute + pdf_from_cdf | 0.38 |
| pdf_from_cdf only (precompute no-op) | 0.032 |

**Root cause of 18x regression:** The `mutable` keyword on heap-allocated `std::vector<double>` fields prevents the compiler from treating the class as alias-free across const method calls. The compiler cannot hoist `alpha_[i]` reads out of the GH quadrature inner loop (LICM disabled) and cannot auto-vectorize the suffix derivative loop because the mutable heap pointer may alias with other arguments. All 4000+ inner-loop evaluations per `operator()()` reload heap-allocated mutable vectors on every iteration.

**Decision: REVERT precompute machinery; RETAIN `pdf_from_cdf`.**

- `_glmm_links.h`: kept `pdf_from_cdf` on all four link structs.
- `fast_ordinal_clmm.cpp`: `OrdinalLikelihood` reverted to on-the-fly `get_alpha_bounds(par, ...)` with stack-local exp. No `mutable` fields.
- `_glmm_engine.h`: `model.precompute(par)` calls removed.

**Lesson:** NEVER use `mutable std::vector<double>` (or other mutable heap-allocated) fields in model classes used inside the GH quadrature inner loop of `GLMMObjective`. If threshold precomputation is revisited, use a caller-allocated stack buffer passed explicitly to `log_prob_derivs` — not `mutable` class state.

---

### v2.2 — Poisson GLMM: Re-Profile as an Independent Problem

#### Timing (retained tree)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `poisson_glmm` | estimate-only fit | 0.0030 |
| `poisson_glmm` | full (with Hessian) | 0.0035 |

Very small gap between estimate-only and full-path: Hessian is nearly free at this problem size; the optimizer iterations dominate.

#### Perf highlights (3000 reps, full path)

| % samples | Symbol |
|---:|---|
| 39.69% | `PoissonGLMMObjective::operator()` |
| (remainder) | GEMV and standard Eigen/libm math |

No Hessian symbol, allocation symbol, or special function appeared at notable sample weight. The objective loop is the entire story.

#### Decision

**No change.** Poisson GLMM is already well-structured for this problem scale. The reason the shared logistic/Weibull/clogit rewrite failed (Phase 2) is consistent: the existing loop is tight enough that a shared-pattern abstraction adds indirection without reducing work. Future work, if warranted, would require SIMD vectorization across observations within a group within a node — a much more invasive change, only justified if Poisson GLMM becomes a dominant bottleneck at larger scale.

---

### v2.3 — ZINB: Special-Function Profiling and Hoisting

#### Timing (retained tree, pre-hoisting)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `zinb` | estimate-only fit | 0.0230 |
| `zinb` | full (with Hessian) | 0.0220 |
| `zap` | estimate-only fit | 0.0010 |
| `zap` | full (with Hessian) | 0.0010 |

Estimate-only and full-path timings nearly identical — the bottleneck is not the Hessian accumulation step; it runs equally in both paths.

#### Perf highlights (800 reps, ZINB full path)

| % samples | Symbol |
|---:|---|
| 20.96% | `Rf_dpsifn.part.0.constprop.0` |
| 15.67% | `__ieee754_pow_fma` |
| 11.89% | `Rf_chebyshev_eval` |
| 11.34% | `__ieee754_log_fma` |
| 4.48% | `ZeroInflatedNegBin::hessian` |
| 4.03% | `ZeroInflatedNegBin::operator()` |
| 2.63% | `pow@@GLIBC_2.29` |
| 1.66% | `Rf_gammafn` |
| 1.64% | `malloc` |
| 1.55% | `__libc_malloc` |
| 1.43% | `cfree` |

More than 60% of samples are in special-function evaluation: `Rf_dpsifn` (digamma, 21%), `pow` (18%), `Rf_chebyshev_eval` (12%), `log` (11%), `gammafn` (1.7%). The `hessian` + `operator()` control-flow accounts for only ~8.5%.

This explains why the previous Hessian restructuring experiments were slower: they did not touch the dominant special-function cost and the restructuring added overhead.

#### Special-function hoisting — RETAINED

**Call count analysis (per `operator()()` call, n=100, p=7, 62% zeros, 38 non-zero obs, 29 distinct positive y values):**

| Call | Before | After | Savings |
|---|---:|---:|---:|
| `pow(theta/A, theta)` | 62 | 0 | 62 (replaced) |
| `exp(theta*(log_theta-log(A)))` | 0 | 62 | — |
| `log(A)` | 0 | 62 | — (but eliminates pow) |
| `digamma(theta)` | 38 | 1 | 37 |
| `lgamma(theta)` | 38 | 1 | 37 |
| `lgamma(y+1)` | 38 | 0 | 38 (precomputed at construction) |
| `digamma(y+theta)` | 38 | 29 | 9 |
| `lgamma(y+theta)` | 38 | 29 | 9 |
| redundant `log(theta/A)` | 38 | 0 | 38 (reuses `log_theta - log_A`) |

**Per `hessian()` call (additional):**

| Call | Before | After | Savings |
|---|---:|---:|---:|
| `digamma(theta)` | 38 | 1 | 37 |
| `trigamma(theta)` | 38 | 1 | 37 |
| `digamma(y+theta)` | 38 | 29 | 9 |
| `trigamma(y+theta)` | 38 | 29 | 9 |

**Implementation changes in `fast_zinb.cpp`:**

1. Constructor: added `m_y_slot[]` (obs → distinct-y table index), `m_distinct_y[]` (unique positive y values), `m_lgamma_y1[]` (`lgamma(y+1)` precomputed once at construction — constant across all optimizer calls).
2. `operator()()`: hoisted `digamma(theta)` and `lgamma(theta)` before the obs loop; built per-call `digamma_yptheta[]` and `lgamma_yptheta[]` tables over distinct y; replaced `pow(theta/A, theta)` with `exp(theta * (log_theta - log(A)))` — `log_theta` is the parameter itself (exact), `log(A)` precomputed once per obs and reused for both `log(theta/A)` and `log(mu/A)`.
3. `hessian()`: same hoisting of `digamma(theta)` and `trigamma(theta)`; per-call tables for `digamma(y+theta)` and `trigamma(y+theta)`.

**Benchmark result (n=100, 60 reps, estimate_only=TRUE):**

| Version | Median (s) | Min (s) |
|---|---:|---:|
| Baseline | 0.0040 | 0.0030 |
| Optimized | 0.0030 | 0.0020 |
| Speedup | **1.33x** | **1.50x** |

Parity: neg_loglik diff vs baseline = 4.5e-5 (PASS — within optimizer tolerance).

**Decision: RETAIN.**

---

### v2.4 — Weibull Frailty: Re-Profile on Retained Tree

#### Timing (retained tree, n=200, n_grp=50, p=4, K=20 GH nodes)

| Kernel | Path | Median (s) | Min (s) |
|---|---:|---:|---:|
| `weibull_frailty` | estimate-only fit | 0.0080 | 0.0070 |
| `weibull_frailty` | full (with Hessian) | 0.0080 | 0.0070 |

Estimate-only and full-path timing are nearly identical: the numerical Hessian (12 extra `operator()` calls at p+2=6 params) is marginal at this problem size, same pattern as Poisson GLMM v2.2.

#### Perf highlights (800 reps, estimate-only path)

| % samples | Symbol |
|---:|---|
| 41.11% | `WeibullFrailtyLikelihood::operator()` |
| 26.29% | `__ieee754_exp_fma` |
| 14.04% | `exp@@GLIBC_2.29` |
| 1.36% | `exp@plt` |
| 0.80% | `__ieee754_log_fma` |
| 0.94% | `malloc` + `cfree` |

**Combined exp cost: ~41.7%.** No GEMV, no special functions, no log_sum_exp visible as separate symbols. Full-path profile is structurally identical — Hessian adds no new dominant symbol.

#### Interpretation

Two distinct cost buckets:
1. **`operator()` control flow (41%)** — loop indexing, matrix element reads, branch overhead.
2. **Scalar `exp` (42%)** — all three `exp` variants sum to ~42%.

The initial sweep noted "vectorized Eigen packet math, Eigen GEMV, vectorized exponentials." On the retained tree, the picture is cleaner: GEMV is absorbed into the operator() bucket (too fast to appear separately at this sample rate), and no special functions appear at all.

#### Exp call count (per `operator()` call, n=200, n_grp=50, K=20)

| Source | Count |
|---|---:|
| Inner forward loop: `exp(wik)` per (obs, node) | n×K = 4,000 |
| `log_sum_exp_wf` per group: K exps × G groups | G×K = 1,000 |
| Posterior weights pass: K exps × G groups | G×K = 1,000 |
| **Total** | **6,000** |

#### Actionable next step: multiplicative exp decomposition

Since `wik = (log_y[i] - eta[i] - uk) / sigma_eps = base_wi[i] + delta_k`, where  
`base_wi[i] = (log_y[i] - eta[i]) / sigma_eps` (obs-only, computed once per optimizer step)  
`delta_k = -uk / sigma_eps` (node-only, computed once per node),

the inner loop's `exp(wik) = exp(base_wi[i]) * exp(delta_k)` — a multiply, not an exp call.

This reduces the inner forward loop from n×K=4,000 exp calls to n+K=220 exp calls, with the remaining n+K ops being multiplies. Applied correctly, this should reduce total exp calls from 6,000 to ~2,220 (37% of current), eliminating the dominant cost.

A separate LSE/posterior-weights fusion (analogous to `log_sum_exp_and_weights` in `_glmm_engine.h`) would eliminate the second `G×K` exp pass, saving another 1,000 calls.

**Why Phase 3 failed:** The previous Phase 3 math rewrite benchmarked 2.9× slower without a clear root cause. The current profile gives better guidance: Phase 3 likely broke the compiler's vectorized exp batch (`__ieee754_exp_fma` uses AVX2 batched evaluation), replacing it with piecemeal scalar multiplies or changing memory access patterns. The multiplicative decomposition must be implemented so the `exp(base_wi)` array pass and `exp(delta_k)` scalar are each visible to the vectorizer as contiguous runs — not interleaved with control flow.

**Guard:** Benchmark immediately after each change against the retained baseline (0.0080s median). Revert if any sub-step is slower.

---

### v2.5 — Gaussian LMM: Profile at Relevant Scale

#### Timing (n=2000, n_grp=1000 all-pairs, p=4)

| Kernel | Path | Per-call (µs) | Ratio |
|---|---:|---:|---:|
| `gaussian_lmm` | estimate-only | 335 | 1× |
| `gaussian_lmm` | full (with Hessian) | 2,570 | 7.7× |

Scaling reference across problem sizes:

| n | estimate-only (µs/call) | full path (µs/call) |
|---:|---:|---:|
| 500 | 85 | 620 |
| 2,000 | 335 | 2,435 |
| 5,000 | 625 | 6,410 |
| 10,000 | 1,270 | 13,135 |

Approximately linear in n for both paths.

#### Perf highlights — estimate-only path (1000 reps, n=2000)

| % samples | Symbol |
|---:|---|
| 11.16% | `R_gc_internal` |
| 9.56% | `neg_ll_and_grad` |
| 6.33% | `bcEval_loop` |
| 5.12% | `fast_gaussian_lmm_cpp` (wrapper) |
| 4.65% | `Rf_allocVector3` |

**R overhead dominates.** The actual C++ `neg_ll_and_grad` is only ~10% of total wall time; the rest is R GC, R bytecode interpreter, and Rcpp wrapper overhead. This confirms that the estimate-only path has no meaningful C++ optimization target — each call is too fast (85 µs at n=500, 335 µs at n=2000) for the arithmetic to be the bottleneck.

#### Perf highlights — full path / Hessian (500 reps, n=2000)

| % samples | Symbol |
|---:|---|
| 12.65% | `lmm_analytic_hessian` |
| 7.95% | `Eigen::PlainObjectBase::resize` |
| 6.39% | `R_gc_internal` |
| 6.13% | `cfree` |
| 5.63% | `Eigen::generic_product_impl` (MatMul) |
| 4.00% | `bcEval_loop` |
| 3.74% | `Eigen::product_evaluator` |
| 3.23% | `malloc` |

**Combined allocation overhead: ~21%** (`resize` 7.95% + `cfree` 6.13% + `malloc` 3.23% + `product_evaluator` 3.74%). The Hessian function allocates approximately 15 tiny `Eigen::MatrixXd` objects per group (Identity, Ones, V, LDLT result P, dV\_e, dV\_b, d²V\_ee, d²V\_bb, dP\_e, dP\_b, and intermediate products in the 2×2 (a,b) loop). At G=1000 groups, this is ~15,000 heap allocs per Hessian call.

#### Interpretation

The two paths have distinct bottleneck structures:

1. **Estimate-only**: dominated by R overhead (GC, bytecode), not C++ arithmetic. The `neg_ll_and_grad` loop is clean — one GEMV + a scalar per-group accumulation with no exp calls, no special functions, no allocation. No C++ optimization available.

2. **Full path (Hessian)**: ~21% of cycles burned in heap allocation/deallocation. Root cause: `lmm_analytic_hessian` uses `Eigen::MatrixXd` (dynamic allocation) for 1×1 and 2×2 group matrices. Since KK designs only have m=1 (singleton) and m=2 (pair) groups, all group-level linear algebra has closed-form scalar/2×2 expressions — no dynamic allocation needed.

#### Actionable optimization

Replace the per-group `Eigen::MatrixXd` Hessian loop with closed-form scalar expressions for m=1 and m=2:
- **m=1**: V = v\_e + v\_b (scalar), V⁻¹ = 1/(v\_e+v\_b), all derivatives are scalars.
- **m=2**: V = [[v\_e+v\_b, v\_b],[v\_b, v\_e+v\_b]], analytically invertible (det = v\_e(v\_e+2v\_b)). All P, dP, and matrix traces reduce to O(1) arithmetic.

Expected payoff: eliminate 15,000 heap allocs/call → cut full-path time from 2,570 µs to ~2,000 µs (~20%). At the sweep level, gaussian\_lmm\_full is 2.6% of total — the net gain is ~0.5% of sweep time. **Not implemented at this time.**

#### Decision

**No changes.** The estimate-only path is R-overhead-bound; the full-path Hessian optimization would save ~20% of a 2.6%-share kernel (~0.5% overall). The allocation fix is technically clean but globally marginal given current priorities.

---

## Current Retained State

| File | Retained changes |
|---|---|
| `EDI/src/fast_ordinal_clmm.cpp` | `get_alpha_bounds()` replaces old full-threshold reconstruction per call; no tiny `dp` allocation per row |
| `EDI/src/_glmm_engine.h` | Workspace reuse for inner derivative buffers; no `model.precompute()` calls |
| `EDI/src/_glmm_links.h` | `pdf_from_cdf` static method on all four link structs |
| `EDI/src/fast_logistic_glmm.cpp` | Direct group accumulation; persistent scratch buffers; Phase 3 math tightening; Hessian blockwise cleanup; `estimate_only` early return |
| `EDI/src/fast_weibull_frailty.cpp` | Phase 1/2 direct-accumulation changes; persistent scratch buffers; multiplicative exp decomposition `exp(wik)=exp(base_wi[i])×exp(delta_k)` + fused LSE+posterior-weights pass (Phase 3 math tightening reverted) |
| `EDI/src/fast_clogit_plus_glmm.cpp` | Direct accumulation; persistent scratch buffers; Phase 3 tightening; Hessian blockwise cleanup; `estimate_only` early return |
| `EDI/src/fast_poisson_glmm.cpp` | `estimate_only` early return only (shared-pattern and Hessian rewrites reverted) |
| `EDI/src/optimal_design_search.cpp` | Cached matrix diagonals; raw storage access in swap-scan loops |
| `EDI/src/pair_dist_helpers.cpp` | Raw column-major pointer rewrite replacing Rcpp accessor overhead |
| `EDI/src/fast_zinb.cpp` | Distinct-y construction precompute; digamma/lgamma/trigamma hoisting; pow→exp(log) replacement |
| `EDI/src/fast_negbin_regression.cpp` | Distinct-y tables for lgamma/digamma(y+θ); hoist lgamma(θ), digamma(θ), log(θ); explicit NB log-PMF replacing R::dnbinom_mu; fast_digamma; preallocated table vectors |
| `EDI/src/fast_hurdle_negbin.cpp` | TruncatedNegBinCount: distinct-y tables (lgamma/digamma/trigamma(y+θ)); hoist lgamma(θ), digamma(θ), log(θ) per step; explicit NB log-PMF + truncation correction; fast_digamma throughout; preallocated table vectors |
| `EDI/src/fast_beta_regression.cpp` | DigammaFunctor: R::digamma → fast_digamma; R::lgammafn → std::lgamma in all unaryExpr lambdas and scalar calls; hoist digamma(φ) per operator() step; hoist trigamma(φ) before expected_hessian inner loop |
| `EDI/src/fast_zero_augmented_poisson.cpp` | operator(): (1) per-row rank-1 grad updates → scalar weight vectors + GEMV; (2) scalar exp-reuse: compute `e_neg = exp(-eta_z)` once per obs for both sigmoid `p = 1/(1+e_neg)` and softplus `lse = eta_z + log1p(e_neg)` (eta_z>0) / `log1p(1/e_neg)` (eta_z≤0), eliminating one `exp` call per positive observation |
| `EDI/src/fast_zinb.cpp` | operator(): same GEMV refactor; additionally moved phi=pow(...) inside yi<=0 branch (avoids pow for ~60% of obs); replaced `log(theta/(theta+mu))` with `log_theta - log_denom` for positive obs (saves 1 log call per positive obs) |
| `EDI/src/fast_log_binomial_regression.cpp` | IRLS backtracking refactor: (1) `ll_curr` cached across iterations (avoids one `loglik_constrained_binomial` call per accepted IRLS step); (2) precompute `delta_eta = X_free * direction` per iteration so each backtracking probe is O(n) vector-add + `loglik_from_eta` (no GEMV); (3) preallocate `eta_try(n)` buffer outside while-loop to avoid per-halving heap allocation; added `loglik_from_eta` helper taking precomputed eta directly |
| `EDI/src/_helper_functions.h` | Shared contiguous-group layout helper introduced in Phase 2 |

**Unchanged from baseline:** `fast_adjacent_category_logit.cpp`, `fast_gaussian_lmm.cpp`, `fast_coxph_regression.cpp`.

---

## Summary Benchmark Table (All Phases)

| Kernel | Path | Original baseline (s) | Final retained (s) | Net change |
|---|---:|---:|---:|---|
| `ordinal_clmm` | estimate-only fit | 0.0470 | 0.0180 | −62% |
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0060 | −87% |
| `logistic_glmm` | Hessian | 0.0070 | 0.0010 | −86% |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | −14% |
| `clogit_plus_glmm` | estimate-only fit | 0.0185 | 0.0200 | +8% (estimate-only now skips Hessian but different problem scale) |
| `clogit_plus_glmm` | Hessian | 0.01015 | 0.00515 | −49% |
| `poisson_glmm` | estimate-only fit | 0.0115 | 0.0010 | −91% (early return + smaller scale) |
| `pair_distance_matrix` | direct build | 0.0080 | 0.0020 | −75% |
| `d_optimal_search` | direct search | 0.0065 | 0.0055 | −15% |
| `a_optimal_search` | direct search | 0.0130 | 0.0120 | −8% |
| `zinb` | estimate-only fit | 0.0040 | 0.0030 | −25% (1.33x) |
| `weibull_frailty` | estimate-only fit (v2.4+exp-reduction) | 0.0080 | 0.0050 | −37% |
| `hurdle_negbin` (TruncatedNegBinCount) | estimate_only (benchmark n=1000) | 13.63ms | 1.85ms | **−86%** (7.4x); distinct-y lgamma/digamma tables + explicit NB log-PMF + fast_digamma |
| `beta_regression` | estimate_only (benchmark n=1000) | 7.93ms | 1.77ms | **−78%** (4.5x); std::lgamma + fast_digamma in DigammaFunctor + hoist scalars |
| `zero_inflated_poisson` (ZIP) | estimate_only (micro n=1000) | 4.94ms | 0.99ms | **−80%** (5x); per-row rank-1 grad updates → weight vectors + GEMV |
| `zero_inflated_poisson` (ZIP) | estimate_only (profiling n=1000, seed=99) | 5.96ms (post-GEMV committed) | 4.97ms | **−17%**; scalar exp-reuse: single `exp(-eta_z)` serves both sigmoid and softplus (TODO-18); reverted erroneous Eigen precompute that added ~28% more exp calls under `EIGEN_DONT_VECTORIZE` |
| `zinb` | estimate_only (micro n=1000) | 2.73ms | 1.15ms | **−58%** (2.4x); same GEMV + hoist phi inside zero branch + log(θ/(θ+μ))→log_theta−log_denom |
| `log_binomial` (IRLS) | estimate (benchmark n=1000) | 9.40ms | 2.16ms | **−77%** (4.4x); ll_curr caching across IRLS iters + precompute delta_eta for O(n) backtracking probes (eliminates GEMV per backtracking halving) |

---

## Test Coverage

| Test file | Coverage |
|---|---|
| `EDI/tests/testthat/test-glmm-cpp-equivalence.R` | Logistic GLMM and Poisson GLMM equivalence vs `lme4::glmer` |
| `EDI/tests/testthat/test-weibull-frailty.R` | Weibull frailty inference path |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:55` | `fast_neg_bin_with_var_cpp` vs `MASS::glm.nb` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:289` | Adjacent-category logit vs `VGAM::vglm(acat(parallel=TRUE))` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:329` | `fast_zinb_cpp` vs `glmmTMB` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:348` | `fast_zero_augmented_poisson_cpp` vs `glmmTMB` |

**Gaps:** No canonical-package `testthat` coverage for `clogit_plus_glmm`, d-optimal search, or pair-distance helper.

---

## TODO: Future Optimization Work

### High priority

**TODO-1: Ordinal CLMM — threshold precomputation via caller-allocated stack buffer** ✓ DONE
File: `EDI/src/fast_ordinal_clmm.cpp`, `EDI/src/_glmm_engine.h`
Implemented: `GLMMObjective::value()` and `operator()()` now call `model.fill_alpha(par, alpha_buf.data())` once per optimizer step (allocating a `std::vector<double>(nm)` outside the inner loops), then pass `alpha_buf.data()` into `log_prob` and `log_prob_derivs`. `get_alpha_bounds` was rewritten to O(1) array lookup (zero exp calls per observation). The `dp` loop in `log_prob_derivs` was also updated to use alpha diffs (`alpha[j] - alpha[j-1]`) instead of `exp(par[j])`. All four link functions (logit, probit, cauchit, cloglog) converge correctly. Re-profile with `perf` on retained tree to measure actual speedup vs. the ~5× theoretical estimate.

**TODO-2: ZINB — faster digamma approximation** ✓ DONE
File: `EDI/src/fast_zinb.cpp`
Implemented: `fast_digamma(x)` — uses recurrence `ψ(x) = ψ(x+1) - 1/x` to shift x ≥ 8, then applies A&S 6.3.18 asymptotic expansion as a 5-term Horner-form polynomial. Accurate to ≤ 4e-12 relative error across x ∈ [0.1, 1000]; falls back to `R::digamma` for x ≤ 0. Replaces both `R::digamma(theta)` and `R::digamma(y + theta)` calls in `operator()`. Re-profile with `perf` to confirm expected 3–5× reduction in digamma time.

**TODO-3: D-optimal / A-optimal search — algorithmic pruning** ✓ DONE
File: `EDI/src/optimal_design_search.cpp`
Implemented sorted-candidate pruning for both D-optimal and A-optimal:
- **D-optimal**: precompute `max_P = P.cwiseAbs().maxCoeff()` once. Per while-loop iteration, sort `t_idxs` ascending by `A[i] = -2*Pw[i]+p_diag[i]` and `c_idxs` ascending by `B[j] = 2*Pw[j]+p_diag[j]`. Inner-j break: `A[i]+B[j] >= best_delta + 2*max_P` (since `delta >= A[i]+B[j]-2*max_P`). Outer-i break: `A[i]+B_min >= best_delta + 2*max_P` where `B_min = B[c_idxs[0]]`.
- **A-optimal**: same pattern with combined score `C[i] = A_H[i] + obj_curr*A_P[i]`, `D[j] = B_H[j] + obj_curr*B_P[j]`, and prune threshold `2*(max_H + obj_curr*max_P)`.
Treatment counts validated (exactly n_T per simulation). All 500-sim quality distributions match reference. Re-profile with `perf` to measure actual speedup.

### Medium priority

**TODO-4: Weibull frailty — targeted exp-reduction profiling** ✓ DONE
File: `EDI/src/fast_weibull_frailty.cpp`
Fresh `perf` trace gathered on retained Phase 1/2 tree (v2.4 section above). Bottleneck: 42% scalar `exp` + 41% `operator()` control flow; no GEMV, no special functions. Exp count: 6,000 per `operator()` call (4,000 inner loop + 2,000 LSE/weights). Identified actionable next step: multiplicative decomposition `exp(wik) = exp(base_wi[i]) * exp(delta_k)` reduces inner loop from n×K=4,000 to n+K=220 exp calls. Phase 3 likely broke vectorizer batch eval; new approach separates the n-length and K-length exp arrays cleanly to preserve SIMD batching.

**TODO-6: Negbin regression — profile the actual bottleneck** ✓ DONE
File: `EDI/src/fast_negbin_regression.cpp`
Applied the same hoisting/tabulation strategy as ZINB (TODO-2):
- Hoist `digamma(theta)`, `lgamma(theta)`, `log(theta)` out of the observation loop (was repeated n times per optimizer step).
- Build per-distinct-y tables for `lgamma(y+theta)` and `digamma(y+theta)` at construction (like ZINB).
- Replace `R::dnbinom_mu(yi, theta, mu_i, true)` with explicit formula using table lookups — avoids slow R function dispatch.
- Use `fast_digamma` (moved to `_helper_functions.h` so both ZINB and negbin share it) for all digamma calls.
- Added `noalias()` to score_beta accumulation.
Hessian diagonal matches finite-difference to 6e-8; fits converge correctly. Re-profile with `perf` to confirm `R::digamma` and `R::dnbinom_mu` are no longer in the top symbols.

### Lower priority

**TODO-8: Ordinal CLMM — reduce log_sum_exp calls** ✓ DONE
File: `EDI/src/_glmm_engine.h`
Implemented `log_sum_exp_and_weights(x, weights)`: computes lse(x) and fills `weights[k] = exp(x[k] - lse)` in a single pass. Used in `GLMMObjective::operator()` in place of `log_sum_exp` + per-k `exp(log_terms[k] - ll_g)` recomputation. Saves nn exp calls per group per optimizer step (nn=20 nodes, previously computed twice: once inside log_sum_exp, once in the gradient loop).

**TODO-9: ZINB/ZAP — allocator pressure** ✓ DONE
File: `EDI/src/fast_zinb.cpp`, `EDI/src/fast_negbin_regression.cpp`
Added `m_lgamma_yptheta` and `m_digamma_yptheta` as preallocated member vectors (sized at construction when the distinct-y count is known) to `ZeroInflatedNegBin` and `NBLogLik`. The per-call `std::vector<double>(nd)` allocations in `operator()` are replaced with in-place fills of the preallocated buffers — eliminates 2 heap allocations per optimizer step for each model.

**TODO-10: Full re-benchmark after all retained changes** ✓ DONE
Tooling: `benchmark/benchmark_model_fits.R` → `package_metadata/benchmark_model_fits.md`
Re-run 2026-06-30 (latest: 21:48 JST). 73 paths benchmarked; no regressions. `InferenceContinQuantileRegr` shows ~1.0x (p>0.05, not significant, untouched code). HurdleNegBin: 13.63ms → 1.85ms (−86%). BetaRegr: 7.93ms → 1.77ms (−78%). All other optimized kernels at or above prior baseline.
Full re-run 2026-07-01 (includes TODO-15/16 GEMV fixes and TODO-17 log-binomial backtracking): LogBinomial estimate 9.40ms → 2.16ms (−77%); ZIP estimate 5.07ms → 4.01ms (−21%); ZINB estimate 2.11ms → 1.69ms (−20%); HurdlePoisson estimate 1.93ms → 1.59ms (−18%); HurdleNegBin estimate → 1.72ms (−87% vs 13.63ms); BetaRegr estimate → 1.56ms (−80% vs 7.93ms).

**TODO-11: Gaussian LMM — profile at relevant scale** ✓ DONE
File: `EDI/src/fast_gaussian_lmm.cpp`
Profile gathered (v2.5 section above). Findings: estimate-only path is R-overhead-bound (neg_ll_and_grad = 9.6% of wall time at n=2000, 335 µs/call); full-path Hessian burns ~21% on heap allocation from ~15,000 tiny Eigen::MatrixXd allocs/call (G=1000 groups × ~15 allocs). Closed-form scalar Hessian for m=1 and m=2 would cut full-path by ~20%, but gaussian_lmm_full is only 2.6% of sweep total — net gain ~0.5%. No changes made.

**TODO-15: ZAP/ZIP — replace per-row rank-1 grad updates with weight vectors + GEMV** ✓ DONE
File: `EDI/src/fast_zero_augmented_poisson.cpp`
Replaced all `m_X.row(i).transpose() * scalar` and `m_Xzi.row(i).transpose() * scalar` rank-1 gradient updates in `operator()` with scalar weight accumulation into `w_cond[n]` and `w_zi[n]` vectors, followed by single `m_X.transpose() * w_cond` and `m_Xzi.transpose() * w_zi` GEMV calls after the observation loop. Applies to both ZIP and hurdle branches. Root cause: row access on column-major matrix is stride-n (non-contiguous); 2×n scattered writes per optimizer step were preventing BLAS vectorization. Micro-benchmark (n=1000): 4.94ms → 0.99ms (−80%, 5x). Correctness verified on ZIP and hurdle data.

**TODO-16: ZINB — GEMV + phi hoisting + log_denom reuse** ✓ DONE
File: `EDI/src/fast_zinb.cpp`
Same per-row rank-1 → GEMV refactor for `m_Xc.row(i).transpose()` and `m_Xz.row(i).transpose()` in `operator()`. Additionally: (1) moved `phi = pow(theta/(theta+mu), theta)` inside `yi<=0.0` branch — avoids ~60% of pow calls (positive obs never use phi); (2) replaced `phi = pow(...)` with `phi = exp(theta*(log_theta - log_denom))` to reuse the precomputed `log_denom`; (3) for positive obs, replaced `theta*log(theta/(theta+mu)) + yi*log(mu/(theta+mu))` with `theta*(log_theta - log_denom) + yi*(ec - log_denom)` — saves 2 log calls per positive obs. Micro-benchmark (n=1000): 2.73ms → 1.15ms (−58%, 2.4x). Correctness verified.

**TODO-13: HurdleNegBin — special-function hoisting in TruncatedNegBinCount** ✓ DONE
File: `EDI/src/fast_hurdle_negbin.cpp`
Applied the same distinct-y precomputation + fast_digamma pattern as negbin (TODO-6) to `TruncatedNegBinCount::operator()` and `hessian()`: added `m_y_slot[]`, `m_distinct_y[]`, `m_lgamma_y1[]`, `m_lgamma_yptheta[]`, `m_digamma_yptheta[]`, `m_trigamma_yptheta[]` members; constructor builds distinct-y table via `std::unordered_map`; operator() hoists `lgamma_r`, `digamma_r`, `log_r`, fills lgamma/digamma tables, replaces `R::dnbinom_mu` with explicit NB log-PMF (`m_lgamma_yptheta[slot] - lgamma_r - m_lgamma_y1[slot] + r*(log_r - log_denom) + y*(log_mu_i - log_denom) - log(trunc_denom)`), replaces `R::digamma(y+r) - R::digamma(r)` with table lookup; hessian() hoists digamma_r, trigamma_r, log_r, fills digamma/trigamma tables. Benchmark (n=1000): 13.63ms → 1.85ms (−86%, 7.4x speedup). Correctness verified: coefficients within noise of known truth.

**TODO-14: BetaRegr — std::lgamma + fast_digamma + scalar hoisting** ✓ DONE
File: `EDI/src/fast_beta_regression.cpp`
Three targeted changes: (1) `DigammaFunctor::operator()`: `R::digamma(x)` → `fast_digamma(x)` — applies to all per-element digamma calls in `operator()` and `hessian()` unaryExpr passes; (2) all `R::lgammafn(x)` in unaryExpr lambdas and scalar use → `std::lgamma(x)` (avoids R error-handling dispatch overhead per element); (3) hoist `const double digamma_phi = fast_digamma(phi)` before `d_neg_ll_d_phi` computation in `operator()`, and `const double trigamma_phi = R::trigamma(phi)` before inner loop in `expected_hessian()` (was recomputed per observation). Benchmark (n=1000): 7.93ms → 1.77ms (−78%, 4.5x speedup). Correctness verified: coefficients within noise of known truth.

**TODO-17: LogBinomial — IRLS backtracking refactor** ✓ DONE
File: `EDI/src/fast_log_binomial_regression.cpp`
Three changes applied to `fit_constrained_binomial_cpp_impl` (used by both log-link and identity-link paths):
(1) **Cache `ll_curr` across iterations**: initialize once before IRLS loop; carry forward `ll_curr = ll_new` after each accepted step — eliminates one `loglik_constrained_binomial(X, y, beta)` call per accepted IRLS iteration (each call does an O(n*p) GEMV + O(n) exp/log1p).
(2) **Precompute `delta_eta = X_free * direction` once per IRLS iteration**: each backtracking probe becomes O(n) vector-add `eta_try = eta + step * delta_eta` + O(n) `loglik_from_eta` scalar loop (no GEMV). Previously: each backtracking `loglik_constrained_binomial(X, y, beta_new)` recomputed X*beta_new as a full GEMV.
(3) **Preallocate `eta_try(n)` before while-loop**: avoids one heap allocation per backtracking halving.
Added `loglik_from_eta` helper: evaluates the log-likelihood directly from a precomputed eta vector (identical scalar loop to `loglik_constrained_binomial` but skips the X*beta GEMV).
perf profile before (n=1000, 100 IRLS iters): 23.46% `log1p`, 19.19% XtWX GEMM, 12.38% `loglik_constrained_binomial` overhead. After: log1p reduced to 17.19% (caching eliminates one loglik pass per accepted step). Benchmark (n=1000): 9.40ms → 2.16ms (−77%, 4.4x speedup). Correctness verified vs `glm(family=binomial(link="log"))`.

**TODO-18: ZIP/Hurdle — scalar exp-reuse for sigmoid + softplus** ✓ DONE
File: `EDI/src/fast_zero_augmented_poisson.cpp`
Root cause: `ZeroAugmentedPoisson::operator()` computed `exp(-eta_z)` for sigmoid, then `log1pexp(eta_zi[i])` called `exp()` a second time internally for the same observation (for positive observations). Since `log1pexp(x) = x + log1p(exp(-x))` for x>0, the `exp(-x)` used there equals the `e_neg` already computed for sigmoid.
Fix: compute `e_neg = std::exp(-eta_z)` once per observation; derive both sigmoid (`p = 1/(1+e_neg)`) and softplus (`lse = eta_z + std::log1p(e_neg)` for eta_z>0; `lse = std::log1p(1.0/e_neg)` for eta_z≤0 — one division instead of one exp). Eliminates one exp() call per positive observation (~68% of obs) at the cost of one division for ~50% of those. Pure scalar C++, no SIMD dependency; portable to all platforms and safe with `EIGEN_DONT_VECTORIZE`.
Also **reverted** an earlier attempt that used Eigen `array().min(700.0).exp()` for bulk precompute: that approach (a) segfaults WITHOUT `EIGEN_DONT_VECTORIZE` due to R's 8-byte-aligned heap vs AVX2's 32-byte alignment requirement; (b) with `EIGEN_DONT_VECTORIZE` set (which prevents the crash), Eigen's `array().exp()` degrades to scalar glibc exp per element — adding ~28% more exp calls (exp(-lambda) for ALL n obs instead of just the ~32% zeros). Note: `EIGEN_DONT_VECTORIZE` was added in commit `972b77e8` precisely to prevent these alignment segfaults; it provides no performance benefit and prevents Eigen's internal pexp SIMD kernel from activating.
Also **reverted** a batch-exp experiment using MKL's `vdExp` via `dlsym(RTLD_DEFAULT, "vdExp")` (available in R processes with MKL as BLAS backend; falls back to scalar on any other platform). Despite MKL vdExp being ~4× faster than scalar exp for large contiguous arrays, it benchmarked at 7.47ms — worse than the 4.97ms scalar baseline. Root cause: for n=1000, the overhead of 5 separate O(n) memory passes (prep loop + 2 batch_exp calls + main loop reading 2 arrays) exceeds the SIMD benefit; memory bandwidth is cheap at this scale and scalar exp is the bottleneck in isolation, not the loop. Batch-exp would likely win for n≥10000 where exp arithmetic outweighs memory bandwidth.
Benchmark (profiling data, n=1000, seed=99): 5.96ms (committed post-GEMV Eigen-precompute version) → 4.97ms (−17%).

**TODO-12: Clogit-plus-GLMM — add canonical-package testthat coverage** ✓ DONE (pre-existing)
File: `EDI/tests/testthat/test-clogit-plus-glmm-cpp-equivalence.R`
File was added in commit 4fea4627. Covers: concordant-only vs `lme4::glmer` (nAGQ=20), discordant-only vs `survival::clogit` (with and without covariates), combined model convergence + score-at-zero fingerprint, and structural smoke tests.

---

## July 2026 Sweep — New Optimization TODO List

Derived from parallel perf-annotate + perf-stat analysis of all 57 annotated kernels (July 2026). Priorities: **CRITICAL** > **HIGH** > **MEDIUM** > **LOW**.

### CRITICAL priority

**TODO-19: ZINB — preallocate `m_eta_c`, `m_eta_z`, `m_w_c`, `m_w_z` as class members** ✓ DONE
File: `EDI/src/fast_zinb.cpp:35–148`
`ZeroAugmentedPoisson` received this fix (TODO-15); `ZeroInflatedNegBin` was missed. `operator()` allocates 4 n-element `VectorXd` per optimizer call (2 GEMVs capturing into `const VectorXd` + 2 `::Zero(m_n)` weight vectors). Profile: `rep stosb` memset = 85.88% of `zinb_est` samples; PLT `operator new` stub = 53.17% local. Fix: add `m_eta_c(y.size()), m_eta_z(y.size()), m_w_c(y.size()), m_w_z(y.size())` as class members initialized in constructor; use `.noalias() =` for GEMVs and `.setZero()` for weight init inside `operator()`. Safe: `operator()` is non-`const`, no LICM/alias risk. Expected: 3–5× speedup on `zinb_est`.

**TODO-20: HurdlePoisson GLMM hessian — preallocate G×K working buffers + fix Xg matrix copy** ✓ DONE
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:141,154,234,244,275–304`
`hurdle_p_var` takes 4× longer than `hurdle_p_est` (66.6s vs 16.3s) — hessian dominates ~75% of var-path time. Root cause: per-call heap allocs of `E_Hik`, `E_GiGiT` (MatrixXd(total,total)), `G_avg` (VectorXd(total)), plus G×K allocs of `res_k`, `d2e_k`, `G_ik`, `H_ik`. For G=200, K=7: ~6,200 allocs per hessian call = 21% of wall time in `malloc`/`cfree`/`calloc`. Add preallocated member fields to `HurdlePoissonGLMMObjective`: `m_exp_bvals(K)`, `m_eta0(max_grp_sz)`, `m_res_k(max_grp_sz)`, `m_d2e_k(max_grp_sz)`, `m_G_ik(total)`, `m_H_ik(total,total)`, `m_E_Hik(total,total)`, `m_E_GiGiT(total,total)`, `m_G_avg(total)`. In `hessian()` use `.setZero()` per-group. Also fix: `const Eigen::MatrixXd Xg = dat.X_s.middleRows(gs, sz)` (line 244) copies the submatrix G times per call — change to `const auto Xg = dat.X_s.middleRows(gs, sz)` (zero-copy expression view). Apply same fix in `operator()` line 154 for `eta0`. Also move `std::vector<double> exp_bvals(K)` (lines 141, 234) to `m_exp_bvals` member.
Implemented with preallocated `m_b_vals`, `m_log_terms`, `m_eta0`, per-node Hessian buffers, per-group expectation buffers, and in-place small-p crossproduct helpers to avoid hidden `weighted_crossprod()` / `Xg.transpose() * d2e_k` temporaries. Correctness: `test-glmm-cpp-equivalence.R` passed (11/11). Targeted full-inference GLMM benchmark using the `benchmark_model_fits.R` adaptive timing setup with the profiler's `fast_hurdle_poisson_glmm_cpp(..., estimate_only=FALSE)` expression: 1.936842ms → 1.905263ms median (N_WALD=200; modest 1.6% improvement at this scale).

### HIGH priority

**TODO-21: Logistic + Probit IRLS — replace manual XtWX triple-loop with cache-friendly score + weighted crossproduct** ✓ DONE
Files: `EDI/src/fast_logistic_regression.cpp:154–177`, `EDI/src/fast_probit_regression.cpp:170–184`
The manual `for(i) for(j) for(k)` XtWX accumulation uses stride-n column reads (column-major X, j-loop steps across rows) which prevents vectorization and accounts for ~81% of `logistic_est` and ~70% of `probit_est` EDI.so samples. Replace with: (1) build weighted residual `diff[i] = w[i]*(y[i]-mu[i])`, (2) `score_free.noalias() = X_free.transpose() * diff` (GEMV, already AVX2), (3) `XtWX = weighted_crossprod(X_free, w)` (BLAS DGEMM path). Apply to both `use_weights` and non-`use_weights` branches.
Implemented a fused column-contiguous score + weighted-crossproduct helper in both logistic and probit IRLS. The direct generic `weighted_crossprod(X_free, w)` version benchmarked slower at the benchmark-model-fits scale, so the final implementation preserves one fused score/Fisher pass while avoiding row-wise stride-n access. Correctness: `test-logistic-cleanup.R`, `test-fast_glm_outputs.R`, `test-warm-start-weights.R`, `test-incidence-probit.R`, and `test-rcpp-fitting-equivalence.R` passed. Targeted estimate-only benchmark (N_GLM=1000, same adaptive timing setup as `benchmark_model_fits.R`): logistic_est 0.165482ms → 0.162029ms (−2.1%); probit_est 0.591988ms → 0.597938ms (+1.0%, effectively flat/slightly slower at p=6).

**TODO-22: Add `fast_lgamma` to `_helper_functions.h`** ✓ DONE
File: `EDI/src/_helper_functions.h`
After TODO-14 replaced `R::digamma` with `fast_digamma` in BetaRegression, `__lgamma_r_finite` now dominates `beta_var` (~40–50% of samples) and is major in `beta_est`. Pattern mirrors `fast_digamma`: Stirling asymptotic expansion `lgamma(x) ≈ (x-0.5)*log(x) - x + 0.5*log(2π) + 1/(12x) - 1/(360x³) + ...` for x ≥ 8; recurrence `lgamma(x) = lgamma(x+1) - log(x)` to shift up; rational polynomial for moderate x. Error ≤ 1e-10 relative; fallback to `std::lgamma` for x ≤ 0. Replace all `std::lgamma(x)` in `BetaRegression::operator()` (lines 62–65). Shared by BetaRegression, ZOIB, and any future kernel.
Implemented `fast_lgamma` using a Lanczos rational approximation for moderate positive inputs and a Stirling expansion for x ≥ 8; nonpositive/nonfinite inputs fall back to `std::lgamma`. Replaced the three `std::lgamma` calls in `BetaRegression::operator()`. Correctness: beta-related coverage in `test-fast_glm_outputs.R`, `test-argument-permutations.R`, `test-rcpp-fitting-equivalence.R`, and `test-rcpp-fitting-real-data.R` passed after a clean rebuild. Targeted adaptive benchmarks: `beta_est` 1.486301ms → 1.191781ms (−19.8%); `beta_var` 0.405680ms → 0.346220ms (−14.7%).

**TODO-23: BetaRegression — fuse 4 array passes into single scalar loop** ✓ DONE
File: `EDI/src/fast_beta_regression.cpp:60–88`
Current: 4 separate `unaryExpr` passes over n elements (2 lgamma + 2 digamma), materializing 4 intermediate `VectorXd` temporaries with heap allocation. Replace with a single scalar obs loop accumulating `neg_ll` and `w_grad[i]` simultaneously using `fast_lgamma` + `fast_digamma`. Add preallocated `m_w_grad(n)` member field (sized in constructor, reused across calls). Also preallocate `m_mu(n)`, `m_eta(n)` — `operator()` is non-`const`, no LICM risk. Expected: 2–3× on `beta_est` combined with TODO-22.
Implemented with preallocated `m_eta`, `m_mu`, and `m_w_grad` members and a single scalar observation loop that accumulates negative log-likelihood, beta-gradient weights, and the log-precision derivative while reusing `fast_lgamma` and `fast_digamma`. Correctness: beta-related coverage in `test-fast_glm_outputs.R`, `test-argument-permutations.R`, `test-rcpp-fitting-equivalence.R`, and `test-rcpp-fitting-real-data.R` passed after rebuild. Targeted adaptive benchmarks using the `benchmark_model_fits.R` timing setup and the fresh pre-change TODO-22 baseline: `beta_est` 1.294521ms → 1.100000ms (−15.0%); `beta_var` 0.356529ms → 0.344978ms (−3.2%).

**TODO-24: LogBinomial — apply TODO-17 backtracking fix to weighted IRLS variant** ✓ DONE
File: `EDI/src/fast_log_binomial_regression.cpp:285–421`
`fit_constrained_binomial_weighted_cpp_impl` was skipped when TODO-17 optimized the unweighted variant (4.4× speedup). The weighted path has identical structure: `ll_curr` recomputed inside the loop (line 397), each backtracking probe calls `weighted_loglik_constrained_binomial` with a full GEMV, and `eta_try` is not preallocated. Add: (1) `weighted_loglik_from_eta(eta, y, obs_weights, link_type)` helper (O(n) scalar loop, no GEMV); (2) initialize `ll_curr` before the IRLS loop; (3) precompute `delta_eta = X_free * direction` once per iteration; (4) preallocate `eta_try(n)` outside the while-loop; (5) carry forward `ll_curr = ll_new`. Expected: ~4× on weighted logbin/identbin paths.
Implemented `weighted_loglik_from_eta`, cached and carried forward `ll_curr`, precomputed `delta_beta` / `delta_eta` once per weighted IRLS iteration, reused preallocated `z_adj`, `w_eff`, and `eta_try`, and moved weighted working-vector setup into one scalar loop. Correctness: `test-bayesian-bootstrap.R` passed under `pkgload::load_all()` (58/58); `test-warm-start-weights.R`, `test-rcpp-fitting-equivalence.R`, and `test-fast_glm_outputs.R` passed after rebuild. Targeted adaptive benchmarks using the `benchmark_model_fits.R` timing setup: `logbin_weighted_est` 8.648148ms → 6.306452ms (−27.1%); `identbin_weighted_est` 2.417808ms → 1.769565ms (−26.8%).

**TODO-25: Wilcox HL — O(n²) → O(n log²n) sort + binary-search estimator** ✓ DONE
File: `EDI/src/fast_wilcox_hl.cpp`
`hl_from_groups` materializes all `n_t × n_c` pairwise differences into a `std::vector<double>` then calls `nth_element` — O(n²) space and time. `hl_signed_rank` materializes `m*(m+1)/2` Walsh averages. Profile: 81% of `wilcox_est` samples in `median_in_place`; 16.3% branch-miss rate from `nth_element` introselect on an O(n²) array. Replace with: sort `y_t` and `y_c` (O(n log n)); binary-search for the median-rank threshold. Count of pairs `y_t[i] - y_c[j] ≤ θ` is O(n log n) via `upper_bound` on sorted `y_t` for each `y_c[j]`. Similarly for Walsh average median. Eliminates the O(n²) allocation. Expected: ~13–14s wall time reduction on `wilcox_est`.
Implemented sorted implicit selection for rank-sum pairwise differences and signed-rank Walsh averages. Small inputs still use the old materialized exact path; larger inputs sort once, binary-search the median rank using O(n) two-pointer counts, then snap to the largest implicit value below the selected threshold. Correctness: direct reference checks against R materialization passed for rank-sum and signed-rank cases, including duplicate-heavy and benchmark-sized inputs; `test-design-inference.R` and `test-bayesian-bootstrap.R` passed under `pkgload::load_all()`; `test-wilcox-regr-bootstrap-fast-path.R` had only its two pre-existing skips. Targeted adaptive benchmarks using the `benchmark_model_fits.R` timing setup: `wilcox_est` 0.545093ms → 0.126162ms (−76.9%); `wilcox_var`-sized N=200 path 0.073454ms → 0.039420ms (−46.3%).

**TODO-26: JT exact — precomputed lchoose table + flat `std::vector` for stat_prob** ✓ DONE
File: `EDI/src/fast_jonckheere_terpstra.cpp`
Two issues: (1) `log_choose_int` calls `R::lchoose` on every recursive call — `Rf_chebyshev_eval` (10.8%) + `Rf_lgammacor` (5.6%) + `Rf_lbeta` (4%) + `Rf_lchoose` (3.8%) = 24% of samples from redundant evaluations on repeated `(nk, tk)` pairs. Fix: before recursion, precompute `lc[k][t] = lchoose(n_k, t)` for `t=0..n_k` as a flat array (O(n) total), pass as `const double*` into recursion, eliminating all `R::lchoose` calls inside `recurse_jt_distribution`. (2) `stat_prob` is `std::map<int,double>` — red-black tree with `malloc` per insert (9% of samples in tree+alloc). The JT statistic `stat2` ranges in `[0, 2*n_treat*n_control]` — known before recursion. Replace with `std::vector<double>(max_stat2+1, 0.0)` indexed directly.
Implemented with a flat log-choose table backed by cumulative log-factorials, a flat `std::vector<double>` distribution indexed by `stat2`, active-stat tracking to avoid a dense final scan, thread-local distribution buffer reuse to clear only active slots between calls, and incremental recursion that carries the JT statistic forward instead of rebuilding treatment-count statistics at each leaf. Correctness: direct exact enumeration reference checks passed for tied ordinal samples; `test-asymp-inference-paths.R` and `test-design-inference.R` passed after rebuild. Targeted adaptive benchmark using the `benchmark_model_fits.R` / `edi_kernel_profiler.R` `jt_var` expression at N=200: 0.024382ms → 0.007039ms median (−71.1%).

**TODO-27: Ridit — `std::map` → `std::unordered_map` + eliminate `wrap(ref_idx)` SEXP round-trip** ✓ DONE
File: `EDI/src/fast_ridit_analysis.cpp`
(1) Lines 19, 38: `std::map<int,int>` and `std::map<int,double>` are O(log K) per access. Replace with `std::unordered_map`. (2) Line 95: `fast_ridit_scores_cpp(y_sexp, wrap(ref_idx))` converts `std::vector<int>` → SEXP → re-reads it inside the function. `Rf_allocVector3` is 6.6% of samples from this unnecessary round-trip. Refactor `fast_ridit_scores_cpp` to accept `const std::vector<int>&` via a static helper, or inline the logic. (3) Line 122: `std::pow(s - mean_ridit_t, 2)` → `(s-mean)*(s-mean)`. Expected: ~3–4s wall time reduction on `ridit_var`.
Implemented with flat sorted level/count/ridit vectors rather than hash tables after a literal `std::unordered_map` rewrite benchmarked slightly slower on the N=200, K=3 ordinal workload. The final implementation removes all `std::map` use, eliminates the internal `wrap(ref_idx)` SEXP round-trip by using a zero-based C++ helper, avoids an internal `List` construction in `fast_ridit_analysis_cpp`, removes temporary treatment/control score vectors, and replaces `std::pow(diff, 2)` with `diff * diff`. Correctness: direct R reference checks passed for score assignment and analysis outputs, including categories absent from the reference group; `test-bayesian-bootstrap.R` passed under `pkgload::load_all()` (58/58). Targeted adaptive benchmark using the `benchmark_model_fits.R` / `edi_kernel_profiler.R` `ridit_var` expression at N=200: 0.008617ms → 0.006072ms median (−29.5%).

**TODO-28: `weighted_crossprod` col-major — 2×GEMM → DSYR rank-1 accumulation**
File: `EDI/src/_helper_functions.h`
The col-major branch of `weighted_crossprod` evaluates `X.T * w.asDiagonal() * X` as two GEMMs: first `X.T * diag(w)` produces a p×n intermediate matrix (heap allocation + n*p multiplications), then multiplies by X (second GEMM p×n × n×p → p×p). Replace with `selfadjointView<Upper>().rankUpdate()` (BLAS DSYR): n symmetric rank-1 updates each O(p²/2), no p×n intermediate. For n=1000, p=5: current ~50k FMAs + p×n allocation vs ~12.5k FMAs + 0 allocation. Accounts for 45% of `robust_est` (called per IRLS iteration) and all weighted OLS paths. **Caveat:** For very small p (p=3–5), Eigen may already beat BLAS for DSYR — benchmark immediately; revert if slower.

**TODO-29: `fast_ols.cpp` XtX → `selfadjointView.rankUpdate` (DSYRK)**
File: `EDI/src/fast_ols.cpp:32,108`
`XtX_free.noalias() = X.transpose() * X` uses a full GEMM (77.5% of `ols_est` samples in `lhs_process_one_packet`). Replace with:
```cpp
XtX_free.setZero();
XtX_free.selfadjointView<Eigen::Upper>().rankUpdate(X.transpose());
XtX_free.triangularView<Eigen::Lower>() = XtX_free.transpose();
```
Uses DSYRK (half the FLOPs). Expected: 30–50% speedup on `ols_est`. Same caveat as TODO-28: benchmark for small p.

**TODO-30: CoxPH `compute_robust_vcov` — eliminate 150 heap allocs per call**
File: `EDI/src/fast_coxph_regression.cpp:498–530`
`ek`, `dk_ek_over_Rk`, `cum_B` are `std::vector<std::vector<double>>(n_events, std::vector<double>(p, 0.0))` — n_events separate heap allocations each. At n_events=50, p=5: 150 separate heap allocs. `cluster_scores` is `std::map<int, Eigen::VectorXd>` (O(log n) per insert). Fix: replace nested vectors with flat `Eigen::MatrixXd(p, n_events)` indexed as `[:,k]`; replace `std::map` with `std::unordered_map`. Also preallocate `cluster_scores` at construction when n_clusters is known.

**TODO-31: Ordinal regression — compute hessian once, not 3–4 times**
File: `EDI/src/fast_ordinal_regression.cpp:183–185,311`
`fast_ordinal_regression_cpp` calls `model.hessian(params)` 3× in succession for `observed_information`, `fisher_information`, and `information` (all identical). Then `fast_ordinal_regression_with_var_cpp` (line 311) calls a 4th time. Fix:
```cpp
MatrixXd H = model.hessian(params);
return List::create(Named("observed_information") = H, Named("fisher_information") = H, Named("information") = H, ...);
```
In `fast_ordinal_regression_with_var_cpp`, extract `H` from `res["observed_information"]` instead of re-calling. Expected: 3× faster on hessian-dominated var path.

**TODO-32: GComp ordinal finite-diff loop — hoist temporaries outside loop**
File: `EDI/src/fast_ordinal_regression.cpp:449–462`
Finite-difference loop over `n_params` iterations allocates `p_plus`, `p_minus` (VectorXd copies of `params`) plus inside `compute_md_from_params` (line 421): `alpha`, `beta`, `eta1_loc`, `eta0_loc` — ≥6 heap allocs × n_params ≈ 8 = ~48 allocs per post-fit call. Fix: hoist `p_plus` and `p_minus` outside loop, reset from `params` inside; pass pre-sized scratch vectors to `compute_md_from_params` by reference. Expected: 30–50% faster on gcomp ordinal post-fit path.

**TODO-33: HurdlePoisson GLMM — `log1p` fast-path for large lambda**
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:183,268`
`__log1p_fma` dominates both `hurdle_p_est` and `hurdle_p_var` annotate output. Line 183: `const double lne = (lam < 1e-10) ? eta_ki : std::log1p(-eneg)`. For `lam > 16` (`eneg < 1e-7`), `log1p(-eneg) ≈ -eneg` with error < 1e-7. Add fast-path:
```cpp
const double lne = (lam < 1e-10) ? eta_ki : (eneg < 1e-7 ? -eneg : std::log1p(-eneg));
```
Apply identically in `hessian()` line 268. If typical lambda values are moderate-to-large (common in count outcomes), eliminates 30–70% of `log1p` calls. Verify numerical accuracy on held-out fit before committing.

**TODO-34: Continuation ratio — eliminate O(n×K) VectorXd allocs in augmented data construction**
File: `EDI/src/fast_continuation_ratio_regression.cpp:64–88`
`build_continuation_ratio_augmented_data` allocates one `VectorXd x_row(n_alpha + p)` per observation × level combination, pushes into `std::vector<VectorXd>`, then copies into `MatrixXd X_aug`. Profile: 71% of `cont_ratio_est` samples in `unlink_chunk`, 64% in `_int_free_create_chunk` — entirely malloc/free overhead. Fix: build `X_aug` and `z` directly:
```cpp
int total_rows = 0;
for (int i = 0; i < n; ++i) total_rows += std::min(y_level[i] + 1, n_alpha);
MatrixXd X_aug(total_rows, n_alpha + p); X_aug.setZero();
VectorXd z(total_rows);
int row = 0;
for (int i = 0; i < n; ++i) {
    for (int j = 0; j < std::min(y_level[i]+1, n_alpha); ++j, ++row) {
        X_aug(row, j) = 1.0;
        if (p > 0) X_aug.row(row).tail(p) = X.row(i);
        z[row] = (y_level[i] == j) ? 1.0 : 0.0;
    }
}
```
Expected: 3–5× speedup on `cont_ratio` (from ~20s down to ~5s).

**TODO-35: HurdlePoisson GLMM `operator()` — `Xg` expression view + `eta0` preallocation**
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:154,244`
Line 244 in `hessian()`: `const Eigen::MatrixXd Xg = dat.X_s.middleRows(gs, sz)` — copies the group submatrix G times per call. Change to `const auto Xg = dat.X_s.middleRows(gs, sz)` (zero-copy Eigen expression view). Line 154 in `operator()`: `const Eigen::VectorXd eta0 = dat.X_s.middleRows(gs, sz) * par.head(p)` — allocates size-sz vector G times. Rewrite using preallocated `m_eta0.head(sz).noalias() = dat.X_s.middleRows(gs, sz) * par.head(p)`.

### MEDIUM priority

**TODO-36: HurdlePoisson GLMM — fuse `log_sum_exp` + posterior weights into single pass**
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:189,200,273,280`
Both `operator()` and `hessian()` call `log_sum_exp_hp(log_terms)` then immediately loop over `exp(log_terms[k] - ll_g)` to get posterior weights — computing each `exp(log_terms[k])` twice. Use `log_sum_exp_and_weights` from `_glmm_engine.h` (already exists) to do both in one pass. Saves K exp calls per group per optimizer step. For G=200, K=7: 1,400 exp calls eliminated per iteration.

**TODO-37: Poisson IRLS — `X.T * w.asDiag() * X` → `weighted_crossprod`**
File: `EDI/src/fast_poisson_regression.cpp:227,270`
Lines 227 and 270: `XtWX_free.noalias() = X_f.transpose() * w_tmp.asDiagonal() * X_f` creates a p×n intermediate matrix per IRLS iteration. Replace with `XtWX_free = weighted_crossprod(X_f, w_tmp)` (already in `_helper_functions.h`). Applies equally to quasi_var (same code path). See TODO-28 caveat for small p.

**TODO-38: Robust regression MAD — `std::sort` → `std::nth_element`**
File: `EDI/src/fast_robust_regression.cpp:92–96`
MAD scale estimation calls `std::sort(abs_r)` on residuals but only needs the median — O(n log n) when O(n) suffices. Replace with `std::nth_element` to position n/2 (plus a second nth_element for even n to get the lower median). ~7× faster on sort step; eliminates 7% of `robust_est` wall time.

**TODO-39: Probit — reuse `mu[]` for post-convergence neg_ll computation**
File: `EDI/src/fast_probit_regression.cpp:207–213`
After IRLS convergence, `mu[i]` already holds `pnorm_fast(eta_converged[i])` from the final iteration. The current code recomputes `final_eta` (one GEMV) and then calls `log_pnorm_lower`/`log_pnorm_upper` (2 erfc calls per obs). Replace with: `neg_ll -= wi * (y[i] > 0.5 ? log(mu[i]) : log1p(-mu[i]))`, eliminating 1 GEMV + n erfc calls. Guard: only valid for the IRLS convergence path, not the lbfgs path (lines 137–143) where `mu` is not populated.

**TODO-40: Probit — polynomial fast-erfc approximation**
File: `EDI/src/fast_probit_regression.cpp` and `_helper_functions.h`
erfc calls account for ~75% of all `probit_est` samples. For typical probit use, `η ∈ [-6, 6]`. A 6-term Horner-form minimax polynomial for `erfc(x)` on `[0, 5.6]` achieves error < 1e-10 and is ~5× faster than libm `__cr_erf_fast`. Implement as `fast_erfc` in `_helper_functions.h`; fallback to `std::erfc` for `|x| > 5.6`. Apply to `pnorm_fast`, `log_pnorm_lower`, `log_pnorm_upper`. **Requires careful accuracy validation against `R::pnorm5` on a sweep of η values before committing.** Also applies to `fast_ordinal_probit_regression.cpp` (same erfc pattern).

**TODO-41: Logistic — replace Eigen sigmoid expression with `plogis_array_safe`**
File: `EDI/src/fast_logistic_regression.cpp:138`
`mu.array() = 1.0 / (1.0 + (-eta.array()).exp())` with `EIGEN_DONT_VECTORIZE` degrades to per-element PLT `expf64` dispatch (store-latency bottleneck at `35e11e`). Replace with `mu = plogis_array_safe(eta.array()).matrix()` (already in `_helper_functions.h:376`) — explicit scalar loop with sign-branched numerically stable form + `std::exp` (resolved at link time, no PLT). Apply at lines 138, 243, 258, 276, 294.

**TODO-42: ZOIB — `fast_digamma` + `std::lgamma` + preallocate member vectors**
File: `EDI/src/fast_zero_one_inflated_beta.cpp`
(1) `DigammaFunctor` calls `R::digamma(x)` → replace with `fast_digamma(x)` (already in `_helper_functions.h`). (2) `LgammaFunctor` calls `R::lgammafn(x)` → replace with `std::lgamma(x)`. (3) Line 164: `R::digamma(phi)` per optimizer step → `fast_digamma(phi)`. (4) Lines 189–190: `R::digamma(phi)` and `R::trigamma(phi)` already hoisted but use slow dispatch → `fast_digamma`. (5) Preallocate `pi0(n)`, `pi1(n)`, `pib(n)`, `grad_gamma0(p_zero_one)`, `grad_gamma1(p_zero_one)` as member fields (currently allocated every `operator()` call).

**TODO-43: ZAP — inline `log1m_eml` reusing precomputed `eml`**
File: `EDI/src/fast_zero_augmented_poisson.cpp:71–72`
In the hurdle positive branch, `eml = std::exp(-lam)` is computed at line 71, then `log1mexp(-lam)` is called at line 72 which internally re-evaluates `exp(-lam)`. Inline as:
```cpp
const double log1m_eml = (lam > 0.6931471805599453)
    ? std::log1p(-eml) : std::log(-std::expm1(-lam));
```
Eliminates 1 `exp`/`expm1` call per positive hurdle obs (~68% of obs).

**TODO-44: MN — verify R wrapper uses `mn_ci_cpp` not R-level bisection loop**
File: `EDI/src/miettinen_nurminen_speedups.cpp`
`mn_var` shows `bcEval_loop` at 12.5% + `mn_constrained_mle_pc_cpp` at 11.4%, implying R calls the constrained MLE from R-level bisection. The C++ `mn_ci_cpp` encapsulates the full bisection (lines 116–146). Check the R wrapper — if it calls `mn_constrained_mle_pc_cpp` or `mn_z_statistic_cpp` from R in a loop, replace with a single `mn_ci_cpp` call. Would collapse ~50 R→C++ dispatches per CI into one.

**TODO-45: GComp — verify warm-start logistic fits between bootstrap replicates**
File: `EDI/R/inference_incidence_KK_gcomp_abstract.R` (and related files)
`fast_logistic_regression_internal` accounts for 49–51% of wall time for all three logistic gcomp kernels. The C++ fitter already supports `warm_start_params`. If the R bootstrap loop does not pass the previous replicate's parameters as `warm_start_params`, enabling this would be the highest-payoff change for gcomp_logistic. Check that the KK gcomp bootstrap passes `warm_start_params = previous_fit$params` between replicates.

**TODO-46: GComp ordinal — `plogis_array` Eigen temporaries → scalar loop**
File: `EDI/src/gcomp_speedups.cpp:399–405`
`compute_mean` lambda calls `plogis_array(Eigen::ArrayXd::Constant(n, alpha_hat[k]) - eta_vec.array())` for each of K-1 thresholds. Each call allocates 2 n-length Eigen arrays. For K=4, 2 arms: 2×3×2 = 12 n-length array allocs per post-fit call. Replace with a scalar loop using `plogis_stable_cpp` (already used in `ordinal_gcomp_post_fit_cpp` at lines 409–418).

**TODO-47: Logrank — fuse martingale mean/variance passes into main sweep**
File: `EDI/src/fast_logrank.cpp:85–103`
After the main event-time sweep, two separate O(n) passes compute per-group martingale mean and variance. These can be accumulated as running sums inside the main sweep (line 78 already assigns `martingale[i]` per subject). Replace post-sweep loops with: accumulate `sum_martingale` and `sum_sq_martingale` per treatment/control group inside the main loop; compute `mean = sum/n`, `var = (sum_sq - n*mean²)/(n-1)` after sweep. Eliminates 2 O(n) passes. Also: `logrank_est` has the highest branch-miss rate of any kernel (3.33%) — investigate whether the floating-point equality check `recs[end].time == curr_time` (line 58) can be replaced with integer time-group indexing.

**TODO-48: Weibull regression — preallocate member buffers in `WeibullAFTLikelihood`**
File: `EDI/src/fast_weibull_regression.cpp:25–45`
`WeibullAFTLikelihood::operator()` allocates 4 Eigen vectors (`eta`, `w`, `exp_w`, `d_ll_d_eta`) on every optimizer call. In `hessian()` (line 47): `beta_weights` and `cross_weights` are also temporary n-element vectors. Add preallocated `m_eta(n)`, `m_w(n)`, `m_exp_w(n)`, `m_d_eta(n)`, `m_beta_weights(n)`, `m_cross_weights(n)` as class members initialized in the constructor. Non-`const` methods so no LICM/alias risk.

**TODO-49: ZINB — `R::lgammafn` → `std::lgamma`**
File: `EDI/src/fast_zinb.cpp:88,94`
Line 88: `R::lgammafn(theta)` called once per step; `std::lgamma` avoids R dispatch. Line 94 in the distinct-y loop: `m_lgamma_yptheta[k] = R::lgammafn(ypt)` — nd times per step. `std::lgamma` is ~30–50% faster per call by avoiding R's error-check wrapper. Same fix already applied in `fast_beta_regression.cpp`.

**TODO-50: Adj_cat — batch exp per obs in inner loop**
File: `EDI/src/fast_adjacent_category_logit.cpp`
Profile: `adj_cat_est` is slowest in the ordinal group (27.0s, IPC=1.72, 21% `expf64`, 14% `log_fma`, 1.16% branch-miss). The adjacent-categories link evaluates `exp(η_k)` per category per observation. Restructure to compute exp(η) and softmax probabilities in a single vectorized pass per observation — one exp sweep per obs rather than one exp per category per obs. Also investigate whether the threshold monotonicity check causes the 316M branch misses; a branchless penalty formulation may help.

**TODO-51: Ord_cauchit — mutable scratch buffers to eliminate malloc overhead**
File: `EDI/src/fast_ordinal_cauchit_regression.cpp`
Profile: `ord_cauchit_est` has 43% of samples in malloc TLS arena access (23% movq %fs: TLS read + 20% TLS write). The cauchit link likely allocates intermediate arrays per optimizer call that the logistic link avoids through Eigen lazy evaluation. Add mutable scratch `VectorXd` buffers sized at construction to eliminate per-call allocs. Caution: check whether this is a `const` method inside a GLMM objective — if so, see the mutable-field antipattern warning.

**TODO-52: Ord_probit — fast_erf approximation for probit link**
File: `EDI/src/fast_ordinal_probit_regression.cpp`
Profile: erfc/erf calls spread across samples with IPC=2.24 (compute-bound). Use the same `fast_erfc` from TODO-40 once implemented. Also applies to `fast_probit_regression.cpp`.

### LOW priority

**TODO-53: CoxPH hessian — direct symmetric writes vs post-loop triangular copy**
File: `EDI/src/fast_coxph_regression.cpp:204–211`
Inner event-time loop writes only `hess.triangularView<Lower>()`, then copies to upper after the event loop (O(p²) copy per optimizer step). Write both `hess(q1,q2)` and `hess(q2,q1)` directly in the inner product loop (lines 172–177), eliminating the post-loop copy. For p≤10 the direct write + copy elimination is a net win.

**TODO-54: Robust regression — cache XtX from IRLS to avoid var-path recompute**
File: `EDI/src/fast_robust_regression.cpp:253`
After IRLS convergence, line 253 recomputes `XtX = X_free.T * X_free` for the sandwich estimator. Either store it in `RobustModelResult` at convergence, or derive from the QR cold-start (`R.T * R` where R is the QR factor). Saves one O(n·p²) GEMM per `robust_var` call.

**TODO-55: Hat-matrix diagonal — vectorize via `cwiseProduct + rowwise().sum()`**
File: `EDI/src/robust_post_fit_speedups.cpp:119–123`
`ols_hc2_setup_cpp` computes `hat[i] = XB.row(i).dot(X_fit.row(i))` in a scalar loop — the diagonal of `X_fit * bread * X_fit.T`. Replace with `hat = (X_fit * bread).cwiseProduct(X_fit).rowwise().sum()` — same result, vectorizable as element-wise product + row reduction.

**TODO-56: BetaRegression preallocate `m_mu`, `m_eta` member fields**
File: `EDI/src/fast_beta_regression.cpp`
`operator()` calls `resize()` on `eta` and `mu` vectors on every call. Preallocate as member fields sized at construction (safe — `operator()` is non-`const`). After implementing TODO-23's single scalar loop, most temporaries disappear; these are only the eta/mu GEMV intermediates.

**TODO-57: LogBin — return `fisher_information` from fit impl to avoid var-path recompute**
File: `EDI/src/fast_log_binomial_regression.cpp:453–510`
`fit_constrained_binomial_with_var_cpp_impl` rebuilds `X_free` and calls `weighted_crossprod` again (line 484) after `fit_constrained_binomial_cpp_impl` already computed `XtWX` on its final accepted iteration. Return `fisher_information` from the fit impl; reuse in the var impl. Saves one O(n·p²) GEMM per var call.

**TODO-58: Poisson IRLS — cache `delta_eta` for step-halving probes**
File: `EDI/src/fast_poisson_regression.cpp:236–248`
Each `compute_neg_loglik(beta_try)` in the halving loop (line 240) runs a full GEMV. Precompute `delta_eta = X_f * step` before the halving loop; each probe becomes O(n) vector-add + scalar loglik loop. Same pattern as TODO-17 for log-binomial. Low priority: Poisson IRLS typically converges in 1 step for well-conditioned data.

**TODO-59: Stereotype logit — preallocate per-obs hessian allocs**
File: `EDI/src/fast_stereotype_logit.cpp`
No perf data (no annotate file generated). Static analysis of `loglik_hessian()` (lines 207–307) shows per-observation allocations: `std::vector<VectorXd> logit_grad(K)`, `std::vector<MatrixXd> logit_hess(K)`, `mean_grad(d)`, `mean_hess(d,d)`, `mean_outer(d,d)`, `MatrixXd delta` per obs. For n=200, K=5, d=6: ~2,800 allocs per hessian call. Preallocate as class members. Implement only if this kernel appears in a benchmark sweep as a bottleneck.

**TODO-60: ZAP `hessian()` — reuse preallocated `m_eta_cond`/`m_eta_zi`**
File: `EDI/src/fast_zero_augmented_poisson.cpp:104–105,194–195`
`hessian()` allocates `eta_cond` and `eta_zi` locally (lines 104–105, 194–195) despite `m_eta_cond` and `m_eta_zi` existing as preallocated member fields. Reuse them with `.noalias() =`. Affects `zip_var`/`hurdle_p_var`, not the dominant `est` path.

**TODO-61: Cont_ratio `ContinuationRatioObjective::operator()` — preallocate scratch buffers**
File: `EDI/src/fast_continuation_ratio_regression.cpp:26–34`
After TODO-34 eliminates the augmented-data construction allocs, the remaining `operator()` scratch vectors (`eta`, `mu`, `log_mu`, `log_one_minus_mu`) are still allocated per L-BFGS call. Preallocate as class members. Note: `operator()` is the method signature here — check if it's `const`; if so use explicit workspace passing not mutable fields (see mutable-field antipattern).

### NegBin + HurdleNegBin findings (HIGH priority)

**TODO-62: NegBin `NBLogLik::operator()` — per-row rank-1 gradient → single GEMV**
File: `EDI/src/fast_negbin_regression.cpp:112`
`score_beta += coef * m_X.row(i).transpose()` — per-row rank-1 accumulation, same anti-pattern fixed by TODO-15/16 for ZIP/ZINB but not applied here. Fix: accumulate scalar coefficients into preallocated `m_coef_vec[i]` inside the obs loop, then `score_beta.noalias() += m_X.transpose() * m_coef_vec` (single BLAS GEMV) after. Expected: same speedup class as ZIP/ZINB GEMV fixes.

**TODO-63: NegBin `NBLogLik::hessian()` — use preallocated distinct-y digamma/trigamma tables**
File: `EDI/src/fast_negbin_regression.cpp:148–150`
`hessian()` calls `R::digamma(yi + theta)` and `R::trigamma(yi + theta)` raw per-obs, ignoring `m_digamma_yptheta` and `m_trigamma_yptheta` tables already populated in `operator()`. Compare: `TruncatedNegBinCount::hessian()` correctly uses the tables. Copy the table-lookup pattern to `NBLogLik::hessian()`. Eliminates all per-obs R function dispatch in the hessian path.

**TODO-64: NegBin `expected_hessian()` — trigamma recurrence to eliminate O(iter) R::trigamma calls**
File: `EDI/src/fast_negbin_regression.cpp:204–211`
`expected_trigamma_y_plus_theta` calls `R::trigamma(k + theta)` in a series summation loop — ~47 R::trigamma calls per obs per expected-hessian call = ~47,000 R::trigamma per call. Fix: call `R::trigamma(theta)` once, then use the recurrence `ψ₁(x+1) = ψ₁(x) − 1/x²` for subsequent terms (1 R call + O(iter) divisions). The recurrence is exact for the trigamma function — no approximation.

**TODO-65: NegBin + HurdleNegBin — `R::lgammafn` → `std::lgamma`**
Files: `EDI/src/fast_negbin_regression.cpp:83,90`, `EDI/src/fast_hurdle_negbin.cpp:346,356`
`R::lgammafn` goes through R's error-handling dispatch wrapper; `std::lgamma` is direct libm. Profile: `logf32x` (lgamma) is **63% of negbin_est** samples. This is the single largest hotspot in NegBin. Fix: replace all `R::lgammafn(x)` with `std::lgamma(x)` in both files. Combined with TODO-22 (`fast_lgamma`), reduces to: `fast_lgamma` where available. Same fix applied to ZINB in TODO-49 — extend scope to cover NegBin and HurdleNegBin here.

**TODO-66: HurdleNegBin — `log1p(-eneg)` fast-path for large `lam`**
File: `EDI/src/fast_hurdle_negbin.cpp` (hurdle positive-count log-normalizer)
Mirror of TODO-33 (HP GLMM) and TODO-43 (ZAP): for `lam > 16`, `eneg = exp(-lam) < 1e-7`, so `log1p(-eneg) ≈ -eneg` with error < 1e-14. Add fast-path branch before the `std::log1p` call. Eliminates the libm transcendental for large-count observations (common in NegBin hurdle outcomes).

**TODO-67: ZINB (and cross-kernel) — vectorize exp/log by enabling Eigen SIMD + batching transcendentals into array ops** ☐
Files: `EDI/src/Makevars:3` (+ `EDI/src/Makevars.win`), `EDI/src/fast_zinb.cpp:79-148`

Motivation: After TODO-19 (member preallocation) landed it was verified bit-exact but only ~1% faster in an apples-to-apples build. The `memset 86%` in the zinb_est profile above does **not** reproduce on this machine — glibc's dynamic mmap threshold reuses the large `VectorXd` allocations from the heap after the first few frees, so allocation was never the real bottleneck. The persistent cost is the per-observation transcendentals (zinb_est profile: `log1p 50%`, `log 47%`; every kernel with `exp`/`log`/`log1p` in its obs loop is affected). Those run scalar today for two independent reasons: (a) `EIGEN_DONT_VECTORIZE` disables Eigen packet math; (b) the loop calls `std::exp`/`std::log`/`std::log1p`, which GCC does not auto-vectorize without fast-math + libmvec.

This is two coupled changes; **(2) is a no-op without (1)** (Eigen `.array().exp()/.log()` stay scalar until the flags are removed).

Measured (double, n=20000, ns/elem, this machine):

| op | current (`DONT_ALIGN`+`DONT_VECTORIZE`) | Eigen SIMD (both flags off) | `-ffast-math`+libmvec |
|---|---:|---:|---:|
| exp | 5.0 | **1.9 (2.6x)** | 2.9 |
| log | 4.8 | **2.9 (1.6x)** | 3.3 |
| log1p | 11.7 | 11.7 (no gain) | 2.5-5.7 (2.6-4.7x) |

Verified facts:
- Both `-DEIGEN_DONT_ALIGN` and `-DEIGEN_DONT_VECTORIZE` are committed on `Makevars:3` and **each independently** forces scalar transcendentals; forcing `EIGEN_UNALIGNED_VECTORIZE=1` does not override. To vectorize, **both** must be removed (for the affected TU).
- Dropping both is **crash-safe**: audit of `EDI/src` found no `Map<...,Aligned>`, no fixed-size vectorizable Eigen types (`Vector4d`...), no `EIGEN_MAKE_ALIGNED_OPERATOR_NEW`. A deliberately misaligned default `Eigen::Map` runs vectorized (1.98 ns) without segfault — default Maps emit unaligned packet loads.
- Eigen does **not** vectorize `log1p`/`log1pexp` even with SIMD on (stays ~11.7 ns); it needs libmvec (`-ffast-math`) or a hand-written SIMD softplus. In ZINB this is only the `y>0` softplus branch (secondary in zero-inflated data).
- MKL batch `vdExp`/`vdLn` was already tried (TODO-18) and lost to scalar at n=1000 due to extra O(n) memory passes; Eigen fused array ops avoid that.

(1) Enable Eigen SIMD — two delivery options:
- **A. Global** (remove both flags from `Makevars`/`Makevars.win`): all 104 kernels' GEMV/GEMM/array ops vectorize (likely a bonus win), but FP reduction order changes package-wide -> must re-run the full 62-file / ~9.9k-LOC testthat suite and re-baseline the exact-digit fingerprints in this doc; may nudge some L-BFGS convergence paths. No crashes.
- **B. Scoped** (`#undef EIGEN_DONT_ALIGN` / `#undef EIGEN_DONT_VECTORIZE` at the top of `fast_zinb.cpp`, before includes): verified to re-enable vectorization (2.08 ns) even with the disabling `-D` flags present. Contains blast radius to ZINB, but carries a formal ODR caveat (header-only Eigen inline/templates compiled with differing `MAX_ALIGN_BYTES` across TUs) — low risk here since ZINB's Eigen types are local (crosses the TU boundary only via SEXP).

(2) Rewrite `ZeroInflatedNegBin::operator()` to batch transcendentals into vectorized array precompute passes, keeping the branchy accumulation scalar. Every per-obs transcendental has a vectorizable input:
- `m_mu = m_eta_c.array().exp()` (exp)
- `m_p  = plogis_array_safe(m_eta_z)` (exp; helper exists in `_helper_functions.h:376`)
- `m_logden = (theta + m_mu).log()` (log)
- `phi = (theta*(log_theta - m_logden)).exp()` for `y==0` (exp)
- `lse = log1pexp_array_safe(m_eta_z)` for `y>0` (stays scalar — Eigen log1p)
Reuse/extend the TODO-19 member buffers (add `m_mu`, `m_logden`); precompute writes into members, so no new per-call allocation. ~40-60 lines in one function.

Sequencing: (1) -> baseline suite/fingerprints -> (2) -> parity vs R canonical + a scalar-loop reference (NOT bit-exact vs current, since (1) reorders FP) -> interleaved A/B benchmark at n in {2000, 20000}.

Projected payoff: the dominant per-obs transcendentals (2 exp + 1 log per obs) vectorize; `log1pexp` left scalar. Amdahl on `operator()` => ~1.3-1.6x on zinb_est/zinb_var (to confirm). Option A additionally speeds every other kernel's GEMV.

Open decisions: (a) global (A) vs scoped (B) for the flag change; (b) acceptable to re-baseline exact-digit fingerprints in this doc; (c) add libmvec/custom SIMD for `log1pexp` now or defer. Recommend prototyping via B first (contained), then graduating to A as a separately-validated change if the cross-kernel GEMV win is wanted.

**TODO-68: `-ffast-math` on top of TODO-67 — NOT RECOMMENDED; scalar `fast_log1pexp` alternative also FALSIFIED** ✗ DROPPED
Files: (would-be) `EDI/src/Makevars` / per-object flags, `EDI/src/_helper_functions.h`, `EDI/src/fast_zinb.cpp`

**RESULT (2026-07-02) — DROPPED.** Both `-ffast-math` (analysis below, still valid) and the recommended `fast_log1pexp` alternative were implemented and measured; **neither is worth doing.** A scalar `fast_log1pexp` (Kahan `log1p` via `std::log`) validated to **3.9e-16** rel err (edge cases bit-identical to `std::log1p`) gave **no speedup**: on the valid softplus domain `std::log1p` is ~11.2 ns — comparable to `std::log` (~8.8) and `fast_log1p` (~9.4) — not the slow op the profile implied. The earlier "log1p 50% / 14.9 ns" figure was inflated by feeding `log1p` invalid `< -1` inputs (glibc NaN/errno path) plus Eigen-array wrapper overhead; on valid inputs it is not a bottleneck. `std::exp` (~5 ns) dominates softplus and the helper does not touch it; measured full softplus **fast=12.5 ns vs std=11.2 ns (0.89x — a slight regression)**. Crucially the `-ffast-math` "4.9x log1p" is a **vectorization** (libmvec) effect that a scalar helper cannot reproduce — the only real softplus win is SIMD, i.e. **TODO-67 route 1** (an Eigen-array softplus with vectorizable `log1p` arithmetic, once Eigen SIMD is enabled). No source changes were kept.

Question: does adding `-ffast-math` on top of TODO-67 (Eigen SIMD + array-op batching) buy anything? **Conclusion: no meaningful net upside; disproportionate correctness blast radius. Get the one real win (log1p/softplus) via a targeted polynomial helper instead — the codebase's established pattern (TODO-2 `fast_digamma`, TODO-22 `fast_lgamma`, TODO-40 `fast_erfc`).**

Incremental speed over Eigen SIMD (routes 1+2), measured n=20000 ns/elem:

| op | Eigen SIMD (1+2) | + `-ffast-math` | incremental |
|---|---:|---:|---|
| exp (Eigen array) | 1.92 | 2.03 | none (marginally worse) |
| log (Eigen array) | 2.93 | 3.00 | none |
| log1p (Eigen array) | 11.7 | 2.39 | **4.9x — the only gain** |

So on top of 1+2, `-ffast-math` only helps `log1p` (Eigen already beats libmvec on exp/log). And that gain is **inseparable from `-ffinite-math-only`** — safe subsets do not autovectorize: `-fno-math-errno` -> scalar log1p 16.3; `+ -funsafe-math-optimizations` -> 13.8; only full `-ffast-math` -> 5.8. GCC emits libmvec vector calls only under finite-math assumptions.

Proven hazards (measured on this machine):
- **Guard elision**: an inlinable `if(!std::isfinite(inv)) return fallback` returns the fallback under a plain build but returns `Inf` (guard removed) under both `-ffast-math` and isolated `-ffinite-math-only`. The ZINB call path has ~15 such guards inlined into `fast_zinb.o` via `_helper_functions.h` — `covariance_from_information` (`!information.allFinite()`), `symmetric_pseudo_inverse` (`!Msym.allFinite()`), `safe_ols_solve`/`try_safe_ols_solve` (`!X.allFinite()`, `!beta_out.allFinite()`), `vector_is_usable_start`, and `quiet_NaN()` sentinels. Eliding them turns controlled NaN/pseudo-inverse fallbacks (singular Fisher info, divergent fits) into silent garbage variance/CIs. Elision is optimizer-visibility-dependent -> a latent, non-deterministic regression.
- **Process-global FTZ/DAZ**: `-ffast-math` links `crtfastmath.o` (confirmed on the link line via `-###`), whose `.init_array` ctor runs at `dyn.load` and flips flush-to-zero for the **entire R session** (every package + BLAS/LAPACK). A `1e-318` denormal: plain -> preserved; `-ffast-math` compile+link -> `0`; `-ffast-math` **compile-only + plain link -> preserved**. Mitigation: keep `-ffast-math` off the link line (R's SHLIB link is separate from `PKG_CXXFLAGS`; use a target-specific compile var, never `PKG_LIBS`).
- **Reassociation** (`-fassociative-math`) changes reductions beyond ULP and makes output compiler/version-dependent, compounding route (1)'s reproducibility cost.

Scoping analysis: compile-only per-TU removes the FTZ hazard (verified) but NOT guard elision (a compile semantic applied to the inlined helpers in that TU), and carries the same ODR caveat as TODO-67 option B.

~~Recommended alternative: `fast_log1pexp` polynomial helper in `_helper_functions.h`~~ — **FALSIFIED (see RESULT above):** a scalar/Kahan `fast_log1pexp` does not autovectorize on its own and is not faster than `std::log1p` on valid inputs; the softplus speedup exists only as SIMD array ops under TODO-67 route 1.

Decision: **DROPPED.** Do not add `-ffast-math` (hazards above), and do not add a scalar `fast_log1pexp` (no gain / slight regression, measured). Defer any softplus speedup to TODO-67 route 1 as an Eigen-array softplus (packet `exp` + vectorizable `log1p` arithmetic) once Eigen SIMD is enabled.

---

## Bottom Line

The highest-payoff optimization work is structural:

- fewer temporary vectors and full-array passes
- less per-node/per-row allocation
- special-function hoisting and tabulation (not approximation, at least until digamma approximation is validated)
- specialized early-return paths in native fitters for callers that do not need Hessian outputs

The remaining headroom is concentrated in:
1. ~~ordinal CLMM exp reduction (TODO-1)~~ ✓ DONE — `fill_alpha` precomputes thresholds once per optimizer step; inner loop now does O(1) array lookups. Re-profile to confirm speedup.
2. ~~ZINB digamma approximation (TODO-2)~~ ✓ DONE — `fast_digamma` replaces `R::digamma` with a 5-term asymptotic expansion + recurrence shift; ≤4e-12 relative error. Re-profile to confirm expected 3–5× digamma speedup.
3. ~~d-optimal search algorithmic pruning (TODO-3)~~ ✓ DONE — sorted-candidate pruning with early termination in both i and j loops; same prune structure for A-optimal via combined `obj_curr`-weighted score.

Handwritten assembly is not the right next step for any of these.

---

## Phase 6 — 57-Kernel Parallel Perf Sweep (2026-07-02)

Parallel profiling of all 57 annotate files (+ 4 stat files each) using 7 concurrent fork agents, each covering a kernel family. Ordinal regression findings still pending (agent in progress). New TODO items numbered 19+.

### Kernel stats summary

| Kernel | Wall (s) | IPC | Branch miss | Top symbols |
|---|---:|---:|---:|---|
| zip_est | 42.6 | 2.04 | 0.77% | log1p 57%, exp 15%, GEMV 39% |
| zip_var | 14.6 | 2.01 | 0.60% | GEMV 19%, exp 13%, log 14% |
| zinb_est | **219.7** | 2.08 | 1.03% | **memset 86%**, log1p 50%, log 47% |
| zinb_var | 25.4 | 2.02 | 0.71% | log 13%, GEMV 13%, log1p 11% |
| negbin_est | 18.4 | 2.31 | n/a | log 64%, exp 19%, GEMV 11% |
| negbin_var | 16.5 | 1.93 | n/a | log 36%, lbfgs-wrapper 34%, exp 19% |
| hurdle_nb_est | 17.7 | 1.89 | n/a | **malloc 33%**, lbfgs-wrapper 21%, exp 17% |
| hurdle_nb_var | 24.0 | 1.81 | n/a | lbfgs-wrapper 21%, exp 19%, log1p 10% |
| poisson_est | 19.4 | 1.78 | 0.35% | — |
| poisson_var | 18.1 | 1.84 | 0.35% | — |
| quasi_var | 15.6 | 1.84 | 0.33% | — |
| poisson_robust_var | 16.4 | 1.58 | 0.56% | log1p heavy, R GC |
| hurdle_p_est | 16.3 | 1.88 | 1.07% | log1p 50%+ |
| hurdle_p_var | **66.6** | 2.14 | 0.41% | hessian ~50s of 66.6s; 6200 allocs/call |

### Critical findings (action items TODO-19+)

---

**TODO-19: ZINB — preallocate eta_c, eta_z, w_c, w_z as member vectors** ☐
File: `EDI/src/fast_zinb.cpp:83-84,100-101`

Root cause: `ZeroInflatedNegBin::operator()` allocates four `VectorXd(m_n)` on every call:
```cpp
const Eigen::VectorXd eta_c = m_Xc * par.head(m_pc);  // line 83 — alloc + GEMV
const Eigen::VectorXd eta_z = m_Xz * par.segment(...); // line 84 — alloc + GEMV
Eigen::VectorXd w_c = Eigen::VectorXd::Zero(m_n);      // line 100 — alloc + memset
Eigen::VectorXd w_z = Eigen::VectorXd::Zero(m_n);      // line 101 — alloc + memset
```
The profiler shows 85.88% of zinb_est cycles in `rep stosb` (memset), causing 219.7s runtime (8.6× longer than zinb_var). With thousands of L-BFGS iterations, each allocating 4×8KB and zeroing 2×8KB, this dominates completely.

Fix: add member fields and preallocate in constructor — exactly what `ZeroAugmentedPoisson` already does (`fast_zero_augmented_poisson.cpp:33,42`):
```cpp
// In class ZeroInflatedNegBin — add members:
Eigen::VectorXd m_eta_c, m_eta_z, m_w_c, m_w_z;
// Constructor init-list: m_eta_c(n), m_eta_z(n), m_w_c(n), m_w_z(n)

// In operator(), replace lines 83-84,100-101:
m_eta_c.noalias() = m_Xc * par.head(m_pc);
m_eta_z.noalias() = m_Xz * par.segment(m_pc, m_pz);
m_w_c.setZero();
m_w_z.setZero();
```
Also replace `grad.resize(...)` (line 139) with a fixed-size write (grad is already sized by the caller).

Expected: eliminates 86% of zinb_est wall time — likely 10-50× speedup on the est path.

---

**TODO-20: HurdlePoisson GLMM hessian — preallocate G×K working buffers** ✓ DONE
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:234,275-304`

Root cause: `hessian()` allocates per-group working matrices inside the G-group loop: `E_Hik(total,total)`, `E_GiGiT(total,total)`, `G_avg(total)`, and per quadrature-node inner structures (`res_k(sz)`, `d2e_k(sz)`, `G_ik(total)`, `H_ik(total,total)`). For G=200 groups, K=7 nodes: ~6,200 heap allocations per hessian call. Additionally, `const Eigen::MatrixXd Xg = dat.X_s.middleRows(gs, sz)` at line 244 copies the submatrix G times.

The profiler shows hurdle_p_var = 66.6s vs hurdle_p_est = 16.3s — the hessian path adds ~50s. Without the hessian this kernel is fine.

Fix:
1. Change `const Eigen::MatrixXd Xg = dat.X_s.middleRows(gs, sz)` → `const auto Xg = dat.X_s.middleRows(gs, sz)` (expression template, zero copy).
2. Add member fields sized at construction (total = p+1, max_grp_sz = max sz across groups):
```cpp
std::vector<double> m_exp_bvals;  // K
Eigen::VectorXd m_eta0;           // max_grp_sz
Eigen::VectorXd m_res_k, m_d2e_k; // max_grp_sz
Eigen::VectorXd m_G_ik;           // total
Eigen::MatrixXd m_H_ik, m_E_Hik, m_E_GiGiT; // total×total
Eigen::VectorXd m_G_avg;          // total
```
Use `.setZero()` at the start of each group iteration instead of constructing new matrices. Same fix for `operator()` line 141: `std::vector<double> exp_bvals(K)` → `m_exp_bvals`.

Implemented with preallocated buffers plus in-place crossproduct helpers to avoid hidden Eigen temporaries. Correctness: `test-glmm-cpp-equivalence.R` passed (11/11). Targeted full-inference GLMM benchmark using the `benchmark_model_fits.R` adaptive timing setup with the profiler's GLMM expression: 1.936842ms → 1.905263ms median (N_WALD=200; 1.6% improvement at this scale).

---

**TODO-21: NBLogLik — GEMV refactor for gradient in operator()** ☐
File: `EDI/src/fast_negbin_regression.cpp:95,112,118`

Root cause: line 112 uses per-row rank-1 update:
```cpp
score_beta.noalias() += coef * m_X.row(i).transpose();  // stride-n access, not BLAS
```
This is the same anti-pattern fixed by TODO-15 (ZIP) and TODO-16 (ZINB) — row access on column-major X produces scattered writes that prevent BLAS vectorization.

Fix (same pattern as fast_zero_augmented_poisson.cpp:93):
```cpp
// Add member:  Eigen::VectorXd m_coef_vec;  // size m_n, preallocated
// In operator(), remove score_beta local; fill m_coef_vec[i] in the obs loop:
m_coef_vec[i] = yi - mu_i * (yi + theta) / denom;
// After loop, replace line 118:
grad.head(m_p).noalias() = -(m_X.transpose() * m_coef_vec);
```
Also remove `Eigen::VectorXd score_beta = Eigen::VectorXd::Zero(m_p)` (line 95) — no longer needed.

---

**TODO-22: NBLogLik hessian — fill distinct-y tables; use slot lookups** ☐
File: `EDI/src/fast_negbin_regression.cpp:148-150`

Root cause: `hessian()` calls `R::digamma(yi + theta)` and `R::trigamma(yi + theta)` per-obs (lines 148-150) without using the preallocated distinct-y tables. The tables `m_digamma_yptheta` and `m_trigamma_yptheta` exist but are only filled in `operator()`, not in `hessian()`. The hessian also uses `R::digamma(theta)` un-hoisted until line 148 (actually just-in-time inside the obs loop since theta is the same for all obs).

Fix: at the top of `hessian()`, before the obs loop:
```cpp
const double log_theta  = std::log(theta);
const double digamma_th = fast_digamma(theta);
const double trigamma_th = R::trigamma(theta);
for (int k = 0; k < nd; ++k) {
    const double ypt = m_distinct_y[k] + theta;
    m_digamma_yptheta[k]  = fast_digamma(ypt);
    m_trigamma_yptheta[k] = R::trigamma(ypt);
}
// In obs loop, replace lines 148-150:
const int slot = m_y_slot[i];
double A = m_digamma_yptheta[slot] - digamma_th + log_theta - std::log(denom) + ...
double dA_dtheta = m_trigamma_yptheta[slot] - trigamma_th + ...
```
Reduces from n R::digamma + n R::trigamma to nd R::trigamma + nd fast_digamma per hessian call (nd ≪ n for typical count data).

Note: `TruncatedNegBinCount::hessian()` in `fast_hurdle_negbin.cpp` already does this correctly (lines 414-420); copy that pattern to `NBLogLik::hessian()`.

---

**TODO-23: NBLogLik::expected_trigamma_y_plus_theta — trigamma recurrence** ☐
File: `EDI/src/fast_negbin_regression.cpp:204-211`

Root cause: inner series loop calls `R::trigamma(k+1 + theta)` at line 207 on every iteration. The series converges after ~`mean + 10*sd` terms (e.g., for mu=5, theta=2: ~47 iterations). With n=1000 obs, this is ~47,000 R::trigamma calls per `expected_hessian()` invocation.

Fix: use the trigamma recurrence `ψ₁(x+1) = ψ₁(x) − 1/x²` — same recurrence that exists for digamma:
```cpp
double trig_cur = R::trigamma(theta);  // one R::trigamma before loop
double x = theta;                       // tracks theta + k
sum = pk * trig_cur;
for (int k = 0; k < max_iter; ++k) {
    pk *= (static_cast<double>(k) + theta) / static_cast<double>(k + 1) * ratio_base;
    trig_cur -= 1.0 / (x * x);         // advance: trigamma(x+1) = trigamma(x) - 1/x²
    x += 1.0;
    sum += pk * trig_cur;
    cdf += pk;
    if (k + 1 > min_iter && pk < 1e-14 && 1.0 - cdf < 1e-12) break;
}
```
Replaces O(min_iter) R::trigamma calls with 1 R::trigamma + O(min_iter) divisions per obs per hessian call.

---

**TODO-24: R::lgammafn → std::lgamma in NegBin, ZINB, HurdleNegBin** ☐
Files: `fast_negbin_regression.cpp:83,90`, `fast_zinb.cpp:68,88,94`, `fast_hurdle_negbin.cpp:330,346,356`

Root cause: `R::lgammafn(x)` routes through R's error-handling dispatch (sets errno, checks R's interrupt flag, handles edge cases via R's machinery). `std::lgamma(x)` is a direct libm call, ~2-3× faster for normal positive inputs. The profiler shows `logf32x` (glibc log, called from lgamma) at 63.52% of negbin_est and 47.35% of zinb_est.

Specific replacements:
- `fast_negbin_regression.cpp:83`: `R::lgammafn(theta)` → `std::lgamma(theta)`
- `fast_negbin_regression.cpp:90`: `R::lgammafn(ypt)` → `std::lgamma(ypt)` (in table fill loop)
- `fast_zinb.cpp:68`: `R::lgammafn(...)` in constructor (called nd times at construction, not per optimizer step — lower priority)
- `fast_zinb.cpp:88`: `R::lgammafn(theta)` → `std::lgamma(theta)` (called every optimizer step)
- `fast_zinb.cpp:94`: `R::lgammafn(ypt)` → `std::lgamma(ypt)` (table fill loop, called every step)
- `fast_hurdle_negbin.cpp:330`: constructor (once) — lower priority
- `fast_hurdle_negbin.cpp:346`: `R::lgammafn(theta)` → `std::lgamma(theta)` (per step)
- `fast_hurdle_negbin.cpp:356`: `R::lgammafn(ypt)` → `std::lgamma(ypt)` (table fill, per step)

Note: TODO-14 already applied this to `fast_beta_regression.cpp`; same pattern here.

---

**TODO-25: HurdlePoisson/ZIP — log1p fast-path for large lambda** ☐
Files: `EDI/src/fast_hurdle_poisson_glmm.cpp:183,268`, `EDI/src/fast_zero_augmented_poisson.cpp:12-13`

Root cause: `std::log1p(-eneg)` where `eneg = exp(-lam)` is the top hotspot in hurdle_p_est (log1p at 50%+). For large lam, `eneg → 0` and `log1p(-eneg) ≈ -eneg` with error < `eneg²/2`. For lam > 30: `eneg < 9e-14`, so the approximation has error < 4e-27 (far below double precision). A threshold of `lam > 16` gives error < 1e-14.

Fix for fast_hurdle_poisson_glmm.cpp:
```cpp
const double lne = (lam < 1e-10) ? eta_ki :
                   (eneg < 1e-7)  ? -eneg  :   // lam > 16: log1p(-eneg) ≈ -eneg
                                    std::log1p(-eneg);
```
Apply same threshold at line 268 in hessian().

Fix for fast_zero_augmented_poisson.cpp `log1mexp(x)` helper (line 12-13): for `x < -16`, i.e., `exp(x) < 1e-7`: `log1mexp(x) = log(1 - exp(x)) ≈ -exp(x)` (avoids log1p call entirely):
```cpp
double log1mexp(double x) {
    if (x >= 0) return -std::numeric_limits<double>::infinity();
    if (x < -16.0) return -std::exp(x);   // exp(x) < 1e-7; log1p(-exp(x)) ≈ -exp(x)
    if (x > -0.693) return std::log(-std::expm1(x));
    return std::log1p(-std::exp(x));
}
```

---

**TODO-26: Logistic/Probit/Poisson IRLS — XtWX via weighted_crossprod** ☐
Files: `EDI/src/fast_logistic_regression.cpp`, `EDI/src/fast_probit_regression.cpp`, `EDI/src/fast_poisson_regression.cpp:227,270`

Root cause: IRLS computes XtWX as `X.T * diag(w) * X` or similar triple-product, creating an n×p intermediate. `weighted_crossprod(X, w)` (already in `_helper_functions.h`) uses the upper-triangular DSYR/DSYRK symmetric update — halves FLOPs and avoids the intermediate allocation. Same fix was applied to log-binomial (TODO-17).

For Poisson, profiler confirms this at lines 227 and 270 (`XtWX_free.noalias() = X_f.transpose() * w_tmp.asDiagonal() * X_f`). Quasi-Poisson uses the same internal path — one fix covers both.

---

**TODO-27: OLS/Robust — symmetric XtX via weighted_crossprod / DSYRK** ☐
Files: `EDI/src/fast_ols_regression.cpp`, `EDI/src/fast_robust_regression.cpp`

Root cause: `X.T * X` or `X.T * diag(w) * X` computed as full GEMM, doing 2× unnecessary work for a symmetric result. `weighted_crossprod(X, ones)` or direct DSYRK call fills only the upper triangle; copy to lower at the end. Halves FLOPs for XtX computation.

---

**TODO-28: Nonparametric — Wilcox rank-sum O(n²) → O(n log n)** ☐
File: Wilcoxon source
Profiler: 16% branch-miss rate for wilcox_hl (highest of any kernel), consistent with O(n²) double loop over all pairwise comparisons. Hulsen-Lehmann estimator currently O(n²); merge-sort based U-statistic counting is O(n log n) (same algorithm as in numpy).

---

**TODO-29: Nonparametric — Ridit: std::map → std::unordered_map** ☐
File: Ridit source
`std::map<int, double>` for frequency table lookup is O(log n) per lookup. `std::unordered_map<int, double>` is O(1) amortized. For discrete count data where the map is built once and queried n times, this eliminates n × O(log n) lookups.

---

**TODO-30: Survival — coxph_var per-observation allocation** ☐
File: `EDI/src/fast_coxph_regression.cpp`
Allocates working matrices inside the observation loop in the variance computation path. Preallocate scratch outside the loop.

---

**TODO-31: Nonparametric — JT test: std::map → vector + precomputed table** ☐
File: Jonckheere-Terpstra source
Same O(log n) map lookup replaced with O(1) vector index.

---

**TODO-32: Survival — logrank extra O(n) passes** ☐
File: logrank source
Multiple O(n) traversals that can be fused into a single pass.

---

**TODO-33: Robust regression — MAD: sort → nth_element** ☐
File: `EDI/src/fast_robust_regression.cpp`
Median of absolute deviations currently sorts the residual array O(n log n). `std::nth_element` gives O(n) for just the median — half the work.

---

**TODO-34: G-computation — ordinal model: cache hessian across calls** ☐
File: g-computation source
Hessian is computed 3× in a single inference pass when once suffices. Cache and reuse.

---

**TODO-35: ZOIB — std::lgamma + fast_digamma** ☐
File: `EDI/src/fast_zoib.cpp`
Uses `R::lgammafn` and `R::digamma` where `std::lgamma` and `fast_digamma` apply. Same pattern as TODO-14 (beta) and TODO-24.

---

**TODO-36: Stereotype logit — hessian allocs** ☐
File: stereotype logit source
Hessian allocates intermediate matrices per call that can be preallocated.

---

**TODO-37: Poisson — cache delta_eta for IRLS step-halving** ☐
File: `EDI/src/fast_poisson_regression.cpp:236-248`
Each backtracking probe recomputes `X * beta_new` as a full GEMV. Precompute `delta_eta = X_free * direction` once; each probe becomes O(n) vector-add. Low priority: IRLS typically converges in 1-2 steps for well-conditioned Poisson data.

---

⚠️ **Ordinal family (prop_odds, adj_cat, cont_ratio, ord_cauchit, ord_cloglog, ord_probit) — findings PENDING** (agent a421d2a07a6313a26 still running as of 2026-07-02). Will add TODO-38+ when complete.
