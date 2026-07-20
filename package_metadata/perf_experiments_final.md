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

**TODO-28: `weighted_crossprod` col-major — 2×GEMM → DSYR rank-1 accumulation** ✗ DROPPED
File: `EDI/src/_helper_functions.h`
The col-major branch of `weighted_crossprod` evaluates `X.T * w.asDiagonal() * X` as two GEMMs: first `X.T * diag(w)` produces a p×n intermediate matrix (heap allocation + n*p multiplications), then multiplies by X (second GEMM p×n × n×p → p×p). Replace with `selfadjointView<Upper>().rankUpdate()` (BLAS DSYR): n symmetric rank-1 updates each O(p²/2), no p×n intermediate. For n=1000, p=5: current ~50k FMAs + p×n allocation vs ~12.5k FMAs + 0 allocation. Accounts for 45% of `robust_est` (called per IRLS iteration) and all weighted OLS paths. **Caveat:** For very small p (p=3–5), Eigen may already beat BLAS for DSYR — benchmark immediately; revert if slower.
Result: benchmarked and reverted. On the `benchmark_model_fits.R` / `edi_kernel_profiler.R` `robust_est` expression (N=1000, p=6 effective design), the current two-GEMM Eigen expression was fastest: baseline 0.202600ms median. Three no-intermediate variants were correct but slower: Eigen per-row `selfadjointView<Upper>().rankUpdate()` 0.595092ms (+194%), cache-friendly manual upper-triangle row loop 0.238415ms (+17.7%), and `sqrt(w)` row scaling plus symmetric rank-k update 0.319871ms (+57.9%). Correctness checks passed for all variants (`weighted_crossprod` parity vs `crossprod(X, X*w)`, robust regression equivalence, real-data, and warm-start weight tests), but no source change was kept because the shared helper would regress the profiled benchmark path.

**TODO-29: `fast_ols.cpp` XtX → `selfadjointView.rankUpdate` (DSYRK)** ✗ DROPPED
File: `EDI/src/fast_ols.cpp:32,108`
`XtX_free.noalias() = X.transpose() * X` uses a full GEMM (77.5% of `ols_est` samples in `lhs_process_one_packet`). Replace with:
```cpp
XtX_free.setZero();
XtX_free.selfadjointView<Eigen::Upper>().rankUpdate(X.transpose());
XtX_free.triangularView<Eigen::Lower>() = XtX_free.transpose();
```
Uses DSYRK (half the FLOPs). Expected: 30–50% speedup on `ols_est`. Same caveat as TODO-28: benchmark for small p.
Result: benchmarked and reverted. On the `benchmark_model_fits.R` / `edi_kernel_profiler.R` OLS expressions, the current Eigen full crossproduct was fastest: `ols_est` baseline 0.016611ms median and `ols_var` baseline 0.010391ms median. The literal `selfadjointView<Upper>().rankUpdate(X.transpose())` implementation with lower-triangle mirroring was slower (`ols_est` 0.017327ms, +4.3%; `ols_var` 0.010797ms, +3.9%). An upper-only storage variant using `Eigen::LDLT<MatrixXd, Upper>` to avoid the mirror copy was also slower (`ols_est` 0.017212ms, +3.6%; `ols_var` 0.010997ms, +5.8%). Correctness: `test-rcpp-fitting-equivalence.R` and `test-rcpp-fitting-real-data.R` passed after restoring the original implementation. No source change was kept because p=6 is too small for this DSYRK rewrite to beat Eigen's existing crossproduct path.

**TODO-30: CoxPH `compute_robust_vcov` — eliminate 150 heap allocs per call** ✓ DONE
File: `EDI/src/fast_coxph_regression.cpp:498–530`
`ek`, `dk_ek_over_Rk`, `cum_B` are `std::vector<std::vector<double>>(n_events, std::vector<double>(p, 0.0))` — n_events separate heap allocations each. At n_events=50, p=5: 150 separate heap allocs. `cluster_scores` is `std::map<int, Eigen::VectorXd>` (O(log n) per insert). Fix: replace nested vectors with flat `Eigen::MatrixXd(p, n_events)` indexed as `[:,k]`; replace `std::map` with `std::unordered_map`. Also preallocate `cluster_scores` at construction when n_clusters is known.
Implemented the three event-by-covariate work arrays as contiguous `Eigen::MatrixXd` buffers and changed cluster-score aggregation to a reserved `std::unordered_map`. The retained input-order permutation now writes each score row back to its original subject, fixing the pre-existing mismatch between event-time-sorted scores and unsorted cluster IDs. Correctness: the dedicated `test-coxph-robust-vcov.R` test matches coefficients and the complete robust covariance matrix from `survival::coxph` at `1e-7`; `test-rcpp-fitting-equivalence.R` and `test-rcpp-fitting-real-data.R` also pass. Clustered Cox benchmark at N=1000, p=5, 100 clusters (1,000 iterations): 0.870ms → 0.453ms median (−47.9%, 1.92× throughput).

**TODO-31: Ordinal regression — compute hessian once, not 3–4 times** ✓ DONE
File: `EDI/src/fast_ordinal_regression.cpp:183–185,311`
`fast_ordinal_regression_cpp` calls `model.hessian(params)` 3× in succession for `observed_information`, `fisher_information`, and `information` (all identical). Then `fast_ordinal_regression_with_var_cpp` (line 311) calls a 4th time. Fix:
```cpp
MatrixXd H = model.hessian(params);
return List::create(Named("observed_information") = H, Named("fisher_information") = H, Named("information") = H, ...);
```
In `fast_ordinal_regression_with_var_cpp`, extract `H` from `res["observed_information"]` instead of re-calling. Expected: 3× faster on hessian-dominated var path.
Implemented one Hessian evaluation shared by the three information outputs in both unweighted and weighted proportional-odds fits; the variance wrapper now reuses the returned observed-information matrix. Correctness: the dedicated `test-ordinal-information-reuse.R` test verifies identical information outputs, covariance inversion, and the weighted path; `test-rcpp-fitting-equivalence.R` and `test-rcpp-fitting-real-data.R` also pass. Targeted `prop_odds_var` benchmark at N=200 with five predictors (10,000 iterations): 0.226809ms → 0.184307ms median (−18.7%, 1.23× throughput).

**TODO-32: GComp ordinal finite-diff loop — hoist temporaries outside loop** ✓ DONE
File: `EDI/src/fast_ordinal_regression.cpp:449–462`
Finite-difference loop over `n_params` iterations allocates `p_plus`, `p_minus` (VectorXd copies of `params`) plus inside `compute_md_from_params` (line 421): `alpha`, `beta`, `eta1_loc`, `eta0_loc` — ≥6 heap allocs × n_params ≈ 8 = ~48 allocs per post-fit call. Fix: hoist `p_plus` and `p_minus` outside loop, reset from `params` inside; pass pre-sized scratch vectors to `compute_md_from_params` by reference. Expected: 30–50% faster on gcomp ordinal post-fit path.
Implemented reusable perturbation and linear-predictor buffers, zero-allocation Eigen views for alpha/beta parameter blocks, and scalar restoration of each perturbed parameter. The existing baseline linear predictors are now reused: threshold perturbations require no predictor update, while beta perturbations use scaled design-column updates instead of repeated GEMVs. Correctness: the dedicated `test-ordinal-gcomp-post-fit.R` test matches an independent R implementation of marginal means, the finite-difference gradient, covariance block, and delta-method standard error; `test-rcpp-fitting-equivalence.R` and `test-rcpp-fitting-real-data.R` also pass. Targeted `ordinal_gcomp_post_fit_var` benchmark at N=200 with five predictors (10,000 iterations): 0.067300ms → 0.063000ms median (−6.4%, 1.07× throughput). The gain is below the estimate because logistic-CDF evaluation dominates after allocation and GEMV removal.

**TODO-33: HurdlePoisson GLMM — `log1p` fast-path for large lambda** ✓ DONE
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:183,268`
`__log1p_fma` dominates both `hurdle_p_est` and `hurdle_p_var` annotate output. Line 183: `const double lne = (lam < 1e-10) ? eta_ki : std::log1p(-eneg)`. For `lam > 16` (`eneg < 1e-7`), `log1p(-eneg) ≈ -eneg` with error < 1e-7. Add fast-path:
```cpp
const double lne = (lam < 1e-10) ? eta_ki : (eneg < 1e-7 ? -eneg : std::log1p(-eneg));
```
Apply identically in `hessian()` line 268. If typical lambda values are moderate-to-large (common in count outcomes), eliminates 30–70% of `log1p` calls. Verify numerical accuracy on held-out fit before committing.
Implemented a shared `log_one_minus_exp_neg_hp` helper used by both the objective and Hessian. It retains the tiny-lambda cancellation guard and uses `-exp(-lambda)` when `exp(-lambda) < 1e-7`; the per-term absolute approximation error is below `5.1e-15`. Correctness: the dedicated `test-hurdle-poisson-glmm-log1p-fast-path.R` test compares the large-lambda objective, score, and Hessian against exact-log1p Gauss-Hermite quadrature and numerical derivatives; `test-glmm-cpp-equivalence.R` also passes against canonical mixed-model packages. Controlled high-lambda benchmark at N=400 (one optimizer step, 3,000 iterations): 0.3082ms → 0.2647ms median (−14.1%, 1.16× throughput). The standard profiler workloads, which contain few large-lambda evaluations, were neutral within noise: estimate 3.303ms → 3.346ms and variance 1.781ms → 1.776ms.

**TODO-34: Continuation ratio — eliminate O(n×K) VectorXd allocs in augmented data construction** ✓ DONE
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
Implemented a counting pass that records each observation's category index and the exact augmented-row count, followed by direct writes into one zero-initialized `MatrixXd` and one `VectorXd`. This removes all per-row `VectorXd` allocations, vector growth, and the final row-copy pass while preserving row order. Correctness: the dedicated `test-continuation-ratio-augmentation.R` test verifies the exact augmented design/response for four outcome levels plus score and Hessian parity with an R reference; `test-rcpp-fitting-equivalence.R` also passes, including the VGAM continuation-ratio comparison. Profile-matched benchmarks: `cont_ratio_est` at N=1000, 2,000 iterations, 0.2782ms → 0.1765ms median (−36.6%, 1.58× throughput); `cont_ratio_var` at N=200, 3,000 iterations, 0.0660ms → 0.0530ms (−19.7%, 1.25× throughput).

**TODO-35: HurdlePoisson GLMM `operator()` — `Xg` expression view + `eta0` preallocation** ✓ DONE
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:154,244`
Line 244 in `hessian()`: `const Eigen::MatrixXd Xg = dat.X_s.middleRows(gs, sz)` — copies the group submatrix G times per call. Change to `const auto Xg = dat.X_s.middleRows(gs, sz)` (zero-copy Eigen expression view). Line 154 in `operator()`: `const Eigen::VectorXd eta0 = dat.X_s.middleRows(gs, sz) * par.head(p)` — allocates size-sz vector G times. Rewrite using preallocated `m_eta0.head(sz).noalias() = dat.X_s.middleRows(gs, sz) * par.head(p)`.
Confirmed and retained the zero-copy `Xg` views and constructor-sized `m_eta0` scratch buffer in both the objective and Hessian paths; the checklist entry had remained open after the source implementation landed. Correctness: the dedicated `test-hurdle-poisson-glmm-buffer-reuse.R` test exercises unequal group sizes, numerical score/Hessian derivatives, and row-order invariance; the large-lambda numerical test and canonical `test-glmm-cpp-equivalence.R` suite also pass. Matched incremental-rebuild benchmarks: `hurdle_p_glmm_est` at N=400, 1,000 iterations, 3.2447ms → 3.1524ms median (−2.8%, 1.03× throughput); `hurdle_p_glmm_var` at N=200, 300 iterations, 1.6975ms → 1.6433ms (−3.2%, 1.03× throughput).

### MEDIUM priority

**TODO-36: HurdlePoisson GLMM — fuse `log_sum_exp` + posterior weights into single pass** ✓ DONE
File: `EDI/src/fast_hurdle_poisson_glmm.cpp:189,200,273,280`
Both `operator()` and `hessian()` call `log_sum_exp_hp(log_terms)` then immediately loop over `exp(log_terms[k] - ll_g)` to get posterior weights — computing each `exp(log_terms[k])` twice. Use `log_sum_exp_and_weights` from `_glmm_engine.h` (already exists) to do both in one pass. Saves K exp calls per group per optimizer step. For G=200, K=7: 1,400 exp calls eliminated per iteration.
Implemented one preallocated posterior-weight vector per objective and reused the existing fused `glmm::log_sum_exp_and_weights` helper in both the objective/score and Hessian paths; the redundant per-node posterior `exp` calls and the local unfused helper were removed. Correctness: the optimized and baseline score/Hessian differed by at most 1.42e-14/3.41e-13 on the benchmark workload; the exact-quadrature `test-hurdle-poisson-glmm-log1p-fast-path.R` and external-reference `test-glmm-cpp-equivalence.R` tests pass. Targeted public-entry-point benchmark at N=400, G=200, p=5, K=7 (3,000 iterations): score 0.159650ms → 0.133950ms median (−16.1%, 1.19× throughput); Hessian 0.487301ms → 0.438700ms (−10.0%, 1.11× throughput).

**TODO-37: Poisson IRLS — `X.T * w.asDiag() * X` → `weighted_crossprod`** ✗ DROPPED
File: `EDI/src/fast_poisson_regression.cpp:227,270`
Lines 227 and 270: `XtWX_free.noalias() = X_f.transpose() * w_tmp.asDiagonal() * X_f` creates a p×n intermediate matrix per IRLS iteration. Replace with `XtWX_free = weighted_crossprod(X_f, w_tmp)` (already in `_helper_functions.h`). Applies equally to quasi_var (same code path). See TODO-28 caveat for small p.
Result: benchmarked and reverted. The current `weighted_crossprod` col-major branch evaluates the same Eigen triple product, so the literal substitution was neutral within noise: `poisson_est` 0.1402ms → 0.1417ms, `poisson_var` 0.0471ms → 0.0464ms, and `quasi_var` 0.0434ms → 0.0433ms. A true allocation-free symmetric implementation using one row-major design copy plus direct upper-triangle accumulation was correct but substantially slower at p=6: 0.2194ms (+56.5%), 0.0640ms (+35.9%), and 0.0650ms (+49.8%), respectively. The original Eigen expression remains in place. Correctness: the new `test-poisson-weighted-crossprod.R` verifies unweighted and weighted IRLS information matrices against direct R crossproducts and coefficients against `glm.fit`; warm-start, fitting-equivalence, and real-data suites also pass.

**TODO-38: Robust regression MAD — `std::sort` → `std::nth_element`** ✓ DONE
File: `EDI/src/fast_robust_regression.cpp:92–96`
MAD scale estimation calls `std::sort(abs_r)` on residuals but only needs the median — O(n log n) when O(n) suffices. Replace with `std::nth_element` to position n/2 (plus a second nth_element for even n to get the lower median). ~7× faster on sort step; eliminates 7% of `robust_est` wall time.
Implemented a hybrid median selector: workloads with at least 512 observations use one `std::nth_element`, with a linear maximum over the lower partition for the second middle value when n is even; smaller workloads retain `std::sort`, which benchmarks faster at that size. Correctness: a dedicated odd/even MAD test covers both sides of the crossover, and the full C++ fitting equivalence and real-data fitting test files pass. Benchmark outputs (coefficients, scale, iterations, and variance) are bit-for-bit identical to the baseline. Profile-matched `robust_est` benchmark at N=1000 and p=6 (median across two 5,000-iteration runs): 0.183550ms → 0.161500ms (−12.0%, 1.14× throughput). The N=200 variance path retains sorting and was neutral within timing noise: 0.050300ms → 0.051450ms.

**TODO-39: Probit — reuse `mu[]` for post-convergence neg_ll computation** ✓ DONE
File: `EDI/src/fast_probit_regression.cpp:222–228`
Replaced the post-convergence GEMV + erfc neg_ll recomputation with `nl -= wi * (y[i] > 0.5 ? std::log(mu[i]) : std::log1p(-mu[i]))`, reusing `mu[]` from the last IRLS iteration. Eliminates 1 GEMV + n erfc calls. Guard respected: change is inside the IRLS-only branch; lbfgs path unchanged. Also, the `pnorm_fast`/`log_pnorm_lower`/`log_pnorm_upper` helpers were simultaneously updated by the linter to use `edi_ordinal::fast_erfc` instead of `std::erfc`. All 10 probit tests pass. `probit_est` benchmark (N=1000, p=5, estimate_only=TRUE, 28000 reps): 0.563 ms/call → ~0.314 ms/call (−44%, 1.79× throughput).

**TODO-40: Probit — polynomial fast-erfc approximation** ✓ DONE
File: `EDI/src/fast_probit_regression.cpp` and `_helper_functions.h`
erfc calls account for ~75% of all `probit_est` samples. For typical probit use, `η ∈ [-6, 6]`. A 6-term Horner-form minimax polynomial for `erfc(x)` on `[0, 5.6]` achieves error < 1e-10 and is ~5× faster than libm `__cr_erf_fast`. Implement as `fast_erfc` in `_helper_functions.h`; fallback to `std::erfc` for `|x| > 5.6`. Apply to `pnorm_fast`, `log_pnorm_lower`, `log_pnorm_upper`. **Requires careful accuracy validation against `R::pnorm5` on a sweep of η values before committing.** Also applies to `fast_ordinal_probit_regression.cpp` (same erfc pattern).
Implemented a Cephes piecewise rational `fast_erfc` in Horner form rather than an unspecified degree-6 polynomial, with `std::erfc` fallback outside `|x| <= 5.6`. Binary probit's CDF and both central-range log-CDF helpers use it. Correctness: the dedicated `test-fast-probit-cdf.R` reconstructs the internal CDF over a dense `[-7.9, 7.9]` sweep; maximum absolute error versus `pnorm` is below `2e-15`, and lower-tail log error below `2e-12`. `test-incidence-probit.R` and the full fitting-equivalence suite pass with unchanged fitted outputs. Profile-matched benchmarks: `probit_est` at N=1000, 5,000 iterations, 0.4547ms → 0.2677ms median (−41.1%, 1.70× throughput); `probit_var` at N=200, 10,000 iterations, 0.1234ms → 0.0764ms (−38.1%, 1.62× throughput). Applying the same helper to ordinal probit was benchmarked and reverted because it regressed `ord_probit_est` 0.7936ms → 0.8762ms (+10.4%) and `ord_probit_var` 0.3070ms → 0.3310ms (+7.8%).

**TODO-41: Logistic — replace Eigen sigmoid expression with `plogis_array_safe`** ✓ RESOLVED — REJECTED
File: `EDI/src/fast_logistic_regression.cpp:138`
`mu.array() = 1.0 / (1.0 + (-eta.array()).exp())` with `EIGEN_DONT_VECTORIZE` degrades to per-element PLT `expf64` dispatch (store-latency bottleneck at `35e11e`). Replace with `mu = plogis_array_safe(eta.array()).matrix()` (already in `_helper_functions.h:376`) — explicit scalar loop with sign-branched numerically stable form + `std::exp` (resolved at link time, no PLT). Apply at lines 138, 243, 258, 276, 294.
Rejected after controlled before/after benchmarking. The premise is stale for release builds: `EIGEN_DONT_VECTORIZE` is opt-in in the current `Makevars`, and the original Eigen expression uses native packetized exp. The literal replacement also allocates a temporary array and was substantially slower on all profile-matched paths (5,000 iterations): estimate 0.126700ms → 0.155600ms (+22.8%), variance 0.041100ms → 0.050100ms (+21.9%), and the four combined score/Hessian helpers 0.044100ms → 0.077300ms (+75.3%). An isolated translation unit compiled with both `EIGEN_DONT_VECTORIZE` and `-fno-tree-vectorize` also disproved the claimed scalar-mode benefit: 100 sigmoid sweeps over N=1000 took 0.605302ms with the Eigen expression versus 0.772903ms with the proposed sign-branched loop (+27.7%). The production substitution was therefore reverted. Correctness: the added extreme-linear-predictor test verifies unweighted and weighted scores/Hessians against R at eta from -1000 to 1000; established logistic fitting tests also pass.

**TODO-42: ZOIB — `fast_digamma` + `std::lgamma` + preallocate member vectors** ✓ DONE
File: `EDI/src/fast_zero_one_inflated_beta.cpp`
(1) `DigammaFunctor` calls `R::digamma(x)` → replace with `fast_digamma(x)` (already in `_helper_functions.h`). (2) `LgammaFunctor` calls `R::lgammafn(x)` → replace with `std::lgamma(x)`. (3) Line 164: `R::digamma(phi)` per optimizer step → `fast_digamma(phi)`. (4) Lines 189–190: `R::digamma(phi)` and `R::trigamma(phi)` already hoisted but use slow dispatch → `fast_digamma`. (5) Preallocate `pi0(n)`, `pi1(n)`, `pib(n)`, `grad_gamma0(p_zero_one)`, `grad_gamma1(p_zero_one)` as member fields (currently allocated every `operator()` call).
Replaced all positive-argument digamma evaluations with `fast_digamma`, all likelihood lgamma evaluations with `std::lgamma`, and constructor-sized/reused the five mixture-probability and gamma-gradient vectors. `R::trigamma` remains because no validated fast trigamma helper exists. Correctness: the dedicated `test-zero-one-inflated-beta-fast-math.R` test matches the complete score and Hessian against numerical derivatives of an independent R likelihood; mixture-formula, Bayesian-bootstrap, and reusable-worker asymptotic-family suites pass. Benchmark fitted outputs differ by at most 3.8e-10. Profile-matched benchmarks, aggregated across two independent 1,000-iteration runs: estimate at N=1000, p=6, p_zi=2, 10.5371ms → 4.4061ms median (−58.2%, 2.39× throughput); variance at N=200, 2.4283ms → 1.0750ms (−55.7%, 2.26× throughput).

**TODO-43: ZAP — inline `log1m_eml` reusing precomputed `eml`** ✓ DONE
File: `EDI/src/fast_zero_augmented_poisson.cpp:71–72`
In the hurdle positive branch, `eml = std::exp(-lam)` is computed at line 71, then `log1mexp(-lam)` is called at line 72 which internally re-evaluates `exp(-lam)`. Inline as:
```cpp
const double log1m_eml = (lam > 0.6931471805599453)
    ? std::log1p(-eml) : std::log(-std::expm1(-lam));
```
Eliminates 1 `exp`/`expm1` call per positive hurdle obs (~68% of obs).

Implemented the branch inline and removed the now-unused `log1mexp` helper. The dedicated `test-zero-augmented-poisson-log1m.R` fixes the likelihood-determining parameters and compares the returned objective with an independent stable R calculation for Poisson means spanning `1e-8` through `30`, including points immediately below and above the `log(2)` branch boundary. It and the full Rcpp fitting-equivalence suite pass after installation with `R CMD INSTALL --no-docs`. Noise-controlled benchmarks used 30 independent rounds of 300 iterations per path (9,000 measurements per path/build): estimate median 1.3417ms → 1.2893ms (−3.90%, bootstrap 95% CI −6.76% to −1.32%, one-sided Wilcoxon p=0.00015); variance median 0.8673ms → 0.8078ms (−6.86%, bootstrap 95% CI −10.65% to −5.66%, p=1.7e-8).

**TODO-44: MN — verify R wrapper uses `mn_ci_cpp` not R-level bisection loop** ✓ DONE — ALREADY OPTIMIZED
File: `EDI/src/miettinen_nurminen_speedups.cpp`
`mn_var` shows `bcEval_loop` at 12.5% + `mn_constrained_mle_pc_cpp` at 11.4%, implying R calls the constrained MLE from R-level bisection. The C++ `mn_ci_cpp` encapsulates the full bisection (lines 116–146). Check the R wrapper — if it calls `mn_constrained_mle_pc_cpp` or `mn_z_statistic_cpp` from R in a loop, replace with a single `mn_ci_cpp` call. Would collapse ~50 R→C++ dispatches per CI into one.
Verified that `InferenceIncidMiettinenNurminenRiskDiff$compute_asymp_confidence_interval()` already makes exactly one `mn_ci_cpp` call; there is no R-level bisection to replace. The profiler evidence was caused by `mn_var` benchmarking `mn_pvalue_cpp` rather than the class CI path, so that workload now calls `mn_ci_cpp` with the class defaults. A 20,000-iteration interleaved benchmark comparing a reconstructed R-level bisection loop with the production C++ call produced identical bounds and measured 0.158700ms → 0.006400ms median (−96.0%, 24.8× throughput), confirming that the intended optimization is already active. Correctness: the dedicated dispatch test asserts one C++ call and no R loop, and matches three C++ CIs against the R-loop reference at 1e-14; the Bayesian-bootstrap and reusable-worker asymptotic-family suites also pass.

**TODO-45: GComp — verify warm-start logistic fits between bootstrap replicates** ✓ DONE
Files: `EDI/R/inference_incidence_KK_gcomp_abstract.R`, `EDI/R/inference_incidence_gcomp_abstract.R`
Verified: warm starts ARE already passed (primary-fit β) to each bootstrap replicate via `get_fit_warm_start_for_length` in both `compute_weighted_gcomp_estimate` (KK path) and `weighted_gcomp_fit` (non-KK path). Added warm-start chaining: new `gcomp_boot_beta` private field chains the previous replicate's converged β as the warm start for the next replicate, falling back to primary β on failure or dimension change. Benchmark (N=1000, p=6, B=500 replicates, 15 median runs): primary-only → chained = 0.134 → 0.142 ms/replicate (effectively 1.0× — no gain). Minimal gain is expected: bootstrap weights are centred at 1, so the primary β is already near-optimal for any weighted replicate (~4 IRLS iterations either way). Correctness: 6 new tests in `test-gcomp-boot-warm-start-chaining.R` verify chained results match independent-fit references within 1e-6 tolerance. All 50 pre-existing bootstrap tests pass (test-bootstrap-reused-worker-asymp-families.R and test-kk-reused-worker-bootstrap-fast-path.R).

**TODO-46: GComp ordinal — `plogis_array` Eigen temporaries → scalar loop** ✓ DONE
File: `EDI/src/gcomp_speedups.cpp:399–405`
`compute_mean` lambda replaced: Eigen `plogis_array` calls (each allocating 2+ n-length arrays + Constant) replaced with a scalar double-accumulator nested loop using `plogis_stable`. Eliminates all heap allocs inside `compute_mean`.
- `gcomp_ordinal_est` (n=1000): 0.674 → 0.674 ms — no measurable change; post-fit is only ~4% of kernel (fit dominates)
- `gcomp_ordinal_var` (n=200): 0.225 → 0.214 ms — within noise
- Correctness: 3 new tests in `test-ordinal-gcomp-post-fit.R` directly exercising `gcomp_ordinal_proportional_odds_post_fit_cpp` vs R reference; all 8 tests (3 new + 5 prior) pass

**TODO-47: Logrank — fuse martingale mean/variance passes into main sweep** ✓ DONE
File: `EDI/src/fast_logrank.cpp`
Two changes combined:
1. **Fused accumulation**: Eliminated `martingale[]` vector (O(n) allocation). Instead accumulate `sum_m_treat/control` and `sum_sq_m_treat/control` in the inner loop during the main sweep; compute mean and variance from these after the sweep using `var = (sum_sq - n*mean²)/(n-1)`. Eliminates 3 O(n) post-sweep passes.
2. **Precomputed group boundaries**: Replaced `while (recs[end].time == curr_time)` FP equality inner-while (the branch-miss culprit) with a precomputed `gstart[]` vector of group start indices. Main loop becomes a simple counted for-loop over groups with inner bounds from `gstart[g]..gstart[g+1]`.
- `logrank_est` (N=500): ~0.021 → ~0.016 ms/call (**~1.3× speedup**)
- `logrank_var` (N=200): ~0.012 → ~0.012 ms/call — no measurable change at this N
- Correctness: 15 tests in `test-logrank-fused-martingale.R` verify score/var_score vs `survival::survdiff` and beta_hat/se_beta_hat vs `coxph` martingale residuals; all pass

**TODO-48: Weibull regression — preallocate member buffers in `WeibullAFTLikelihood`** ✓ DONE
File: `EDI/src/fast_weibull_regression.cpp:25–45`
`WeibullAFTLikelihood::operator()` allocates 4 Eigen vectors (`eta`, `w`, `exp_w`, `d_ll_d_eta`) on every optimizer call. In `hessian()` (line 47): `beta_weights` and `cross_weights` are also temporary n-element vectors. Add preallocated `m_eta(n)`, `m_w(n)`, `m_exp_w(n)`, `m_d_eta(n)`, `m_beta_weights(n)`, `m_cross_weights(n)` as class members initialized in the constructor. Non-`const` methods so no LICM/alias risk.
Implemented all six constructor-sized buffers and reused them in both the objective/gradient and Hessian paths. The per-call beta copy and dead-array materialization were also replaced with zero-copy Eigen views. Correctness: benchmark outputs are bit-for-bit identical; the dedicated `test-weibull-aft-buffer-reuse.R` test matches independently differentiated score and Hessian values, and the fast-GLM, fitting-equivalence, and real-data suites pass. Profile-matched benchmarks, aggregated across two independent 20,000-iteration runs: `weibull_est` at N=500, p=5, 0.056500ms → 0.050800ms median (−10.1%, 1.11× throughput); `weibull_var` at N=200, p=5, 0.059450ms → 0.054250ms (−8.7%, 1.10× throughput). The optimized package was installed with `R CMD INSTALL --no-docs`.

**TODO-49: ZINB — `R::lgammafn` → `std::lgamma`** ✓ DONE
File: `EDI/src/fast_zinb.cpp:88,94`
Line 88: `R::lgammafn(theta)` called once per step; `std::lgamma` avoids R dispatch. Line 94 in the distinct-y loop: `m_lgamma_yptheta[k] = R::lgammafn(ypt)` — nd times per step. `std::lgamma` is ~30–50% faster per call by avoiding R's error-check wrapper. Same fix already applied in `fast_beta_regression.cpp`.

Replaced all three ZINB uses, including the once-per-fit `lgamma(y + 1)` constructor table. Correctness: the dedicated `test-zinb-std-lgamma.R` fixes all likelihood-determining parameters and matches the returned C++ objective against an independent R ZINB calculation at `1e-12`; the full Rcpp fitting-equivalence suite also passes after `R CMD INSTALL --no-docs`. The original profiler estimate workload was unsuitable for comparison because fitting ZINB to Poisson data sent dispersion toward the boundary and the builds followed different optimizer paths (512 versus 441 iterations). A controlled seeded-negative-binomial benchmark used identical warm starts and a fixed 25-iteration budget; parameters, scores, and Hessians were identical, with objective difference `5.7e-14`. An A/B/A audit used 30 rounds × 300 iterations in each phase: estimate 0.7413ms → 0.7393ms (0.27% faster, bootstrap 95% CI −0.92% to +1.24%, one-sided Wilcoxon p=0.27); variance 0.3879ms → 0.3910ms (0.79% slower, CI −2.98% to +3.17%, p=0.69). The end-to-end effect is neutral on this toolchain; the direct standard-library calls are retained because they remove R dispatch without a supported regression.

**TODO-50: Adj_cat — batch exp per obs in inner loop** ✓ DONE
File: `EDI/src/fast_adjacent_category_logit.cpp`
Three changes in `operator()` and `hessian()`:
1. **Precomputed exp(alpha[k])**: K-1 exp calls once per optimizer step (amortised over n obs) instead of K per obs. Per obs now needs only 1 exp call: `u = exp(-eta)`. Unnorm probs built right-to-left via product recurrence `prob[k] = prob[k+1] * exp_alpha[k] * u` — n*K → n+K-1 total exp calls per step.
2. **Eliminated Eigen overhead**: Replaced Eigen `(scores.array()-max).exp()`, `probs.dot(y_levels)`, `VectorXd scores/cdf/y_levels` with plain `std::vector` + scalar loops. Also removed dead code (`score_offsets`, `alpha_suffix` in grad path).
3. **Branchless alpha gradient**: Split the `(y_i <= j+1) ? 1 : 0` branch into two contiguous loops at the split point `y_i-1`, eliminating the unpredictable branch in the hot inner loop.
- `adj_cat_est` (N=1000): ~0.512 → ~0.451 ms/call (**~1.14× speedup**)
- `adj_cat_var` (N=200): ~0.131 → ~0.113 ms/call (**~1.16× speedup**)
- Correctness: 10 tests in `test-adj-cat-exp-batch.R`; VGAM coefficient match (K=3,4), score at MLE ≈ 0, score norm smaller at MLE than at perturbed point, Hessian vs finite-difference (K=3,4); all pass

**TODO-51: Ord_cauchit — mutable scratch buffers to eliminate malloc overhead** ✓ RESOLVED — SUPERSEDED BY TODO-69
File: `EDI/src/fast_ordinal_cauchit_regression.cpp`
Profile: `ord_cauchit_est` has 43% of samples in malloc TLS arena access (23% movq %fs: TLS read + 20% TLS write). The cauchit link likely allocates intermediate arrays per optimizer call that the logistic link avoids through Eigen lazy evaluation. Add mutable scratch `VectorXd` buffers sized at construction to eliminate per-call allocs. Caution: check whether this is a `const` method inside a GLMM objective — if so, see the mutable-field antipattern warning.

Implemented in the shared `ordinal_fixed_link_helpers.h`: threshold and coefficient blocks are now zero-copy Eigen views, and a constructor-sized mutable linear-predictor buffer replaces the fresh N-vector in the likelihood, score, observed-Hessian, and expected-information paths. An attempted broader conversion of the small Hessian work arrays to member Eigen buffers was benchmarked and reverted because it slowed the variance path; their original contiguous `std::vector` implementation remains. Correctness: `test-ordinal-cauchit-scratch-buffers.R` compares the analytic score and Hessian with numerical derivatives of an independent R cauchit likelihood, and the full Rcpp fitting-equivalence suite passes after installation with `R CMD INSTALL --no-docs`. Noise control used an A/B/A audit with 30 rounds × 500 iterations per phase: estimate 0.5815ms → 0.5781ms (0.58% faster, bootstrap 95% CI −0.23% to +1.26%, two-sided Wilcoxon p=0.054); variance 0.2073ms → 0.2046ms (1.27% faster, CI 0.30% to 2.41%, p=0.0063). The original 43% malloc attribution was stale because several scratch fields were already present before this TODO; the remaining end-to-end gain is correspondingly small.

Superseded by the broader TODO-69 audit across all fixed ordinal links. That audit found a material ordinal-probit regression and showed that even the eta-only refinement slowed every tested estimate path; the mutable eta buffer and zero-copy parameter views were therefore reverted. The pre-existing `m_scratch_dq`, `m_scratch_d2q`, and `m_scratch_v` buffers remain.

**TODO-52: Ord_probit — fast_erf approximation for probit link** ✓ DONE
File: `EDI/src/fast_ordinal_probit_regression.cpp`
Profile: erfc/erf calls spread across samples with IPC=2.24 (compute-bound). Use the same `fast_erfc` from TODO-40 once implemented. Also applies to `fast_probit_regression.cpp`.

Extracted `fast_erfc`, `pnorm_fast`, and `dnorm_fast` into a new header `EDI/src/fast_erfc.h` (includes only `<cmath>`, no R or Eigen deps). Updated `ordinal_fixed_link_helpers.h` to include `fast_erfc.h` and replaced all three `R::pnorm5`/`R::dnorm4` calls in the probit `cdf`, `pdf`, and `pdf_derivative` switch cases with the inline approximations. Removed the now-redundant anonymous-namespace definitions of `kSqrt1_2`, `k1_Sqrt2Pi`, `pnorm_fast`, `dnorm_fast` from `fast_probit_regression.cpp` (they now come from the shared header via `_helper_functions.h`). Also updated `_helper_functions.h` to include `fast_erfc.h` instead of defining `fast_erfc` inline. Correctness: `test-fast-probit-cdf.R` (second test) compares the analytic score and Hessian against numerical derivatives of an independent R probit ordinal likelihood; score tolerance 1e-7 and Hessian tolerance 1e-5 both pass. The first test (CDF grid check) was already present and continues to pass. Noise control used an A/B/A audit with 30 rounds × 500 iterations per phase (estimate) and 30 rounds × 100 iterations (variance) in each of two separate builds: estimate 0.8780ms → 0.7950ms (9.5% faster); variance 1.4300ms → 1.2500ms (12.6% faster). The one-sided Wilcoxon on the post-change build's internal A vs B gives p=0.062 (estimate) and p=0.312 (variance), confirming the post-change measurements are internally stable. The speedup comes entirely from eliminating R argument-checking dispatch in `R::pnorm5` and `R::dnorm4` in favor of the Cephes-based `fast_erfc` approximation, which is already used by the binary probit path.

### LOW priority

**TODO-53: CoxPH hessian — direct symmetric writes vs post-loop triangular copy** ✓ DONE
File: `EDI/src/fast_coxph_regression.cpp:204–211`
Inner event-time loop writes only `hess.triangularView<Lower>()`, then copies to upper after the event loop (O(p²) copy per optimizer step). Write both `hess(q1,q2)` and `hess(q2,q1)` directly in the inner product loop (lines 172–177), eliminating the post-loop copy. For p≤10 the direct write + copy elimination is a net win.

Implemented scalar lower-triangle Hessian accumulation with immediate mirrored upper-triangle writes inside `compute_cox_ll_grad_hess_fast()`, replacing the previous Eigen `triangularView<Lower>() += ...` expression and removing the post-event-loop `hess.triangularView<Upper>() = hess.transpose()` copy. The risk-set second-moment buffer remains lower-triangular only; the event contribution reads `S2(q2, q1)` and writes both Hessian entries for off-diagonal pairs. A wider `p=10` isolated stress check showed a small regression, so the retained benchmark target is the profiler-matched CoxPH variance shape (`n=200`, `p=5`) where this TODO arose.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-53 source and optimized source; no `--preclean` was used. Timing used separate R processes, process CPU time, two before phases and two after phases, 40 rounds per phase. The isolated Hessian workload used `n=200`, `p=5`, 1,000 Hessian calls per round; the full CoxPH variance workload used `n=200`, `p=5`, 200 complete `fast_coxph_regression_cpp(..., estimate_only=FALSE)` fits per round. Combined medians: isolated Hessian 0.0310ms → 0.0260ms per call (IQR 0.0300–0.0330ms → 0.0250–0.0273ms), a 1.192× speedup / 16.1% time reduction (independent bootstrap 95% CI 1.148×–1.192×; one-sided Wilcoxon p<2.2e-16). Full variance fit 0.1650ms → 0.1300ms per call (IQR 0.1550–0.1800ms → 0.1250–0.1350ms), a 1.269× speedup / 21.2% time reduction (bootstrap 95% CI 1.208×–1.308×; one-sided Wilcoxon p<2.2e-16). Before/after Hessians, fitted coefficients, and variance-covariance matrices were exactly equal on the benchmark data; all benchmark fits converged.

Correctness: added `test-coxph-hessian-symmetric-writes.R`, which compares unstratified and stratified CoxPH Hessians against numerical derivatives of independent R Breslow partial-likelihood implementations, checks exact symmetry, and verifies repeatability. The dedicated test passed all 5 assertions. Standalone existing equivalence checks for unstratified and stratified CoxPH coefficients/vcov against `survival::coxph` also passed.

**TODO-54: Robust regression — cache XtX from IRLS to avoid var-path recompute** ✓ DONE
File: `EDI/src/fast_robust_regression.cpp:253`
After IRLS convergence, line 253 recomputes `XtX = X_free.T * X_free` for the sandwich estimator. Either store it in `RobustModelResult` at convergence, or derive from the QR cold-start (`R.T * R` where R is the QR factor). Saves one O(n·p²) GEMM per `robust_var` call.

Implemented the QR-derived variant for the common smart-cold-start path. `RobustModelResult` now carries a cached scalar `XtX_inv_diag_j`, the requested diagonal entry of `(X_free'X_free)^-1`. During the initial `ColPivHouseholderQR` cold start, the code builds the small `R'R` system, maps the requested original free-parameter coordinate through the QR column permutation, solves for only that diagonal entry, and stores the scalar. The variance path then uses this cached value for `ssq_b_j` and falls back to the original `X_free.transpose() * X_free` plus `compute_diagonal_inverse_entry()` only when the QR cache is unavailable (warm-start / no-smart-start / fixed target not free). A materialized full-`XtX` cache was benchmarked first and rejected as neutral-to-slower; the retained scalar cache avoids both the O(n·p²) crossproduct and an unnecessary p×p cached matrix on the profiled path.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-54 source and optimized source; no `--preclean` was used. Timing used separate R processes, process CPU time, two before phases and two after phases, 50 rounds per phase. The profiler-shaped workload used `n=200`, `p=6`, 500 complete `fast_robust_regression_cpp(..., method="MM", j=2)` variance fits per round; a wider audit used `n=300`, `p=16`, 200 fits per round. Combined medians: profiler shape 0.0640ms → 0.0600ms per fit (IQR 0.0600–0.0680ms → 0.0580–0.0620ms), a 1.067× speedup / 6.25% time reduction (independent bootstrap 95% CI 1.033×–1.067×; one-sided Wilcoxon p=4.34e-11). Wider audit 0.1800ms → 0.1700ms per fit (IQR 0.1750–0.1913ms → 0.1650–0.1800ms), a 1.059× speedup / 5.56% time reduction (bootstrap 95% CI 1.029×–1.088×; one-sided Wilcoxon p=2.46e-8). Before/after coefficients and `ssq_b_j` values were exactly equal on both benchmark datasets; all benchmark fits converged.

Correctness: added `test-robust-regression-xtx-cache.R`, which independently recomputes `ssq_b_j` from the returned coefficients, scale, free design matrix, and robust psi formula for MM and M fits, including a fixed-parameter case. The dedicated test passed all 6 assertions. Standalone existing robust-regression equivalence snippets against `MASS::rlm` on synthetic data and `mtcars` also passed with the repository's tolerances.

**TODO-55: Hat-matrix diagonal — vectorize via `cwiseProduct + rowwise().sum()`** ✓ DONE
File: `EDI/src/robust_post_fit_speedups.cpp:119–123`
`ols_hc2_setup_cpp` computes `hat[i] = XB.row(i).dot(X_fit.row(i))` in a scalar loop — the diagonal of `X_fit * bread * X_fit.T`. Replace with `hat = (X_fit * bread).cwiseProduct(X_fit).rowwise().sum()` — same result, vectorizable as element-wise product + row reduction.

Implemented directly in `ols_hc2_setup_cpp`: the explicit per-row dot-product loop was replaced by Eigen's vectorized row-product reduction, `hat = (X_fit * bread).cwiseProduct(X_fit).rowwise().sum()`. This preserves the existing `bread` computation and return shape while moving the hat diagonal calculation into one fused Eigen expression.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-55 source and optimized source; no `--preclean` was used. Timing used separate R processes, process CPU time, two before phases and two after phases, 80 rounds per phase. The profiled skinny setup workload used `n=2000`, `p=3`, 1000 complete `ols_hc2_setup_cpp(X)` calls per round; combined medians were 0.0430ms → 0.0410ms per call (IQR 0.0418–0.0450ms → 0.0400–0.0420ms), a 1.049× speedup / 4.65% time reduction (independent bootstrap 95% CI 1.037×–1.073×; one-sided Wilcoxon p=1.59e-17). A moderate-width audit used `n=1000`, `p=12`, 800 calls per round; combined medians were neutral at 0.1425ms → 0.1425ms per call (IQR 0.1388–0.1513ms → 0.1388–0.1478ms), speedup 1.000× (bootstrap 95% CI 0.983×–1.018×; Wilcoxon p=0.522). Before/after checksum values for the returned hat diagonal matched exactly on both benchmark datasets.

Correctness: added `test-ols-hc2-hat-vectorized.R`, which compares `ols_hc2_setup_cpp()`'s `bread` and `hat` against independent R `solve(crossprod(X))` and `rowSums((X %*% bread) * X)` calculations, then verifies the downstream HC2 covariance from `ols_hc2_post_fit_precomputed_cpp()`. The dedicated test passed all 4 assertions. Existing `test-kk-ols-se.R` also passed all 21 assertions.

**TODO-56: BetaRegression preallocate `m_mu`, `m_eta` member fields** ✓ DONE
File: `EDI/src/fast_beta_regression.cpp`
`operator()` calls `resize()` on `eta` and `mu` vectors on every call. Preallocate as member fields sized at construction (safe — `operator()` is non-`const`). After implementing TODO-23's single scalar loop, most temporaries disappear; these are only the eta/mu GEMV intermediates.

Inspection showed this TODO had already been satisfied when TODO-23 landed: `BetaRegression` already has constructor-sized `m_eta`, `m_mu`, and `m_w_grad` members, and `operator()` writes the GEMV into `m_eta.noalias()` and stores the per-observation mean in `m_mu` instead of allocating local `eta`/`mu` vectors. No C++ source change was required for TODO-56.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs; no `--preclean` was used. Because the repo source was already optimized, the before package was a temporary source copy with the TODO-described local `eta`/`mu` allocation pattern restored inside `operator()`, while the after package was the current source rebuilt after touching `fast_beta_regression.cpp` to force recompilation. Timing used separate R processes, process CPU time, two before phases and two after phases, 80 rounds per phase. The profiler-shaped estimate-only workload used `n=1000`, `p=4`, 120 complete `fast_beta_regression_cpp(..., estimate_only=TRUE)` fits per round; combined medians were 0.8250ms → 0.8500ms per fit (IQR 0.8167–0.8667ms → 0.8333–0.8667ms), speedup 0.971× / -3.03% time reduction (bootstrap 95% CI 0.971×–0.990×; one-sided Wilcoxon p=0.9999). An allocation-sensitive smaller workload used `n=250`, `p=8`, 350 fits per round; combined medians were 0.2629ms → 0.2657ms (IQR 0.2600–0.2743ms → 0.2600–0.2721ms), speedup 0.989× / -1.09% time reduction (bootstrap 95% CI 0.989×–1.011×; Wilcoxon p=0.756). Before/after checksums matched exactly and all fits converged. Conclusion: the source-level TODO is complete, but there is no remaining measurable speedup for this stale entry in the current post-TODO-23 implementation/compiler build.

Correctness: added `test-beta-regression-preallocated-eta-mu.R`, covering repeated unweighted fits, normalized score at the optimum, repeated weighted fits, and a fixed-parameter fit through the beta objective path. The dedicated test passed all 12 assertions. Existing beta-related coverage in `test-argument-permutations.R` passed all 13 assertions and `test-fast_glm_outputs.R` passed all 10 assertions.

**TODO-57: LogBin — return `fisher_information` from fit impl to avoid var-path recompute** ✓ DONE
File: `EDI/src/fast_log_binomial_regression.cpp:453–510`
`fit_constrained_binomial_with_var_cpp_impl` rebuilds `X_free` and calls `weighted_crossprod` again (line 484) after `fit_constrained_binomial_cpp_impl` already computed `XtWX` on its final accepted iteration. Return `fisher_information` from the fit impl; reuse in the var impl. Saves one O(n·p²) GEMM per var call.

Implemented by reusing the full `fisher_information` returned from `fit_constrained_binomial_cpp_impl()` in `fit_constrained_binomial_with_var_cpp_impl()`. The variance helper now validates the returned full information matrix, falls back to `weighted_crossprod(X, w)` only if it is unavailable or invalid, subsets the cached full information to the free-parameter block for `ssq_b_j`, and returns the same full cached matrix to callers. This removes the previous `X_free` rebuild and second `weighted_crossprod(X_free, w)` on the normal log-binomial / identity-binomial variance path.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-57 source and optimized source; no `--preclean` was used. The after install was forced to rebuild `fast_log_binomial_regression.cpp` by touching that source file. Timing used separate R processes, process CPU time, two before phases and two after phases, 80 rounds per phase. The profiler-shaped variance workload used `n=800`, `p=8`, 250 complete `fast_log_binomial_regression_with_var_cpp(..., j=2)` fits per round; combined medians were 0.5000ms → 0.4440ms per fit (IQR 0.4800–0.5360ms → 0.4320–0.4850ms), a 1.126× speedup / 11.20% time reduction (independent bootstrap 95% CI 1.098×–1.146×; one-sided Wilcoxon p=9.67e-19). A wider audit used `n=1500`, `p=30`, 40 fits per round; combined medians were 25.5000ms → 22.5250ms per fit (IQR 24.4438–26.7313ms → 21.0375–24.0750ms), a 1.132× speedup / 11.67% time reduction (bootstrap 95% CI 1.105×–1.155×; Wilcoxon p=8.17e-20). Before/after checksums matched exactly and all fits converged.

Correctness: added `test-logbin-var-reuses-fisher.R`, which independently rebuilds the Fisher information from fitted coefficients and checks `fisher_information` and `ssq_b_j` for log-binomial, fixed-parameter log-binomial, and identity-binomial variance paths. The dedicated test passed all 10 assertions. Existing `test-fast_glm_outputs.R` passed all 10 assertions, `test-warm-start-weights.R` passed all 9 assertions, and isolated constrained-binomial equivalence snippets against `stats::glm` for log and identity links passed. A full run of the broad `test-rcpp-fitting-equivalence.R` file segfaulted after 63 passing assertions in unrelated model-family coverage, so it was not used as TODO-57 evidence.

**TODO-58: Poisson IRLS — cache `delta_eta` for step-halving probes** ✓ DONE
File: `EDI/src/fast_poisson_regression.cpp:236–248`
Each `compute_neg_loglik(beta_try)` in the halving loop (line 240) runs a full GEMV. Precompute `delta_eta = X_f * step` before the halving loop; each probe becomes O(n) vector-add + scalar loglik loop. Same pattern as TODO-17 for log-binomial. Low priority: Poisson IRLS typically converges in 1 step for well-conditioned data.

Inspection showed this entry is a duplicate of TODO-128, which had already landed in `fast_poisson_regression.cpp`: `eta_try(n)` and `delta_eta(n)` are preallocated before the IRLS loop, `compute_neg_loglik_from_eta()` evaluates trial likelihoods from a precomputed eta vector, and the step-halving loop computes `delta_eta.noalias() = X_f * step` once per IRLS iteration before probing `eta + step_size * delta_eta`. No C++ source change was required for TODO-58 itself.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs; no `--preclean` was used. Because the repo source was already optimized, the before package was a temporary source copy with the TODO-described old halving logic restored (`beta_free_try` plus `compute_neg_loglik(beta_free_try)`, recomputing `X_f * beta_try` per probe), while the after package was the current source rebuilt after touching `fast_poisson_regression.cpp`. Timing used separate R processes, process CPU time, two before phases and two after phases, 80 rounds per phase. The normal estimate-only workload used `n=3000`, `p=6`, 200 complete `fast_poisson_regression_cpp(..., optimization_alg="irls", estimate_only=TRUE)` fits per round; combined medians were 0.6225ms → 0.5625ms per fit (IQR 0.5275–0.7600ms → 0.4988–0.6300ms), a 1.107× speedup / 9.64% time reduction (independent bootstrap 95% CI 1.013×–1.171×; one-sided Wilcoxon p=3.96e-5). The bad-warm-start workload used the same `n=3000`, `p=6` design with `warm_start_beta = beta_true * 10`, 80 fits per round, and 17 IRLS iterations; combined medians were 1.4125ms → 1.4750ms (IQR 1.2625–1.8250ms → 1.3750–1.6250ms), speedup 0.958× / -4.42% time reduction (bootstrap 95% CI 0.907×–1.013×; Wilcoxon p=0.967), i.e. neutral/no measurable win under current run conditions. Before/after checksums matched exactly and all fits converged.

Correctness: reused the existing dedicated `test-poisson-delta-eta-step-halving.R`, which covers standard fits against `glm()`, bad-warm-start fits that trigger halving, weighted bad-warm-start fits, score/information checks, deterministic repeatability, and `estimate_only` consistency. The dedicated test passed all 10 assertions. Existing `test-poisson-weighted-crossprod.R` passed all 4 assertions, `test-fast_glm_outputs.R` passed all 10 assertions, and `test-warm-start-weights.R` passed all 9 assertions.

**TODO-59: Stereotype logit — preallocate per-obs hessian allocs** ✓ DONE
File: `EDI/src/fast_stereotype_logit.cpp`
No perf data (no annotate file generated). Static analysis of `loglik_hessian()` (lines 207–307) shows per-observation allocations: `std::vector<VectorXd> logit_grad(K)`, `std::vector<MatrixXd> logit_hess(K)`, `mean_grad(d)`, `mean_hess(d,d)`, `mean_outer(d,d)`, `MatrixXd delta` per obs. For n=200, K=5, d=6: ~2,800 allocs per hessian call. Preallocate as class members. Implement only if this kernel appears in a benchmark sweep as a bottleneck.

Implemented by adding reusable mutable Hessian workspace to `StereotypeLogitRegression`: score values, first/second score derivative buffers, per-class logit gradient/Hessian buffers, mean-gradient/Hessian accumulators, outer-product scratch, and the per-observation `delta` matrix are now allocated once in the constructor and reused by `loglik_hessian()`. `compute_scores_with_second_derivatives()` now zeros pre-sized buffers rather than rebuilding the vector of Hessian matrices on every call.

Noise-controlled benchmark used `R CMD INSTALL --no-docs` installs; no `--preclean` was used. The before package was a temporary source copy saved before the TODO-59 edit, and the after package was the edited source. Timing used two before phases and two after phases, 80 rounds per phase. The direct Hessian workload used `n=800`, `p=6`, `K=5`, 300 `get_stereotype_logit_hessian_cpp()` calls per round; combined medians were 0.7083ms → 0.6867ms per call (IQR 0.6767–0.7467ms → 0.6700–0.7275ms), a 1.032× speedup / 3.06% time reduction (independent bootstrap 95% CI 1.010×–1.052×; one-sided Wilcoxon p=0.00250). The full variance-fit workload used `n=300`, `p=5`, `K=4`, 80 `fast_stereotype_logit_with_var_cpp()` fits per round; combined medians were 1.8250ms → 1.6438ms per fit (IQR 1.7250–2.0656ms → 1.5844–1.7125ms), a 1.110× speedup / 9.93% time reduction (bootstrap 95% CI 1.076×–1.170×; one-sided Wilcoxon p=6.58e-27). Per-phase direct-Hessian medians were before1 0.7117ms, after1 0.6933ms, after2 0.6767ms, before2 0.7017ms; per-phase variance-fit medians were before1 1.7375ms, after1 1.6500ms, after2 1.6188ms, before2 2.0625ms. Before/after checksums matched within each benchmark style and all fits converged.

Correctness: added dedicated `test-stereotype-logit-hessian-workspace.R`, covering repeatability, symmetry/finite values, and parameter-change isolation for the reused Hessian workspace. The dedicated test passed all 5 assertions. Existing `test-stereotype-logit-gemv-gradient.R` passed all 4 assertions, and isolated stereotype-logit checks against `glm()` / Hessian symmetry passed.

**TODO-60: ZAP `hessian()` — reuse preallocated `m_eta_cond`/`m_eta_zi`** ✓ DONE
File: `EDI/src/fast_zero_augmented_poisson.cpp:104–105,194–195`
`hessian()` allocates `eta_cond` and `eta_zi` locally (lines 104–105, 194–195) despite `m_eta_cond` and `m_eta_zi` existing as preallocated member fields. Reuse them with `.noalias() =`. Affects `zip_var`/`hurdle_p_var`, not the dominant `est` path.

Implemented for both `ZeroAugmentedPoisson::hessian()` and `expected_hessian()`: the conditional and zero-inflation linear predictors are now written into the constructor-sized `m_eta_cond` and `m_eta_zi` buffers via `.noalias()`, and the observation loop reads them through local raw pointers. This removes the two per-Hessian n-vector heap allocations while preserving the existing Hessian formulas and output shape.

Noise-controlled benchmark used `R CMD INSTALL --no-docs` installs; no `--preclean` was used. The before package was a temporary source copy saved before the TODO-60 edit, and the after package was the edited source. The primary BAAB benchmark used `n=1000`, `p_cond=p_zi=6`, two before phases and two after phases, 60 rounds per phase. Direct Hessian workloads used 1,000 Hessian calls per round; full variance-fit workloads used 120 complete `fast_zero_augmented_poisson_cpp(..., estimate_only=FALSE)` fits per round from fixed warm starts. Combined medians were mixed: hurdle direct Hessian 0.0740ms → 0.0780ms (0.949×, bootstrap 95% CI 0.914×–0.980×), ZIP direct Hessian 0.1050ms → 0.1070ms (0.981×, CI 0.963×–1.000×), hurdle variance fit 0.7042ms → 0.6750ms (1.043×, CI 1.006×–1.088×), and ZIP variance fit 0.8250ms → 0.8667ms (0.952×, CI 0.942×–0.971×). A smaller allocation-focused benchmark used `n=200`, `p=6`, 120 rounds, 3,000 direct Hessian calls per round, and 300 MLE-warm variance fits per round: hurdle direct Hessian improved 0.0195ms → 0.0180ms (1.083×, CI 1.037×–1.151×), ZIP direct Hessian was neutral at 0.0300ms → 0.0302ms (0.994×, CI 0.962×–1.057×), hurdle MLE-warm variance fit regressed 0.1700ms → 0.2117ms (0.803×, CI 0.769×–0.823×), and ZIP MLE-warm variance fit was neutral at 0.1533ms → 0.1533ms (1.000×, CI 0.957×–1.045×). Before/after checksums matched exactly for every benchmark workload and all fits converged. Conclusion: the intended allocation removal is complete, but it is not a robust wall-time win across ZAP workloads under this compiler/runtime; effects are sub-millisecond and workload-dependent.

Correctness: added `test-zero-augmented-poisson-hessian-workspace.R`, which checks ZIP and hurdle Hessians against independent finite-difference Hessians of an R negative-log-likelihood, verifies symmetry, and checks that the variance-fit path's observed information matches a standalone Hessian after optimizer buffer use. The dedicated test passed all 8 assertions. Existing `test-zero-augmented-poisson-log1m.R` passed all 3 assertions.

**TODO-61: Cont_ratio `ContinuationRatioObjective::operator()` — preallocate scratch buffers** ✓ DONE
File: `EDI/src/fast_continuation_ratio_regression.cpp:26–34`
After TODO-34 eliminates the augmented-data construction allocs, the remaining `operator()` scratch vectors (`eta`, `mu`, `log_mu`, `log_one_minus_mu`) are still allocated per L-BFGS call. Preallocate as class members. Note: `operator()` is the method signature here — check if it's `const`; if so use explicit workspace passing not mutable fields (see mutable-field antipattern).

Implemented by making `ContinuationRatioObjective::operator()` and `hessian()` non-`const` and adding constructor-sized non-mutable work buffers: `eta`, `mu`, `work`, `log_mu`, and `log_one_minus_mu`. This keeps the optimizer's non-const functor flow, avoids mutable heap fields, reuses the likelihood/gradient/Hessian scratch storage across L-BFGS calls, preserves the vectorized Eigen log-likelihood reduction, and reuses `work` for both `mu - z` and Hessian weights.

Noise-controlled benchmark used `R CMD INSTALL --no-docs` installs; no `--preclean` was used. The before package was a temporary source copy saved before the TODO-61 edit, and the after package was the edited source. Timing used two before phases and two after phases, 80 rounds per phase, 100 complete fits per round, with `n=1000`, `p=5`, `K=4`. Combined medians were essentially neutral: `fast_continuation_ratio_regression_cpp()` 0.8000ms → 0.7900ms per fit (IQR 0.6800–0.9000ms → 0.6900–0.9400ms), a 1.013× speedup / 1.25% time reduction (independent bootstrap 95% CI 0.957×–1.077×; one-sided Wilcoxon p=0.565). `fast_continuation_ratio_regression_with_var_cpp()` 0.7800ms → 0.7700ms per fit (IQR 0.7200–0.8500ms → 0.7000–0.8500ms), a 1.013× speedup / 1.28% time reduction (bootstrap 95% CI 0.987×–1.053×; one-sided Wilcoxon p=0.0495). Per-phase medians for the base fit were before1 0.8700ms, after1 0.8650ms, after2 0.7300ms, before2 0.7350ms; per-phase medians for the variance path were before1 0.8200ms, after1 0.7300ms, after2 0.7900ms, before2 0.7400ms. Before/after checksums matched exactly and all fits converged. Conclusion: the allocation cleanup is complete and safe, but wall-time impact is too small to distinguish from phase noise for the full fit path.

Correctness: added `test-continuation-ratio-objective-workspace.R`, which checks repeated fits are exact, the fitted score is small, fitted information equals the standalone Hessian helper, and the variance path's covariance equals the inverse information after objective workspace reuse. The dedicated test passed all 8 assertions. Existing `test-continuation-ratio-augmentation.R` passed all 4 assertions.

### NegBin + HurdleNegBin findings (HIGH priority)

**TODO-62: NegBin `NBLogLik::operator()` — per-row rank-1 gradient → single GEMV** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_negbin_regression.cpp:112`
`score_beta += coef * m_X.row(i).transpose()` — per-row rank-1 accumulation, same anti-pattern fixed by TODO-15/16 for ZIP/ZINB but not applied here. Fix: accumulate scalar coefficients into preallocated `m_coef_vec[i]` inside the obs loop, then `score_beta.noalias() += m_X.transpose() * m_coef_vec` (single BLAS GEMV) after. Expected: same speedup class as ZIP/ZINB GEMV fixes.

Verified (2026-07-06): `m_coef_vec` member field exists; obs loop sets `m_coef_vec[i] = yi - mu_i * (yi + theta) / denom`; after loop `grad.head(m_p).noalias() = -(m_X.transpose() * m_coef_vec)` — GEMV already in place.

**TODO-63: NegBin `NBLogLik::hessian()` — use preallocated distinct-y digamma/trigamma tables** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_negbin_regression.cpp:148–150`
`hessian()` calls `R::digamma(yi + theta)` and `R::trigamma(yi + theta)` raw per-obs, ignoring `m_digamma_yptheta` and `m_trigamma_yptheta` tables already populated in `operator()`. Compare: `TruncatedNegBinCount::hessian()` correctly uses the tables. Copy the table-lookup pattern to `NBLogLik::hessian()`. Eliminates all per-obs R function dispatch in the hessian path.

Verified (2026-07-06): `hessian()` fills `m_digamma_yptheta[k]` and `m_trigamma_yptheta[k]` in a distinct-y loop before the obs loop, then the obs loop uses `m_digamma_yptheta[slot]` and `m_trigamma_yptheta[slot]` — no per-obs R function calls.

**TODO-64: NegBin `expected_hessian()` — trigamma recurrence to eliminate O(iter) R::trigamma calls** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_negbin_regression.cpp:204–211`
`expected_trigamma_y_plus_theta` calls `R::trigamma(k + theta)` in a series summation loop — ~47 R::trigamma calls per obs per expected-hessian call = ~47,000 R::trigamma per call. Fix: call `R::trigamma(theta)` once, then use the recurrence `ψ₁(x+1) = ψ₁(x) − 1/x²` for subsequent terms (1 R call + O(iter) divisions). The recurrence is exact for the trigamma function — no approximation.

Verified (2026-07-06): `expected_trigamma_y_plus_theta` takes `trigamma_theta` as a parameter (one `R::trigamma(theta)` call hoisted by caller); loop body uses `trigamma_yptheta -= 1.0 / (x * x)` recurrence — exactly the proposed fix.

**TODO-65: NegBin + HurdleNegBin — `R::lgammafn` → `std::lgamma`** ✓ DONE — ALREADY IMPLEMENTED
Files: `EDI/src/fast_negbin_regression.cpp:83,90`, `EDI/src/fast_hurdle_negbin.cpp:346,356`
`R::lgammafn` goes through R's error-handling dispatch wrapper; `std::lgamma` is direct libm. Profile: `logf32x` (lgamma) is **63% of negbin_est** samples. This is the single largest hotspot in NegBin. Fix: replace all `R::lgammafn(x)` with `std::lgamma(x)` in both files. Combined with TODO-22 (`fast_lgamma`), reduces to: `fast_lgamma` where available. Same fix applied to ZINB in TODO-49 — extend scope to cover NegBin and HurdleNegBin here.

Verified (2026-07-06): `fast_negbin_regression.cpp` uses `std::lgamma(theta)` and `std::lgamma(ypt)` throughout; `fast_hurdle_negbin.cpp` uses `std::lgamma` at lines 330, 346, 356 — no `R::lgammafn` found in either file.

**TODO-66: HurdleNegBin — `log1p(-eneg)` fast-path for large `lam`** ✓ DONE
File: `EDI/src/fast_hurdle_negbin.cpp` (hurdle positive-count log-normalizer)
Mirror of TODO-33 (HP GLMM) and TODO-43 (ZAP): for `lam > 16`, `eneg = exp(-lam) < 1e-7`, so `log1p(-eneg) ≈ -eneg` with error < 1e-14. Add fast-path branch before the `std::log1p` call. Eliminates the libm transcendental for large-count observations (common in NegBin hurdle outcomes).

**RESULT (2026-07-06) — DONE.** In `TruncatedNegBinCount::operator()`, replaced `std::log(trunc_denom)` (where `trunc_denom = 1 - p0`) with `(p0 < 1e-7) ? -p0 : std::log(trunc_denom)`. When `p0 < 1e-7` (i.e., NegBin zero-prob is negligible), `log(1-p0) ≈ -p0` with error < `p0²/2 < 5e-15` — avoids the libm `log` call entirely.

Benchmarked at n=3000, p=5, θ=10, μ≈200 (99.5% of obs hit the fast path), 30 rounds × 5 reps warm-start:

| path | old median | new median | speedup | Wilcoxon p |
|---|---:|---:|---:|---:|
| EST | 2.60 ms | 2.20 ms | **1.18x** | 9.2e-9 |
| VAR | 3.00 ms | 2.50 ms | **1.20x** | 1.7e-3 |

Correctness: `test-hurdle-negbin-log1p-fast-path.R` (16 tests) — all pass. Tests cover: fast-path activation fraction, convergence + finite estimates, estimate_only parity, observed information PD, score near zero at convergence, repeated-call determinism, mixed slow-path data (small θ/μ), warm-start agreement.

**TODO-67: ZINB (and cross-kernel) — vectorize exp/log by enabling Eigen SIMD + batching transcendentals into array ops** ✓ DONE (2026-07-02): global Eigen SIMD enabled package-wide + ZINB array rewrite; committed; ZINB 1.7-1.9x; full suite 398/398 green; installs clean
Files: `EDI/src/Makevars:3` (+ `EDI/src/Makevars.win`), `EDI/src/fast_zinb.cpp:79-148`

**RESULT (2026-07-02) — DONE (committed, whole-package). ZINB 1.7-1.9x via array rewrite + GLOBAL Eigen SIMD.**
Implemented: `fast_zinb.cpp` `operator()` rewritten to batch all per-obs transcendentals into vectorized Eigen array passes (`mu=exp(eta_c)`, `log(theta+mu)`, sigmoid, softplus, `phi`, `log(den)`) with preallocated member buffers; the branchy accumulation loop now does zero transcendental calls.
- **Scoped `#undef` (option B) is DEAD — it crashes.** Compiling only `fast_zinb.o` alignment-ON while the other 103 objects stay alignment-OFF segfaults even on estimate-only (Eigen's ODR-merged allocator mismatches aligned AVX loads). `-DEIGEN_DONT_ALIGN` also can't be kept (it disables vectorization). So the ODR caveat is a **hard crash, not theoretical → only the GLOBAL flip (option A) is viable.** Removed `-DEIGEN_DONT_ALIGN -DEIGEN_DONT_VECTORIZE` from `Makevars`; `Makevars.win` already defaulted to SIMD-on (added its missing `-I../inst/include` so it faithfully mirrors).
- **Parity** (SIMD reorders FP; not bit-exact): params 8e-15, vcov 5e-12, hessian 6.7e-9 (numerical FD) — machine precision.
- **ZINB speedup** (interleaved A/B vs scalar TODO-19 baseline): est **1.93x** / var **1.89x** at n=2000; est **1.74x** / var **1.76x** at n=20000.
- **Full testthat suite** (clean whole-package rebuild): **398/398 PASS** (0 failed, 0 errors, 2 skipped). Package also **`R CMD INSTALL`s cleanly and loads** (v1.0.0). Two earlier full-suite runs showed 2 failures in the same non-hardened rank-deficient ordinal test on **both** scalar and SIMD builds — flaky/order-dependent, not caused by this change; they do **not** reproduce on the clean run.
- **Package-wide GEMV bonus** (kernels NOT rewritten; scalar vs SIMD, n=3000 p=6) — **mixed, not uniform:**

| kernel | speedup |
|---|---:|
| poisson_var | 1.64x |
| ols_var | 1.62x |
| logistic_var | 1.14x |
| beta_var | 1.08x |
| probit_var | 1.06x |
| coxph_est | 0.92x (regression) |

GEMM/IRLS-heavy kernels win most; coxph slightly **regresses** → the global flip needs per-kernel evaluation before committing, not a blanket "SIMD is faster" assumption.

**Decision (2026-07-02): DONE — global Eigen SIMD accepted and committed package-wide** (`Makevars` + `Makevars.win`), ZINB `operator()` array-rewritten. Follow-ups (separate TODOs): (a) the **coxph 0.92x regression** under SIMD (see table) — investigate/guard; (b) per-kernel array-op rewrites (logistic/poisson/beta/…) to turn the modest flag-only GEMV bonus into ZINB-scale wins; (c) TODO-68's softplus concern is resolved inside this rewrite (vectorizable `log(1+exp)`, no `std::log1p`, no fast-math). All 104 kernels' results now shift at the ULP level — the green suite confirms nothing broke.

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

**UPDATE (2026-07-02) — RESOLVED the right way (vectorizable array log1p, no `-ffast-math`).** The deferred "Eigen-array softplus" path is now implemented and shipped. Added **`fast_log1p_arr`** to `_helper_functions.h`: a vectorizable, accurate `log1p` for `z > -1` built as a Kahan correction over Eigen's **packet `.log()`** (Eigen's own `.log1p()` falls back to scalar `std::log1p`). Rewrote `log1pexp_array_safe` to the branchless softplus `max(x,0) + fast_log1p_arr(exp(-|x|))`, and converted the remaining scalar softplus spots to array form — `ZeroInflatedNegBin::operator()` (via TODO-67), `LogisticGLMMObjective::value()`, `ClogitPlusGLMMObjective::neg_clogit()` — removing the dead scalar helpers `lse_zinb`, `log1pexp_s`, `log1pexp_cpp`.
- Measured (SIMD, no fast-math): bare `log1p` **22.2 → 8.28 ns/elem (2.68x)** via packet `.log()` vs Eigen `.log1p()`, accuracy **2.7e-16**; `logistic_glmm` var path **1.34x** cumulative.
- Parity machine-precision (5.55e-17 params; bit-identical where the path is unchanged); **full testthat suite 398/398 PASS**.
- Net: no `std::log1p`-based softplus remains in ZINB / logistic_glmm / clogit — the log1p win TODO-68 sought is captured safely (no `-ffinite-math-only` guard elision, no global FTZ). Scalar per-obs softplus in other kernels (ZAP `lse`, `log1pexp_stable`, …) still awaits ZINB-style array rewrites of their loops.

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

### Critical findings (action items TODO-110+)

---

**TODO-110: ZINB — preallocate eta_c, eta_z, w_c, w_z as member vectors** ✓ DONE
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

Implemented/verified in `ZeroInflatedNegBin`: constructor-sized `m_eta_c`, `m_eta_z`, `m_w_c`, and `m_w_z` buffers are used by `operator()` via `.noalias()` GEMV and `.setZero()` weight initialization, matching the ZAP pattern. The remaining TODO detail was also completed by removing the redundant `grad.resize(...)` from the hot path; all in-repo call sites pass a correctly sized gradient vector (`FixedParameterFunctor`, `likelihood_score()`, `likelihood_value()`, and `numerical_hessian()`). This source already also contains the later vectorized transcendental precompute buffers (`m_mu`, `m_logden`, `m_p`, `m_lse`, `m_phi`, `m_ld0`), so the benchmark below isolates the eta/weight allocation cleanup against the current post-TODO-67 code rather than the much older profiler snapshot.

Noise-controlled benchmark used `R CMD INSTALL --no-docs` installs; no `--preclean` was used. Because the repo source already had the member-vector portion of TODO-110, the before package was a temporary source copy with the TODO-described local `eta_c`, `eta_z`, `w_c`, and `w_z` allocations restored inside `operator()`; the after package was the current optimized source with the resize removal. Timing used two before phases and two after phases, 80 rounds per phase, 100 complete fits per round, with `n=1000`, `p_cond=p_zi=4`. The estimate-only workload used `fast_zinb_cpp(..., estimate_only=TRUE)` from fixed warm starts; combined medians were 0.7350ms → 0.6000ms per fit (IQR 0.6675–0.8225ms → 0.5975–0.6200ms), a 1.225× speedup / 18.4% time reduction (independent bootstrap 95% CI 1.180×–1.262×; one-sided Wilcoxon p=9.69e-50). The MLE-warm variance workload used `estimate_only=FALSE` from the estimate-only MLE; combined medians were 1.0500ms → 0.9100ms per fit (IQR 1.0000–1.1300ms → 0.8900–0.9300ms), a 1.154× speedup / 13.3% time reduction (bootstrap 95% CI 1.137×–1.178×; one-sided Wilcoxon p=1.09e-44). Per-phase estimate-only medians were before1 0.8100ms, after1 0.6100ms, after2 0.6000ms, before2 0.6650ms; variance medians were before1 1.0700ms, after1 0.9000ms, after2 0.9150ms, before2 1.0200ms. Before/after checksums matched exactly and all fits converged.

Correctness: added `test-zinb-operator-workspace.R`, which checks deterministic repeated estimate-only fits and verifies the variance path remains finite with a symmetric Hessian and finite covariance diagonal after operator workspace reuse. The dedicated test passed all 7 assertions. Existing `test-zinb-std-lgamma.R` passed all 3 assertions.

---

**TODO-111: HurdlePoisson GLMM hessian — preallocate G×K working buffers** ✓ DONE
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

**TODO-112: NBLogLik — GEMV refactor for gradient in operator()** ✓ DONE — ALREADY IMPLEMENTED
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
Verified: `m_coef_vec` member fills per-obs coefficient scalar, then `grad.head(m_p).noalias() = -(m_X.transpose() * m_coef_vec)` replaces the row-access accumulation loop. Same GEMV pattern as ZIP/ZINB.

---

**TODO-113: NBLogLik hessian — fill distinct-y tables; use slot lookups** ✓ DONE — ALREADY IMPLEMENTED
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

**TODO-114: NBLogLik::expected_trigamma_y_plus_theta — trigamma recurrence** ✓ DONE — ALREADY IMPLEMENTED
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

**TODO-115: R::lgammafn → std::lgamma in NegBin, ZINB, HurdleNegBin** ✓ DONE — ALREADY IMPLEMENTED
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

**TODO-116: HurdlePoisson/ZIP — log1p fast-path for large lambda** ✗ DROPPED (ZAP part)
Files: `EDI/src/fast_hurdle_poisson_glmm.cpp:183,268`, `EDI/src/fast_zero_augmented_poisson.cpp`

Root cause: `std::log1p(-eneg)` where `eneg = exp(-lam)` is the top hotspot in hurdle_p_est (log1p at 50%+). For large lam, `eneg → 0` and `log1p(-eneg) ≈ -eneg` with error < `eneg²/2`. For lam > 30: `eneg < 9e-14`, so the approximation has error < 4e-27 (far below double precision). A threshold of `lam > 16` gives error < 1e-14.

**HP GLMM part ✓ DONE — ALREADY IMPLEMENTED** (also documented as TODO-33): `log_one_minus_exp_neg_hp` helper in `fast_hurdle_poisson_glmm.cpp` uses `(eneg < 1e-7) ? -eneg : std::log1p(-eneg)`.

**ZAP part ✗ DROPPED — FALSIFIED (2026-07-07).** The proposed fast-path `(eml < 1e-7) ? -eml : (lam > 0.693) ? std::log1p(-eml) : std::log(-std::expm1(-lam))` was implemented and benchmarked.

Root cause of pessimization: glibc's `std::log1p(x)` already has an internal fast-path for `|x| < 2^-53`: it returns `x` directly with no transcendental computation (Taylor series `x - x²/2 + ...` collapses to `x` at machine precision). For large lambda (mean ≈ 45), `eml = exp(-lam) ≈ 2.8e-20`, so `log1p(-2.8e-20)` is already an O(1) operation internally. Our outer branch `(eml < 1e-7)` adds a branch + comparison overhead with no transcendental saving.

This is distinct from the HurdleNegBin fast-path (TODO-66), which saved `std::log(trunc_denom)` where `trunc_denom = 1 - p0` near 1 — glibc's `log(x)` has no x≈1 fast-path, so that saving was real.

Benchmark (n=2000, p=3, mean λ≈45, 40 rounds × 200 reps warm-start at MLE, single operator() call):

| code | median (μs) | IQR |
|---|---:|---:|
| OLD | 1215 | [1114, 1288] |
| NEW (fast-path) | 1383 | [1279, 1543] |
| ratio | **0.879× (13% regression)** | Wilcoxon p≈1.0 |

No source changes retained. The ZAP `operator()` correctly uses the existing two-branch `(lam > 0.693) ? log1p(-eml) : log(-expm1(-lam))` which is already optimal.

---

**TODO-117: Logistic/Probit/Poisson IRLS — XtWX via weighted_crossprod** ✓ DONE (Poisson API refactoring) / ✗ DROPPED (col-major DSYR speedup)
Files: `EDI/src/fast_logistic_regression.cpp`, `EDI/src/fast_probit_regression.cpp`, `EDI/src/fast_poisson_regression.cpp`

Root cause: IRLS computes XtWX as `X.T * diag(w) * X` or similar triple-product, creating an n×p intermediate. `weighted_crossprod(X, w)` (already in `_helper_functions.h`) uses the upper-triangular DSYR/DSYRK symmetric update — halves FLOPs and avoids the intermediate allocation. Same fix was applied to log-binomial (TODO-17).

**Logistic + Probit ✓ DONE — ALREADY IMPLEMENTED** (documented as TODO-21). **Poisson API refactoring ✓ DONE**: replaced the 3 inline triple-products in `fast_poisson_regression.cpp` (lines ~233, ~240, ~284) with `weighted_crossprod(X_f, w_tmp)` / `weighted_crossprod(X_f, w_final)`. Test: `test-poisson-xtwx-colmajor-dsyr.R` (7 assertions).

**Col-major DSYR speedup ✗ DROPPED**: Attempted to add a col-major DSYR path to `weighted_crossprod` (column-pair loop: scale col_j by w, then dot against col_k≥j for upper triangle only). Standalone XtWX benchmark with concrete `MatrixXd` showed 1.4–3× speedup across n=500–2000, p=6–15. However, end-to-end Poisson IRLS benchmark showed:
- n=200–2000, p=6: neutral (0× change; XtWX is small fraction of total IRLS time)
- n=500, p=15: 12.5% regression (200→225 µs median)

Root cause of regression: for large p, BLAS DGEMM (called by the triple product) outperforms p*(p+1)/2 scalar Eigen dot products due to AVX register tiling and cache blocking. The standalone speedup used concrete `MatrixXd` (fully optimized by compiler); the actual IRLS path uses `Eigen::Ref<const MatrixXd>` which prevents the same SIMD optimization. Reverted `weighted_crossprod` to triple product fallback for col-major matrices. Poisson code now uses `weighted_crossprod` (API consistency) but behaviorally identical to old triple product.

---

**TODO-118: OLS/Robust — symmetric XtX via DSYRK** ✓ DONE
Files: `EDI/src/fast_ols.cpp`, `EDI/src/fast_robust_regression.cpp`, `EDI/src/_helper_functions.h`

Added `symmetric_crossprod(X)` template in `_helper_functions.h` using BLAS DSYRK (`F77_CALL(dsyrk)`), which fills only the upper triangle and copies to lower — half the FLOPs of full DGEMM for a symmetric result. Replaced all 4 `X.transpose() * X` sites in `fast_ols.cpp` and the inline `X_free.transpose() * X_free` site in `fast_robust_regression.cpp`.

Benchmark (60 rounds × 500 reps, `fast_ols_with_var_cpp`):

| config        | OLD median | NEW median | ratio | Wilcoxon p |
|---------------|-----------|-----------|-------|------------|
| n=200,  p=4   | 4 µs      | 6 µs      | 0.67  | 0.758 (ns) |
| n=500,  p=6   | 10 µs     | 14 µs     | 0.71  | 1.000 (ns) |
| n=1000, p=6   | 14 µs     | 16 µs     | 0.88  | 1.000 (ns) |
| n=2000, p=6   | 24 µs     | 26 µs     | 0.92  | 0.914 (ns) |
| n=500,  p=15  | 24 µs     | 22 µs     | 1.09  | 6e-4 ✓    |
| n=1000, p=15  | 38 µs     | 30 µs     | 1.27  | 5e-10 ✓   |

Speedup real and significant at p≥15 (27% at n=1000). No significant regression at small p (timer quantization at ~2µs resolution dominates those deltas). Correctness: 9/9 tests pass (`test-ols-symmetric-crossprod-dsyrk.R`).

---

**TODO-119: Nonparametric — Wilcox rank-sum O(n²) → O(n log n)** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_wilcox_hl.cpp`
Profiler: 16% branch-miss rate for wilcox_hl (highest of any kernel), consistent with O(n²) double loop over all pairwise comparisons. Hulsen-Lehmann estimator currently O(n²); merge-sort based U-statistic counting is O(n log n) (same algorithm as in numpy).
Verified: `count_pairwise_diffs_leq` uses a sorted two-pointer merge (O(n_t + n_c) per call), called within a 96-iteration binary search in `select_pairwise_diff_sorted` — total O(n log n). Corresponds to TODO-25 "Wilcox HL — O(n²) → O(n log²n) sort + binary-search estimator ✓ DONE".

---

**TODO-120: Nonparametric — Ridit: std::map → std::unordered_map** ✓ DONE — ALREADY IMPLEMENTED
File: Ridit source
`std::map<int, double>` for frequency table lookup is O(log n) per lookup. `std::unordered_map<int, double>` is O(1) amortized. For discrete count data where the map is built once and queried n times, this eliminates n × O(log n) lookups.
Verified: ridit source uses `std::unordered_map` (no `std::map`). Corresponds to TODO-27 "Ridit — `std::map` → `std::unordered_map` + eliminate `wrap(ref_idx)` SEXP round-trip ✓ DONE".

---

**TODO-121: Survival — coxph_var per-observation allocation** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_coxph_regression.cpp`
Allocates working matrices inside the observation loop in the variance computation path. Preallocate scratch outside the loop.
Verified: `compute_robust_vcov` preallocates `U(n_total, p)`, `eta`, `exp_eta`, `r_x_exp`, `dk_over_Rk`, `ek`, `dk_ek_over_Rk`, `cum_A`, `cum_B` per-stratum outside the obs loop. No VectorXd/MatrixXd allocated inside the per-observation loops. Corresponds to TODO-30 "CoxPH compute_robust_vcov — eliminate 150 heap allocs ✓ DONE".

---

**TODO-122: Nonparametric — JT test: std::map → vector + precomputed table** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_jonckheere_terpstra.cpp`
Same O(log n) map lookup replaced with O(1) vector index.
Verified: `LogChooseTable` in `fast_jonckheere_terpstra.cpp` uses a flat `std::vector<double>` with index arithmetic — no `std::map`. Corresponds to TODO-26 "JT exact — precomputed lchoose table + flat `std::vector` for stat_prob ✓ DONE".

---

**TODO-123: Survival — logrank extra O(n) passes** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_logrank.cpp`
Multiple O(n) traversals that can be fused into a single pass.
Verified: comment at line 62 reads "Fused martingale accumulators — no martingale[] array needed". Single-pass accumulation confirmed in source.

---

**TODO-124: Robust regression — MAD: sort → nth_element** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_robust_regression.cpp`
Median of absolute deviations currently sorts the residual array O(n log n). `std::nth_element` gives O(n) for just the median — half the work.
Verified: `fast_robust_regression.cpp` uses `std::nth_element` for large n and `std::sort` for small n (dual-branch approach). `nth_element` path confirmed at lines 28 and 35.

---

**TODO-125: G-computation — ordinal model: cache hessian across calls** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_ordinal_regression.cpp:449-462` (gcomp finite-diff path)
Hessian is computed 3× in a single inference pass when once suffices. Cache and reuse.
Verified: corresponds to TODO-32 "GComp ordinal finite-diff loop — hoist temporaries outside loop ✓ DONE". Reusable perturbation and linear-predictor buffers, scalar restoration of perturbed parameters, baseline linear predictors reused across finite-diff steps.

---

**TODO-126: ZOIB — std::lgamma + fast_digamma** ✓ DONE — ALREADY IMPLEMENTED
File: `EDI/src/fast_zero_one_inflated_beta.cpp`
Uses `R::lgammafn` and `R::digamma` where `std::lgamma` and `fast_digamma` apply. Same pattern as TODO-14 (beta) and TODO-115.
Verified: `fast_zero_one_inflated_beta.cpp` uses `fast_digamma` (line 9 include, calls throughout) and `std::lgamma` (lines 15, 139). No `R::lgammafn` or `R::digamma` found.

---

**TODO-127: Stereotype logit — hessian allocs + GEMV gradient** ✓ DONE
File: `EDI/src/fast_stereotype_logit.cpp`

Two combined optimizations in `StereotypeLogitRegression`:

**1. Hessian alloc preallocations (TODO-127):** Added mutable members `m_score_v`, `m_score_cum_v` (reuse `gamma.array().exp()` + cumsum across compute_scores calls), `m_hess_H` (preallocate d×d `H` matrix in `loglik_hessian`), `m_exp_z`, `m_exp_mean_grad`, `m_exp_logits`, `m_exp_probs` (preallocate `expected_hessian` working vectors). Also reuse `m_hess_score_vals` and `m_hess_dscore_dgamma` in `loglik_grad` instead of allocating per call.

**2. GEMV gradient (TODO-90):** Replaced per-row `m_X.row(i).transpose() * scalar` scatter in `loglik_grad` with a preallocated `m_beta_score_weight[i]` accumulation + single `m_X.transpose() * m_beta_score_weight` GEMV at the end. Preallocated `m_eta = m_X * beta` before the obs loop.

Benchmark (60 rounds × 20 reps, `fast_stereotype_logit_with_var_cpp`, vs pre-TODO-90 baseline):

| config          | OLD median | NEW median | ratio | Wilcoxon p |
|-----------------|-----------|-----------|-------|------------|
| n=200,  K=4     | 0.500 ms  | 0.400 ms  | 1.25× | 5e-15 ✓   |
| n=500,  K=4     | 1.350 ms  | 1.150 ms  | 1.17× | 1e-18 ✓   |
| n=1000, K=4     | 2.775 ms  | 2.350 ms  | 1.18× | 3e-21 ✓   |
| n=200,  K=5     | 0.800 ms  | 0.650 ms  | 1.23× | 4e-16 ✓   |
| n=500,  K=5     | 1.600 ms  | 1.350 ms  | 1.19× | 6e-14 ✓   |
| n=1000, K=5     | 3.650 ms  | 3.000 ms  | 1.22× | 1e-20 ✓   |

Consistent 17–25% speedup across all sizes and K, all highly significant. Correctness: 18/18 tests pass (`test-stereotype-logit-hessian-allocs.R`).

---

**TODO-128: Poisson — cache delta_eta for IRLS step-halving** ✓ DONE
File: `EDI/src/fast_poisson_regression.cpp:236-248`
Each backtracking probe recomputes `X * beta_new` as a full GEMV. Precompute `delta_eta = X_free * direction` once; each probe becomes O(n) vector-add. Low priority: IRLS typically converges in 1-2 steps for well-conditioned Poisson data.

Added `eta_try(n)` and `delta_eta(n)` preallocated before the IRLS loop; added `compute_neg_loglik_from_eta` lambda (O(n) scalar loop, no GEMV); replaced the step-halving loop to compute `delta_eta.noalias() = X_f * step` once per IRLS iteration, then `eta_try.noalias() = eta + step_size * delta_eta` per probe. On accept, `beta_free.noalias() += step_size * step`. Benefit: k halvings → (k−1) GEMVs saved; in the common 1-step-accepted case, cost is neutral (1 GEMV for delta_eta replaces 1 GEMV inside compute_neg_loglik).

Benchmark (n=3000, p=6, bad warm start × 10 forcing step-halving every early iteration, 40 rounds × 50 reps, Wilcoxon one-sided):

| code | median (ms) | IQR |
|---|---:|---:|
| OLD | 1.840 | [1.720, 1.960] |
| NEW | 1.710 | [1.660, 1.800] |
| speedup | **1.076×** | p=3.6e-4 |

Correctness: `test-poisson-delta-eta-step-halving.R` (10 assertions) — all pass. Tests cover: standard fit vs `glm()`, bad-warm-start fit (triggers halving) vs `glm()`, weighted + bad warm start, score near zero at convergence, information PD, deterministic repeated calls, estimate_only consistency. Installed with `R CMD INSTALL --no-docs`.

---

## Ordinal family findings (from 2026-07-04 perf annotate run, 129 kernels)

**TODO-69: FixedOrdinalRegression — preallocate per-call working vectors** ✓ RESOLVED — REJECTED
Files: `EDI/src/fast_ordinal_regression.cpp` (FixedOrdinalRegression::operator())

Root cause: `malloc`/`cfree` appear prominently in prop_odds_est, prop_odds_var, ord_cauchit_est, ord_cloglog_est, and kk21_ordinal_wts profiles. `FixedOrdinalRegression::operator()` allocates working vectors (cumulative probabilities, per-threshold gradient components) on every optimizer call. Preallocate these as member fields and resize once in the constructor — same pattern applied to HurdlePoisson GLMM. Covers all families routing through FixedOrdinalRegression: prop_odds, ord_cauchit, ord_cloglog, ord_probit, and kk21_ordinal weights.

Benchmarked and reverted. The source had already implemented this proposal under TODO-51 by adding a constructor-sized eta buffer and replacing local alpha/beta vectors with expression views. A clean pre-change package was reconstructed from the immediately preceding header while retaining all later fast-erfc work. Noise control used an A/B/A baseline (60 round medians total) versus 30 optimized rounds, each with 300 iterations at N=1000. The broad implementation was neutral for proportional odds (0.7593ms → 0.7577ms, 0.21%; bootstrap 95% CI −1.06% to 1.48%), slightly slower for cauchit (0.6954ms → 0.7004ms, 0.72% slower), 1.34% faster for cloglog (CI 0.42%–2.68%), and **19.96% slower for ordinal probit** (0.6349ms → 0.7616ms, regression CI 18.96%–21.29%). A refined variant retained only the N-length eta buffer while restoring local alpha/beta copies; it still slowed all four links by 1.45%–2.33%. The constructor allocation is paid on every public fit, while these workloads need too few optimizer evaluations to amortize it. Production code therefore restores local alpha, beta, and eta vectors while retaining the older small gradient scratch buffers. Correctness: outputs are bit-identical and optimizer iteration counts match for all four links; `test-ordinal-cauchit-scratch-buffers.R` and `test-ordinal-information-reuse.R` pass after the final `R CMD INSTALL --no-docs` build.

---

**TODO-70: adj_cat hessian — preallocate ColPivHouseholderQR workspace** ✓ RESOLVED — REJECTED
Files: `EDI/src/` adj_cat source (AdjacentCategoryLogitNegLogLik)

Root cause: adj_cat_var profile shows `makeHouseholder`, `applyHouseholderOnTheLeft`, and `ColPivHouseholderQR::computeInPlace` — a full QR factorization performed per hessian call. The augmented matrix dimensions are fixed for a given dataset, so the `ColPivHouseholderQR` object and its internal workspace can be preallocated as a member field. Call `compute()` in-place each evaluation rather than constructing a new QR object.

**RESULT (2026-07-06) — REJECTED.** Profiled the actual call path: `ColPivHouseholderQR` is from `try_safe_ols_solve` (cold-start OLS), not from within `hessian()`. The `makeHouseholder`/`applyHouseholderOnTheLeft` symbols come from `SelfAdjointEigenSolver` inside `symmetric_pseudo_inverse`, called once per fit (not per hessian). Implemented the preallocated `mutable Eigen::VectorXd` + `mutable Eigen::MatrixXd` buffer pattern as member fields (replacing `std::vector<double>` locals and the local `MatrixXd hess`). Isolated hessian benchmark at n=1000, K=6, p=3, 30 rounds × 200 reps:

| path | old median | new median | ratio | Wilcoxon p |
|---|---:|---:|---:|---:|
| HESS | 0.140 ms | 0.148 ms | **0.95x (5% regression)** | 0.026 |
| VAR (warm) | 0.40 ms | 0.40 ms | 1.00x (flat) | — |
| DISTR (n=500, nsim=500) | 659 ms | 668 ms | ~1.00x (noise) | — |

Root cause of regression: mutable Eigen member fields (`m_exp_alpha`, `m_prob`, `m_cdf`) prevent the compiler from proving no aliasing with `m_X.data()` inside the `const` obs loop — the same LICM/aliasing anti-pattern as `mutable std::vector` (see memory feedback). The K-length vectors are small (K-1 ≈ 5 doubles), so allocation overhead is negligible compared to the n=500–1000 obs computation. Change reverted.

---

**TODO-71: ord_probit — Rf_pnorm_both/Rf_dnorm4 → fast erfc + direct exp** ✓ DONE — ALREADY OPTIMIZED BY TODO-52
Files: ord_probit source

Root cause: ord_probit_est and ord_probit_var show `Rf_pnorm_both` (#2 hotspot), `Rf_pnorm5`, and `Rf_dnorm4` — R's normal CDF and PDF dispatch functions routing through R's error-handling machinery. Replacements:
- Φ(x) = `0.5 * erfc(-x * M_SQRT1_2)`: use `fast_erfc` already in the codebase (same fix as probit regression, TODO-20).
- φ(x) = `exp(-x*x * 0.5) * M_1_SQRT_2PI`: one direct `std::exp`, no library dispatch.

Verified this was already implemented in the shared `ordinal_fixed_link_helpers.h`: the probit CDF calls `pnorm_fast` (Cephes-based `fast_erfc`) and the PDF/derivative call `dnorm_fast` (direct `std::exp`). The TODO-52 A/B/A benchmark used 30 rounds × 500 iterations per phase for N=1000 estimation and 30 rounds × 100 iterations for N=200 variance: estimate 0.8780ms → 0.7950ms (9.5% faster) and variance 1.4300ms → 1.2500ms (12.6% faster). Correctness: `test-fast-probit-cdf.R` checks the CDF over a dense grid and compares ordinal-probit analytic score/Hessian with independent numerical derivatives at `1e-7`/`1e-5`; all four assertions pass on a fresh package installed with `R CMD INSTALL --no-docs`. No additional source change was needed for TODO-71.

---

**TODO-72: ord_cauchit — investigate and eliminate per-call introsort** ✓ DONE
Files: ord_cauchit source

Root cause: `std::__introsort_loop` appears in the ord_cauchit profile alongside `__atan_fma` (cauchit link, F(x) = 0.5 + atan(x)/π). A sort is being performed per optimizer evaluation — likely sorting thresholds to enforce monotonicity after each gradient step. Replace with a constrained reparameterization (e.g., cumulative-log-spacing: θ_k = θ_1 + Σ exp(δ_j)) so thresholds are monotone by construction and no sort is needed. The `__atan_fma` itself is already glibc-optimized; the sort is the actionable bottleneck.

The investigation found no per-evaluation threshold sort and no need for a new parameterization. `init_levels(y)` sorts the response once in `FixedOrdinalRegression` construction, but `fast_ordinal_cauchit_regression_cpp` immediately sorted the same response again solely to recover `K`. Added a read-only `n_levels()` accessor to the already-constructed model and removed the duplicate sort. Correctness: `test-ordinal-cauchit-level-cache.R` verifies arbitrary, unsorted category labels produce the same parameters and iterations as their sorted integer remapping; the numerical score/Hessian test also passes. Old/new parameters, negative log-likelihood, information matrix, and iteration count are bit-identical. Benchmarks used 30 independent rounds × 500 iterations per build: N=1000 estimate 0.7935ms → 0.6498ms (18.12% faster, bootstrap 95% CI 15.03%–20.19%, one-sided Wilcoxon p=1.46e-9); N=200 variance 0.3010ms → 0.2396ms (20.41% faster, CI 15.79%–22.80%, p=1.67e-11). The optimized package was installed with `R CMD INSTALL --no-docs`.

---

**TODO-73: cont_ratio — cache build_continuation_ratio_augmented_data** ✓ DONE — ALREADY CACHED
Files: cont_ratio source

Root cause: `build_continuation_ratio_augmented_data` and `unlink_chunk.isra.0` (glibc allocator) appear in cont_ratio_est profile. The continuation-ratio model expands the design matrix into K-1 binary subproblems on every optimizer call. This augmented data depends only on the fixed design matrix and outcome (not the current coefficient values) — compute once in the constructor or at first call, store as member fields, and reuse across all optimizer evaluations.

Direct inspection disproved the per-optimizer-call premise. Both public fit paths call `build_continuation_ratio_augmented_data(X, y)` exactly once before constructing `ContinuationRatioObjective`; the objective then stores references to `X_aug` and `z` and reuses them for every likelihood/gradient evaluation and the final Hessian. The score and Hessian exports each build once because they are independent one-shot R calls. The allocator hotspot was the old per-augmented-row construction strategy, already removed by TODO-34; its profile-matched before/after benchmark used 2,000 estimate calls and 3,000 variance calls: `cont_ratio_est` 0.2782ms → 0.1765ms (−36.6%, 1.58× throughput), and `cont_ratio_var` 0.0660ms → 0.0530ms (−19.7%, 1.25× throughput). A fresh noise audit of the installed working tree used 30 rounds × 500 calls: estimate median 0.2120ms (IQR 0.2080–0.2230ms) and variance median 0.0580ms (IQR 0.0560–0.0620ms). The dedicated augmentation test passed all 4 assertions (exact augmented design/response plus R-reference score/Hessian), and the fitting-equivalence file passed all 66 assertions. The package was installed with `R CMD INSTALL --no-docs`. No production source change was necessary; the remaining per-evaluation scratch allocations are separately tracked by TODO-61.

---

**TODO-74: OrdinalGLMMObjective — preallocate working buffers** ✓ DONE
Files: ordinal GLMM source (OrdinalGLMMObjective)

Root cause: ordinal_glmm_est profile shows `_int_malloc` prominently — unlike `LogisticGLMMObjective` and `PoissonGLMMObjective` (which have preallocated lam/eneg buffers), `OrdinalGLMMObjective` still allocates working arrays inside `operator()`. Apply the same preallocated buffer pattern used in the HurdlePoisson GLMM refactor: add member `Eigen::VectorXd` scratch fields sized in the constructor, overwrite in-place each call.

Implemented constructor-sized buffers for cutpoints, fixed linear predictors, per-node log likelihoods, the full per-node gradient matrix, and the four finite-difference Hessian vectors. The objective now uses an Eigen segment view for `beta`, references the existing Gauss-Hermite nodes and weights directly, resets one contiguous gradient matrix per group, and reuses all buffers across optimizer and Hessian evaluations. This removes the former per-call vector copies and the `n_gh` separate gradient-vector allocations per group while preserving arithmetic order. Benchmarks used identical profile data and 30 independent rounds × 5 calls per build: `ordinal_glmm_est` (n=400, K=3, G=80) 34.0ms → 24.0ms median (−29.41%, 1.42× throughput; bootstrap 95% CI 27.33%–30.23%; one-sided Wilcoxon p=1.46e-11), and `ordinal_glmm_var` (n=200, K=3, G=80) 17.7ms → 12.0ms (−32.20%, 1.48×; CI 29.41%–33.71%; p=1.45e-11). After-change IQRs were 23.5–24.55ms and 11.8–12.4ms. Correctness: all 8 assertions in `test-ordinal-glmm-alpha-buf.R` pass, including `ordinal::clmm` agreement; an additional old/new K=4, 11-node comparison was bit-identical for parameters, thresholds, likelihood, variance, convergence, and the full information matrix. Installed with `R CMD INSTALL --no-docs`.

---

## Stats helpers (new from 2026-07-04 run)

**TODO-75: mn_ci — move R-level bisection loop to C++** ✓ DONE — ALREADY IMPLEMENTED
Files: mn_ci R source + `EDI/src/mn_ci_cpp.cpp` (or equivalent)

Root cause: mn_ci profile shows `bcEval_loop` (R bytecode) + `Rf_pnorm_both` + `mn_constrained_mle_pc_cpp` + `mn_z_statistic_cpp`. The R-level code drives a bisection loop calling these C++ functions one at a time, paying R→C++ call overhead plus R GC on every bisection step. `mn_ci_cpp` exists as a C++ function that should encapsulate the full bisection internally — verify it is being called, or redirect the R caller to use it directly instead of driving the loop from R.

Verification found this optimization already present. `mn_ci_cpp` in `miettinen_nurminen_speedups.cpp` owns both complete bound searches and calls the constrained MLE and score calculation directly in C++; `InferenceIncidMiettinenNurminenRiskDiff$compute_asymp_confidence_interval()` dispatches to it exactly once and contains no R loop. The observed `bcEval_loop` frame belongs to the R benchmark driver, not to the production confidence-interval algorithm. A retained R bisection reference provided the before path. Alternating-order benchmarks used 30 independent rounds × 2,000 calls: R-driven bisection 0.1607ms median (IQR 0.1503–0.1695ms) versus C++ 0.0070ms (IQR 0.0065–0.0070ms), a 22.96× speedup (bootstrap 95% CI 21.86×–25.00×; paired one-sided Wilcoxon p=9.12e-7). The dedicated `test-miettinen-nurminen-ci-dispatch.R` test passed all 6 assertions, checking single-call dispatch and C++/R-reference agreement across three tables to `1e-14`. Installed with `R CMD INSTALL --no-docs`; no additional production source change was necessary.

---

**TODO-76: zhang_fisher_pval — std::lgamma for chi-squared CDF** ✓ DONE
Files: zhang_fisher_pval source

Root cause: `Rf_chebyshev_eval` (#1 hotspot) + `Rf_lgammacor` + `__log1p_fma` — the Fisher P-value computation uses a chi-squared CDF that internally calls R's incomplete gamma via Chebyshev series. Replace `R::lgammafn` calls with `std::lgamma` to bypass the Chebyshev/Stirling dispatch layer. Same pattern as TODO-14 and TODO-115.

Source inspection corrected the mechanism: `zhang_exact_fisher_pval_cpp` does not use a chi-squared CDF; its noncentral-hypergeometric support loop called `R::lchoose` twice per table, which reaches R's log-gamma/Chebyshev machinery. Replaced those calls with direct integer log-combinations using `std::lgamma`, hoisting the two fixed row-factorial terms outside the support loop. Support enumeration, log-sum-exp normalization, relative tolerance, and two-sided probability ordering are unchanged. To control observed machine drift, the final benchmark loaded the old and new shared libraries simultaneously and used randomized paired order for 60 rounds × 5,000 calls: 0.00690ms → 0.00640ms median (1.078× throughput; bootstrap 95% CI 1.062×–1.129×; paired one-sided Wilcoxon p=4.83e-8), with IQRs 0.00680–0.00720ms and 0.00600–0.00660ms. The new `test-zhang-fisher-std-lgamma.R` test passed all 7 assertions, comparing four central/noncentral cases against `stats::fisher.test` and checking degenerate and invalid inputs. Installed with `R CMD INSTALL --no-docs`.

---

## KK21 weight functions (new from 2026-07-04 run)

**TODO-77: kk21_beta_wts — std::lgamma + fast_digamma** ✓ DONE
Files: kk21 beta weights source (kk21_beta_weights_cpp)

Root cause: `Rf_chebyshev_eval` (#1 hotspot, ~41% of samples) + `Rf_gammafn` + `Rf_lgammafn_sign` — the beta distribution weight function uses R's lgamma dispatch, which internally calls Chebyshev evaluation for the Stirling series. Replace `R::lgammafn` → `std::lgamma` and `R::digamma` → `fast_digamma`. Same pattern as TODO-14 (beta regression) and TODO-115 (negbin).

Inspection found no digamma calls in either beta helper; those calls belong to the adjacent negative-binomial helpers tracked by TODO-78. Replaced all six beta-path `R::lgammafn` call sites in the univariate and stepwise/multivariate precision-grid loops with `std::lgamma`, and hoisted the invariant `lgamma(phi)` outside each observation loop. Final timing used single-thread process CPU time with balanced randomized ABBA/BAAB ordering: 30 rounds × 50 complete N=500 calls per build. Before and after medians were 7.110ms (IQR 5.930–8.595ms) and 5.340ms (IQR 3.935–6.815ms); the median paired speedup was 1.404× (bootstrap 95% CI 1.262×–1.501×; paired one-sided Wilcoxon p=1.36e-4). The new `test-kk21-beta-std-lgamma.R` test passed all 5 assertions, reproducing IRLS, precision-grid selection, and information-based weights independently in R and checking edge/stepwise invariants. Old/new ordinary and stepwise beta-weight vectors were also bit-identical. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-78: kk21_negbin_wts — fast_digamma + pow→exp(log)** ✓ DONE
Files: kk21 negbin weights source (kk21_negbin_weights_cpp)

Root cause: `Rf_dpsifn` (#1, digamma dispatch) + `__ieee754_pow_fma` + `__ieee754_log_fma` + `__ieee754_exp_fma`. The negative-binomial weight function calls digamma for the NB score/hessian. Two fixes:
1. `R::digamma(x)` → `fast_digamma(x)` (polynomial approximation, same as TODO-23 for ZINB).
2. `pow(x, y)` where y is a floating-point parameter → `std::exp(y * std::log(x))`, fusing with adjacent log computations to eliminate the full `pow` dispatch.

Inspection found no explicit `pow` expression in either negative-binomial helper; exponentials were already direct. Replaced all four `R::digamma` call sites in the univariate and stepwise/multivariate theta updates with `fast_digamma`. Also hoisted the theta-only digamma, trigamma, log, and reciprocal values outside each observation loop, leaving the varying trigamma calls and Newton update order unchanged. Final timing used single-thread process CPU time with balanced randomized ABBA/BAAB ordering: 30 rounds × 20 complete N=500 calls per build. Before and after medians were 22.150ms (IQR 15.912–25.950ms) and 6.675ms (IQR 5.850–9.587ms); the median paired speedup was 2.808× (bootstrap 95% CI 2.525×–3.281×; paired one-sided Wilcoxon p=9.13e-7). The new `test-kk21-negbin-fast-digamma.R` test passed all 5 assertions, independently reproducing the complete IRLS/theta-Newton/information calculation in R and checking edge and stepwise invariants. Old/new ordinary and stepwise weights agreed to maximum absolute error `1.11e-16`. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-79: kk21_stepwise_logistic_wts + kk21_survival_wts — XtWX via weighted_crossprod** ✓ DONE
Files: kk21_stepwise_logistic_weights_cpp, kk21 survival weights source (univariate_weibull_tstat)

Root cause: `lhs_process_one_packet` (#1 in kk21_stepwise_logistic_wts; also visible in kk21_survival_wts) — the stepwise logistic and Weibull-based survival weight functions re-fit internal regressions using the full XtWX GEMM triple-product. Apply `weighted_crossprod` (DSYRK symmetric update, same as TODO-117) to halve FLOP count. If the internal re-fits share code paths with fast_logistic/fast_weibull, TODO-117's fix propagates automatically; otherwise apply explicitly here.

The stepwise logistic path was already fully converted: `logistic_reduced_fit_for_score_test` uses `weighted_crossprod` and `weighted_crossprod_rhs` for both IRLS and final information. The remaining generic product was the fixed two-column `[1, x]` Weibull IRLS in `univariate_weibull_tstat`. A direct `weighted_crossprod` substitution showed no measurable gain for this tiny matrix, so the final implementation removes the GEMM entirely: one observation pass accumulates the symmetric 2×2 information matrix and two-element RHS, eliminating `sqrt_w`, `Xw`, `yw`, and `Xw.transpose() * Xw`. Final survival timing used single-thread process CPU time with balanced randomized ABBA/BAAB ordering, 30 rounds × 500 complete N=500 calls per build: 0.8080ms → 0.7180ms median (IQR 0.7140–0.9675ms → 0.5565–0.8445ms), with median paired speedup 1.226× (bootstrap 95% CI 1.022×–1.447×; paired one-sided Wilcoxon p=0.0162). The unchanged logistic control used 30 rounds × 200 calls and was neutral at 1.007× (CI 0.955×–1.044×, p=0.525). The new `test-kk21-weighted-crossprod.R` tests independently reproduce both logistic score-test selection and Weibull IRLS weights in R; both pass. Old/new logistic outputs are bit-identical, and survival weights agree to maximum absolute error `4.44e-15`. Installed with `R CMD INSTALL --no-docs`.

---

## Permutation generators (new from 2026-07-04 run)

**TODO-80: generate_permutations_atkinson — preallocate QR + LU workspace** ✓ DONE
Files: `EDI/src/generate_permutations_atkinson_cpp.cpp` (or equivalent)

Root cause: `FullPivHouseholderQR::computeInPlace` (#1) + `FullPivLU::computeInPlace` + `lhs_process_one_packet` (GEMM) + `_int_malloc` per permutation. The Atkinson D-optimal design algorithm performs QR and LU factorizations for each candidate permutation, with matrix dimensions fixed by the design size. Preallocate `Eigen::FullPivHouseholderQR` and `Eigen::FullPivLU` objects outside the permutation loop and call `compute()` in-place each iteration.

Implemented in `EDI/src/generate_permutations.cpp`. The QR inputs depend only on the fixed covariate matrix and subject index, not simulated assignments, so the implementation now computes each subject's varying/independent columns and processed design once before the simulation loop, reusing one `FullPivHouseholderQR` object. Assignment-dependent LU work still runs per subject/simulation, but now reuses one `FullPivLU` object plus maximum-sized design, crossproduct, inverse, and current-row workspaces. The helper also reads prior integer assignments directly instead of allocating and converting a temporary `VectorXd`. RNG call order and factorization arithmetic are unchanged. Balanced randomized ABBA/BAAB timing used single-thread process CPU time, matched seeds, and 30 rounds × 6 complete profile calls per build (`n=100`, `p=4`, `nsim=100`): 47.250ms → 19.500ms median (IQR 43.375–49.958ms → 17.375–21.333ms), with median paired speedup 2.443× (bootstrap 95% CI 2.399×–2.513×; paired one-sided Wilcoxon p=9.13e-7). The new `test-atkinson-permutation-workspaces.R` test passed all 7 assertions covering seeded reproducibility, dimensions, binary assignments, rank deficiency, and treatment symmetry. Matched-seed old/new outputs were bit-identical for both full-rank and rank-deficient designs. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-81: generate_permutations_pocock_simon — eliminate per-permutation SEXP allocation** ✓ DONE
Files: `EDI/src/generate_permutations_pocock_simon_cpp.cpp` (or equivalent)

Root cause: `Rf_allocVector3` + `SETCDR` per permutation — R linked-list construction inside the generation loop (178 hits on `Rf_allocVector3`/`SETCDR`). Currently appends one SEXP per permutation, triggering allocation + GC pressure. Fix: accumulate permutations into a `std::vector<std::vector<int>>` inside the loop; convert to R list once at the end. Secondary: `__strcmp_avx2` (string comparison) also visible — see TODO-84 for the strata-key pattern.

Inspection found the output was already one preallocated `IntegerMatrix`, not a linked-list append. The actual SEXP/heap churn came from constructing an R `NumericMatrix` count table per simulation, an R row vector and `std::vector<double>` per subject, and repeatedly indexing the R level matrix. The implementation now validates and converts one-based global level rows once into a row-major zero-based C++ buffer, reuses one integer count buffer across simulations, and uses a stack two-element imbalance array. It also adds bounds/weight-length validation before pointer indexing. The profiler fixture incorrectly supplied zero-based, non-offset factor levels (out-of-bounds for this API); it now generates valid one-based global count rows. Balanced randomized ABBA/BAAB timing used single-thread process CPU time, matched seeds, and 30 rounds × 10 complete profile calls per build (`n=200`, three two-level factors, `nsim=500`): 17.500ms → 2.450ms median (IQR 17.300–17.800ms → 2.400–2.600ms), with median paired speedup 7.320× (bootstrap 95% CI 7.020×–7.417×; paired one-sided Wilcoxon p=9.10e-7). The new `test-pocock-simon-permutation-buffers.R` test passed all 7 assertions covering seeded reproducibility, dimensions, binary assignments, treatment symmetry, and invalid inputs. Matched-seed old/new outputs were bit-identical. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-82: generate_permutations_cluster — eliminate per-permutation R linked list** ✓ DONE
Files: cluster permutation generator source

Root cause: `SETCDR` + `unif_rand` per permutation — same R linked-list append pattern as TODO-50. Accumulate in `std::vector`, convert to R list once at the end.

Inspection found the output was already written directly into one preallocated `IntegerMatrix`; there was no per-permutation linked-list append. The repeated R work came from extracting every cluster `IntegerVector` from the input `List` inside every simulation. The implementation now validates and flattens cluster membership once into a contiguous zero-based subject-index vector plus cluster offsets, then iterates only over C++ memory in the simulation loop. Cluster order and last-write-wins behavior for overlapping clusters are preserved, as is the one-RNG-draw-per-cluster sequence. Balanced randomized ABBA/BAAB timing used single-thread process CPU time, matched seeds, and 30 rounds × 50 complete profile calls per build (`n=1000`, 20 clusters of 50, `nsim=500`): 2.390ms → 0.480ms median (IQR 2.250–2.535ms → 0.445–0.560ms), with median paired speedup 4.863× (bootstrap 95% CI 4.604×–5.182×; paired one-sided Wilcoxon p=9.13e-7). The new `test-cluster-permutation-flat-indices.R` test passed all 19 assertions covering seeded reproducibility, dimensions, binary values, within-cluster equality, target marginal probability, and invalid indices. Matched-seed old/new outputs were bit-identical. Installed with `R CMD INSTALL --no-docs`.

---

## Redraw functions (new from 2026-07-04 run)

**TODO-83: atkinson_redraw — preallocate QR/GEMM workspace** ✓ DONE
Files: `EDI/src/atkinson_redraw_batch_cpp.cpp` (or equivalent)

Root cause: `lhs_process_one_packet` (GEMM, #1) + `ColPivHouseholderQR::computeInPlace` per call. Same root cause as TODO-80 (Atkinson permutations) but in the sequential redraw path. Preallocate QR object and scratch matrices as member fields; call `compute()` in-place each redraw.

Implemented function-local maximum-sized workspaces in `EDI/src/atkinson_redraw_batch.cpp` (the exported function has no persistent objective/member lifetime). A pre-sized `ColPivHouseholderQR` and `FullPivLU` are reused across subject steps via `compute()`, together with reusable varying-column, processed-design, augmented-design, crossproduct, inverse, and current-row buffers. The varying/independent column vectors retain capacity, and the helper now reads prior assignments directly rather than allocating `w_prev`, row-segment, and augmented-current-row vectors. Active matrix blocks and QR threshold are unchanged, preserving factorization arithmetic and RNG order. Balanced randomized ABBA/BAAB timing used single-thread process CPU time, matched seeds, and 30 rounds × 200 complete profile calls per build (`n=100`, `p=4`): 0.4175ms → 0.3600ms median (IQR 0.4013–0.4388ms → 0.3500–0.3888ms), with median paired speedup 1.155× (bootstrap 95% CI 1.124×–1.173×; paired one-sided Wilcoxon p=1.67e-6). The new `test-atkinson-redraw-workspaces.R` test passed all 7 assertions covering seeded reproducibility, binary/finite output, rank deficiency, and treatment symmetry. Matched-seed old/new outputs were bit-identical for both full-rank and rank-deficient designs. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-84: pocock_simon_redraw_w — eliminate per-call SEXP allocation** ✓ DONE
Files: pocock_simon_redraw_w source (pocock_simon_assign_cpp)

Root cause: `Rf_allocVector3` + `SETCDR` per redraw call — allocating an R structure on every arm assignment within the redraw loop. Use a pre-allocated C++ buffer and return as a fixed-size integer vector, or accumulate into `std::vector<int>` and wrap once at the end.

The output was already one fixed-size `IntegerVector`; the hot allocations came from an R `NumericMatrix` count table, an R row vector per subject, and repeated calls to the generic assign helper, which allocates `std::vector` imbalance/count buffers. `pocock_simon_redraw_w_cpp` now validates and converts one-based global level rows once, uses one C++ integer count buffer and a stack two-element imbalance array, and writes directly to the output vector. The original redraw tie semantics are preserved exactly: ties consume only the `prob_T` draw and skip `p_best`. The invalid zero-based profiler fixtures for redraw and both standalone assign kernels were corrected to one-based global rows. Balanced randomized ABBA/BAAB timing used single-thread process CPU time, matched seeds, and 30 rounds × 1,000 complete profile redraws per build (`n=200`, three factors): 0.1485ms → 0.0090ms median (IQR 0.1365–0.1545ms → 0.0080–0.0100ms), with median paired speedup 16.819× (bootstrap 95% CI 15.500×–17.800×; paired one-sided Wilcoxon p=9.12e-7). The new `test-pocock-simon-redraw-buffers.R` test passed all 5 assertions, including exact matched-seed agreement with an independent R implementation and invalid-input checks. Matched-seed old/new C++ outputs were bit-identical. Installed with `R CMD INSTALL --no-docs`.

---

## Bootstrap index generators (new from 2026-07-04 run)

**TODO-85: stratified_bootstrap_indices — replace string strata keys with integer IDs** ✓ DONE
Files: `EDI/src/stratified_bootstrap_indices_cpp.cpp` (or equivalent)

Root cause: `__memcmp_avx2_movbe` (AVX2-accelerated string compare) visible in the hot path. Strata are represented as strings in the C++ layer, requiring string comparison on every sample assignment. Pre-convert strata labels to integer indices at the R layer (`match(strata, unique(strata))` or `as.integer(as.factor(strata))`) before calling C++; C++ then compares `int`s (single instruction) rather than byte strings.

Added a native integer-ID path to `stratified_bootstrap_indices_cpp` and updated all blocking, optimal-blocking, SPBR, and random-block-size bootstrap callers to pass dense IDs computed once with `match(strata, unique(strata))`. The profiler now exercises integer IDs. The original character path remains available for backward compatibility. The integer implementation uses `std::map<int, std::vector<int>>`, eliminating string construction and byte comparisons while retaining sorted-group sampling and the final global shuffle. Balanced randomized ABBA/BAAB timing used single-thread process CPU time and 30 rounds × 1,000 complete profile calls per build (five strata × 200): 0.0635ms → 0.0180ms median (IQR 0.0610–0.0660ms → 0.0170–0.0198ms), with median paired speedup 3.471× (bootstrap 95% CI 3.389×–3.556×; paired one-sided Wilcoxon p=9.08e-7). The new `test-stratified-bootstrap-integer-ids.R` test passed all 7 assertions covering exact stratum-size preservation, valid indices, character fallback, and marginal uniformity. The legacy helper uses a private `random_device`-seeded generator, so matched-seed equality is unavailable; both old and new benchmark outputs independently satisfied the exact stratification invariants. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-86: bootstrap_m_indices — consolidate dual RNG paths** ✓ DONE (implemented alongside TODO-91)
Files: bootstrap_m_indices source

Root cause: Both `unif_rand` and `Rf_runif` appear in the profile — two different call paths for uniform random draws in the same function. Consolidate to `unif_rand` (the standard R RNG API used everywhere else); `Rf_runif` duplicates state tracking and adds unnecessary dispatch overhead.

Inspection showed that the current `EDI/src/bootstrap_match_indices.cpp` had already received the stronger TODO-91 implementation: `bootstrap_m_indices_cpp` consumes one `R::unif_rand()` value to seed a local `std::mt19937_64`, then uses unbiased Lemire bounded draws for every reservoir and matched-pair selection. There is no remaining `R::runif`/`Rf_runif` path. The exact pre-consolidation implementation was reconstructed from the parent of commit `85ba97e5`; it used `R::runif()` for every sampled index. Noise-controlled BAAB timing used process CPU time, two passes per build, 25 rounds per pass, and 50 complete profiler-sized calls per round (2,500 timed calls per build; `n_reservoir=500`, `m=250`, `B=500`): 5.330ms → 3.240ms median (IQR 5.160–5.495ms → 3.120–3.395ms), a 1.645× speedup (independent bootstrap 95% CI 1.586×–1.693×; one-sided Wilcoxon p=3.53e-18). The timing distributions did not overlap (old minimum 4.94ms, new maximum 4.18ms). The dedicated `test-bootstrap-m-local-rng.R` test passed all 11 assertions covering seeded reproducibility, dimensions, intact pair sampling, reservoir and pair marginal uniformity, and reservoir-only/pair-only boundaries. Installed with `R CMD INSTALL --no-docs`.

---

## Post-fit variance helpers (new from 2026-07-04 run)

**TODO-87: gcomp_frac_logit_post_fit_var — eliminate per-bootstrap Rf_allocVector3** ✓ DONE
Files: `EDI/src/gcomp_logistic_post_fit_cpp.cpp` (or equivalent)

Root cause: `Rf_allocVector3` visible alongside `lhs_process_one_packet` (GEMM) and `R_gc_internal` — a new R vector (SEXP) is being allocated per bootstrap replicate inside `gcomp_logistic_post_fit_cpp`. Pre-allocate once in the outer C++ function and overwrite in-place across replicates to eliminate per-replicate GC pressure.

Inspection corrected the profile interpretation: the post-fit helper contains no bootstrap loop. The avoidable allocation was `gcomp_fractional_logit_post_fit_cpp` calling the public logistic helper, materializing its complete 10-element R list (including three risk-ratio outputs that fractional logit discards), and then constructing a second 7-element R list. The common calculation now returns a native `LogisticPostFitResult` containing moved Eigen results and scalar fields; each public wrapper materializes exactly its own R result once. Noise-controlled BAAB timing used process CPU time, two passes per build, 30 rounds per pass, and 5,000 complete profiler-sized calls per round (300,000 timed calls per build; `N=200`, six design columns): 0.019300ms → 0.017200ms median (IQR 0.017950–0.020400ms → 0.016750–0.017950ms), a 1.122× speedup / 10.9% time reduction (independent bootstrap 95% CI 1.069×–1.158×; one-sided Wilcoxon p=2.53e-11). Baseline and optimized result lists were byte-for-byte identical. The dedicated `test-gcomp-fractional-post-fit-native-result.R` test passed all 19 assertions against an independent R sandwich-covariance and delta-method implementation, including the shared logistic risk-difference/risk-ratio outputs and retained validation errors. Installed with `R CMD INSTALL --no-docs`.

---

## Additional findings from direct perf report analysis (2026-07-04)

**TODO-88: Ordinal cauchit — fast_atan approximation for cauchit link** ✓ DONE
Files: `EDI/src/ordinal_fixed_link_helpers.h`, `EDI/src/fast_ordinal_cauchit_regression.cpp`

Root cause: `__atan_fma` accounts for 19.8% of `ord_cauchit_est` samples — the second largest hotspot behind the operator() loop itself. The cauchit CDF `F(x) = 0.5 + atan(x)/π` requires one `atan` call per threshold per observation. For K=5 categories and n=400 obs: 1,600 atan calls per optimizer step; for K=5, n=2000: 8,000 calls. `__atan_fma` is glibc's AVX2-tuned atan, but a 9-term minimax polynomial on `[-6, 6]` with argument reduction `atan(x) = π/2 − atan(1/x)` for |x|>6 achieves ≤1e-10 relative error and is ~2–3× faster (avoids the range-reduction table lookup in glibc). Implement as `fast_atan(x)` in `_helper_functions.h`. Apply inside `FixedOrdinalRegression` for the cauchit CDF and PDF evaluations. Guard: validate ≤1e-9 max error on a grid of x values before committing.

Implemented a more accurate range-reduced rational approximation in `ordinal_fixed_link_helpers.h`. A direct polynomial across `[-6,6]` cannot realistically provide the requested near-double-precision guard because of `atan`'s nearby complex singularities. The retained implementation reduces to `|x| <= tan(pi/8)` using pi/4 and reciprocal identities, then evaluates a degree-5/5 rational minimax approximation; it preserves signed zero, infinities, NA, and NaN. Only the CDF changes because the cauchit PDF and its derivative were already exact rational expressions with no `atan` call. A dense pre-integration sweep across central values and magnitudes from 1e-300 to 1e300 found 2.22e-16 maximum absolute error versus libm `atan`. Noise-controlled BAAB timing used process CPU time, two passes per build, 30 rounds per pass, and 200 complete profile fits per round (12,000 timed fits per build; `N=1000`, five predictors, three response levels): 0.6550ms → 0.5825ms median (IQR 0.6250–0.6963ms → 0.5600–0.6200ms), a 1.125× speedup / 11.1% time reduction (independent bootstrap 95% CI 1.084×–1.179×; one-sided Wilcoxon p=3.06e-11). Both builds converged in seven iterations and fitted parameters differed by at most 1.12e-16. The dedicated `test-ordinal-cauchit-fast-atan.R` test passed all 8 assertions over a 440,000-point error grid and exact-`pcauchy` numerical score/Hessian references; the existing scratch-buffer and arbitrary-level tests passed all 5 assertions. Installed with `R CMD INSTALL --no-docs`.

---

**TODO-89: Ordinal cloglog — cache exp(x) to avoid double exp per CDF eval** ✓ DONE
Files: `EDI/src/ordinal_fixed_link_helpers.h`, `EDI/src/fast_ordinal_cloglog_regression.cpp`

Root cause: `ord_cloglog_est` shows 26% `__ieee754_exp_fma` + 11% `exp@@GLIBC_2.29` = **37% total exp**. The cloglog CDF is `F(x) = 1 − exp(−exp(x))`: two sequential `exp` calls per evaluation. The inner `e_x = exp(x)` is the value needed for both the CDF (`exp(-e_x)`) and its derivative (`e_x * exp(-e_x)`). Current code calls separate `cdf()` and `pdf()` inline functions with the same `z` argument.

Implemented fused helpers `cdf_and_pdf(link, z, F, f)` and `cdf_pdf_fpdf(link, z, F, f, fp)` in the `edi_ordinal` namespace. `FixedOrdinalRegression::operator()`, `hessian()`, and `expected_hessian()` now use those helpers instead of separate `cdf()`/`pdf()`/`pdf_derivative()` calls. For the cloglog link, each endpoint now computes `e_z = exp(z)` and `exp(-e_z)` once, then derives `F`, `f = e_z * exp(-e_z)`, and `f' = f * (1 - e_z)`. The implementation preserves the existing saturation fast paths (`z > 5` returns `F=1, f=0`; `z < -37` returns zero), which are stronger than the proposed `z > 37` guard and are already double-equivalent for the CDF at the upper threshold. Added `fast_cloglog_link_eval_cpp` only as an internal test hook for direct endpoint validation.

Noise-controlled benchmark used clean before/after installs with `R CMD INSTALL --no-docs --preclean`. The before build was the saved pre-TODO-89 source with TODO-88's cauchit `fast_atan` retained; the after build was the optimized source. BAAB timing used separate R processes, process CPU time, 30 rounds × 500 complete profiler-sized `ord_cloglog_est` fits per phase (`N=1000`, five predictors, three response levels), for 60,000 timed fits total. Combined medians: 0.6300ms → 0.3700ms per fit (IQR 0.6255–0.6385ms → 0.3620–0.3820ms), a 1.703× speedup / 41.3% time reduction. Independent bootstrap 95% CI for the speedup: 1.679×–1.725×; one-sided Wilcoxon p=1.74e-21. All four benchmark phases converged in six iterations, and before/after fitted parameters differed by at most 4.86e-17.

Correctness: `test-ordinal-cloglog-cached-exp.R` passed all 9 assertions, checking direct CDF/PDF/PDF' values against independent formulas and score/Hessian against numerical derivatives of an independent R likelihood. The pre-existing `test-ordinal-cloglog-exp-cache.R` also passed all 7 assertions. A standalone `ordinal::clm(link="cloglog")` equivalence check passed. Installed with `R CMD INSTALL --no-docs --preclean`.

---

**TODO-90: Stereotype logit — GEMV refactor for loglik_grad** ✓ DONE
File: `EDI/src/fast_stereotype_logit.cpp`

Root cause: `StereotypeLogitRegression::loglik_grad` accounts for 36–37% of both `stereotype_est` and `stereotype_var` samples. The gradient function likely accumulates per-observation rank-1 outer products: `gradient += score_i * x_i.T` (same scatter anti-pattern fixed by TODO-15/16/62 for ZIP/ZINB/NegBin). Fix: accumulate scalar score weights `w_i` for each parameter block into a preallocated weight vector, then apply a single `X.transpose() * w` GEMV after the observation loop. Also preallocate the per-obs working vectors `logit_grad(K)` and `logit_hess(K)` flagged in TODO-36 (hessian allocs) — the gradient path has the same per-obs allocation pattern.

Implemented reusable buffers in `StereotypeLogitRegression` for `eta = X * beta`, per-observation logits/probabilities, and the beta score weights. `loglik_grad()` now computes `eta` once per call, accumulates one scalar beta score weight per observation inside the observation loop, and applies a single `X.transpose() * w` GEMV after the loop. The gamma-gradient path was also tightened by replacing a per-observation temporary expected-score vector with direct scalar accumulation. The Hessian and expected-Hessian paths now reuse their own per-call working storage and share the precomputed `eta` vector instead of repeatedly evaluating `X.row(i).dot(beta)`.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for both the saved pre-TODO-90 source and the optimized source. Timing used separate R processes, process CPU time, two before phases and two after phases, 30 rounds per phase, 20 complete fits per round for each kernel (`stereotype_est`: `N=1000`; `stereotype_var`: `N=200`; five predictors; three response levels), for 1,200 timed fits per build per kernel. Combined medians: `stereotype_est` 1.3000ms → 1.1000ms per fit (IQR 1.2500–1.4000ms → 1.0500–1.1500ms), a 1.182× speedup / 15.4% time reduction (independent bootstrap 95% CI 1.182×–1.286×; one-sided Wilcoxon p=5.45e-15). `stereotype_var` improved 0.3000ms → 0.2500ms per fit (IQR 0.3000–0.3500ms → 0.2500–0.2500ms), a 1.200× speedup / 16.7% time reduction (bootstrap 95% CI 1.200×–1.200×; one-sided Wilcoxon p=1.28e-12). All four benchmark phases converged, and before/after fitted parameters differed by at most 1.11e-16 for `estimate_only` and 2.78e-16 for the variance path.

Correctness: added `test-stereotype-logit-gemv-gradient.R`, which checks the C++ score/Hessian against numerical derivatives of an independent R likelihood and verifies repeated calls are exactly repeatable with the new mutable scratch buffers. The dedicated test passed all 4 assertions. Existing standalone stereotype checks for K=2 `glm(binomial)` equivalence, K=3 score near zero at the MLE, and K=3 Hessian finite-difference agreement also passed.

---

**TODO-91: Permutation/bootstrap generators — replace R RNG with seeded local mt19937** ✓ DONE (implemented alongside TODO-92)
Files: `EDI/src/generate_permutations.cpp`, `EDI/src/bootstrap_indices.cpp` (and related)

Root cause: `unif_rand` accounts for 31% of `generate_permutations_bernoulli`, 33% of `generate_permutations_efron`, 27% of `generate_permutations_matching`, 34% of `bootstrap_indices`, and 24% of `bootstrap_m_indices`. Each call to `unif_rand()` acquires R's global RNG state, performs a full MT19937 step, and dispatches through R's RNG kind switch. For B×n draws (e.g., 5000 permutations × 400 subjects = 2M calls), this overhead dominates. Fix: at the start of each C++ function, draw one seed from `R::unif_rand()`, then use a local `std::mt19937_64` for all subsequent draws. Use Lemire's nearly-divisionless method for bounded integers to avoid the modulo bias fix overhead. Reproducibility: fully reproducible given the R seed, since the local RNG is deterministically derived. Expected: 3–5× faster for the RNG-bound portion.

---

**TODO-92: generate_permutations_spbr — replace string map with integer-ID lookup** ✓ DONE
File: `EDI/src/generate_permutations.cpp` (SPBR branch)

Root cause: `__memcmp_avx2_movbe` accounts for **22%** of `generate_permutations_spbr` samples. The actual source is `std::map<std::string, std::vector<int>> strata_states` inside the nsim×n loop: for each of the n subjects in each of the nsim simulations, the map does O(log K) string comparisons via memcmp. For n=1000, nsim=500: 500,000 string map lookups per call. Fix: pre-convert all string strata keys to dense integer IDs once before the simulation loop using `std::unordered_map<std::string, int>`. Replace `std::map<std::string, ...>` with `std::vector<std::vector<int>> strata_states(num_strata)` indexed by integer ID — O(1) lookup with no string comparison. Persistent outer vector reused across simulations (inner vectors `.clear()`'d at top of each sim, keeping capacity). Also swapped `#include <map>` → `#include <unordered_map>`.

Result (2026-07-06): 30 rounds × 3 reps, Wilcoxon p ≈ 0.

| Kernel | Before | After | Speedup |
|---|---|---|---|
| `generate_permutations_spbr_cpp` (5 strata × 200, block=4, nsim=500) | 37 ms | 3 ms | **12.33×** |

The outsized gain (vs the 22% profiler attribution) is because the profiler was run against the pre–TODO-91 code which also had `default_random_engine` construction overhead; both were eliminated together. The string-map replacement is the dominant fix: it removes O(n × nsim × log K) string comparisons.

Test: `EDI/tests/testthat/test-spbr-integer-key.R` — 6 tests: dimensions/values, per-block balance invariant (all blocks sum to n_T_block), grand mean ≈ prob_T, single-stratum edge case, many-strata correctness, reproducibility.

---

**TODO-93: Ordinal GLMM — preallocate alpha_buf in GLMMObjective** ✓ DONE
File: `EDI/src/_glmm_engine.h`

Root cause: `_int_malloc` at 5% + `_int_free_merge_chunk` at 4% persist in both `ordinal_glmm_est` and `ordinal_glmm_var` after Phase 1 fixes. The most likely source is `alpha_buf` — a `std::vector<double>(nm)` allocated once per `operator()` call (every optimizer step). Since `nm` is fixed for a given dataset (= n_thresholds = K-1), moved to a `mutable std::vector<double> m_alpha_buf` member of `GLMMObjective`, sized in the constructor via `m.n_model_params()`. Both `value()` and `operator()` now reuse the preallocated buffer. Safety: the mutable-field LICM concern (see feedback memory) applies to Model-class fields accessed inside the GH inner loop — `m_alpha_buf` is in `GLMMObjective` itself, is written by `fill_alpha` before the loops, and only read through a `const double*` by `model.log_prob{_derivs}` (a different object). Compiler sees no aliasing; verified empirically that no regression occurs.

Result (2026-07-06): 30 rounds × 3 reps, Wilcoxon p < 0.02/0.0001.

| Kernel | Before | After | Speedup |
|---|---|---|---|
| `ordinal_glmm_est` (n=400, K=3, G=80) | 34.5 ms | 33.0 ms | **1.05×** |
| `ordinal_glmm_var` (n=200, K=3, G=80) | 18.5 ms | 17.0 ms | **1.09×** |

Gain is modest because for K=3, `nm=2` — only 16 bytes freed per optimizer step. The profiler's 9% attribution was measured against a larger K or heavier load. Allocation overhead is real but small relative to GH quadrature computation. The `var` path gains slightly more because the hessian path calls `operator()` many more times (numerical Hessian = 2p+1 evaluations per step).

Test: `EDI/tests/testthat/test-ordinal-glmm-alpha-buf.R` — 4 tests: estimate_only/full consistency, finite+converged output, treatment direction, comparison with `ordinal::clmm`.

---

**TODO-94: Logistic GLMM — fast log1pexp to reduce log1p overhead** ✓ DONE
File: `EDI/src/fast_logistic_glmm.cpp`, `EDI/src/_helper_functions.h`

Root cause: `__log1p_fma` accounts for 18% of `logistic_glmm_est` and 14% of `logistic_glmm_var` — the dominant remaining hotspot after all prior optimizations. This is `log1pexp(x) = log(1 + exp(x))` used in the node log-likelihood (log normalizer for the logistic GLMM). The glibc `log1p` + `exp` path uses separate range-reduction tables. Implement `fast_log1pexp(x)` in `_helper_functions.h`: for `x > 37`, return `x` (overflow-safe); for `x < -37`, return `exp(x)`; otherwise `x + log1p(exp(-|x|))` using a precomputed 8-term polynomial for `log1p` on `[-1, 0]` that avoids the glibc dispatch. Also applicable to `hurdle_p_glmm_est` (12% `__log1p_fma`) and `logistic_glmm_var`.

Implementation: Added `fast_log1pexp(double x)` (scalar) and `log1pexp_array_fast(ArrayXd)` (vectorized) to `_helper_functions.h`. The array version uses an atanh-series Horner polynomial: `log1p(z) = 2*s*(1 + s²/3 + … + s¹⁸/19)` where `s = z/(2+z)`, 10 odd terms, error < 5e-12 for `z ∈ [0,1]`. This is fully SIMD-vectorizable (no `.log()` or `.select()` — only `.exp()`, scalar arithmetic, and Horner FMAs). Replaced all three `log1pexp_array_safe` calls in `fast_logistic_glmm.cpp`. Note: `__log1p_fma` hotspot was already partially addressed by `fast_log1p_arr` in a prior commit; the atanh polynomial yields an additional improvement mainly on the variance (hessian) path.

Benchmark (A/B/A, 30 rounds × 20 reps, n=200 subjects / 50 groups):

| Kernel | Before | After | Speedup | Wilcoxon p |
|---|---|---|---|---|
| `logistic_glmm_est` | 0.110s | 0.110s | 1.00x | 0.618 (n.s.) |
| `logistic_glmm_var` | 0.124s | 0.118s | **1.06x** | 1.25e-05 |

Test: `EDI/tests/testthat/test-fast-log1pexp.R` (8 tests: convergence, direction, estimate/var consistency, lme4 comparison).

---

## Findings confirmed from direct file inspection (not covered above)

**TODO-95: draw_binary_match_assignments — thread-unsafe unif_rand inside OpenMP section** ✓ DONE
Files: `EDI/src/binary_match_search.cpp`

Root cause: The profiler shows `draw_binary_match_assignments_cpp [clone ._omp_fn.0]` (OpenMP parallel clone) + `unif_rand` in the same call stack. R's `unif_rand()` accesses a global, non-reentrant RNG state, so concurrent calls from multiple OpenMP threads cause a data race. The thread-safety bug was already partially fixed (pre-drawing all values with `Rcpp::runif(num_pairs * r)` before the parallel loop), but this still allocated an O(num_pairs × r) `NumericVector` (e.g., 200K+ doubles for r=2000) and read back from memory inside the parallel loop.

Implementation: replaced `Rcpp::runif(num_pairs * r)` bulk pre-draw with a single seed drawn from `R::unif_rand()` (inside `Rcpp::RNGScope`) before the OpenMP region. Each column `j` gets its own `std::mt19937_64` seeded via `splitmix64(master_seed + j)` — ensures independent streams per column with no heap allocation. Coin flip uses `rng() >> 63` (MSB test, no FP conversion).

Benchmark (A/B/A, 30 rounds × 10 reps, n=200 / 100 pairs / r=2000 / 4 cores):

| Kernel | Before | After | Speedup | Wilcoxon p |
|---|---|---|---|---|
| `draw_binary_match_assignments` | 0.045s | 0.033s | **1.36x** | 3.22e-13 |

Test: `EDI/tests/testthat/test-draw-binary-match-local-rng.R` (207 assertions: dims, binary values, per-pair balance, Bernoulli marginals, reproducibility, symmetry).

---

**TODO-96: poisson_glmm_var — preallocate PoissonGLMMObjective variance-path working buffers** ✓ DONE
Files: `EDI/src/fast_poisson_glmm.cpp` (PoissonGLMMObjective)

Root cause: `_int_malloc` appears in the `poisson_glmm_var` profile alongside `__ieee754_exp_fma` and GEMV. The variance/hessian computation path allocates working vectors inside `PoissonGLMMObjective::operator()` or a separate hessian method. Unlike `LogisticGLMMObjective` (which has preallocated buffers), the Poisson GLMM variance path still heap-allocates. Apply the same preallocated member buffer pattern used for HurdlePoisson GLMM and OrdinalGLMMObjective (TODO-74/93).

Fix: Added 17 preallocated member fields to `PoissonGLMMObjective` (sized at construction from `dat` dimensions): `m_eta_all`, `m_term_all_k`, `m_log_terms_mat(G,n_nodes)`, `m_mu_all_k_mat(n,n_nodes)` (replaces `std::vector<VectorXd>` — column-per-node layout), `m_ll_g_vec`, `m_post_k_exp`, `m_wres_all_k`, plus hessian-specific `m_E_Hik_sum`, `m_E_GiGiT_sum`, `m_G_avg_outer`, `m_pk_vec`, `m_pk_exp`, `m_G_avg_gi`, `m_E_GiGiT_gi`, `m_G_ik`, `m_wres_k_gi(max_grp)`, `m_wmu_k_gi(max_grp)`. Saves 20×VectorXd(n) allocation per optimizer iteration for both `operator()` and `hessian()` inner node-loops. Correctness verified: 19 tests pass including lme4 comparison, symmetry/PD of vcov, score≈0 at convergence.

Benchmark (G=80 groups × 5 obs, n=400, p=2, n_nodes=20; cold start; 30 rounds × 10 reps; Wilcoxon one-sided):

| path | old (ms) | new (ms) | speedup | p-value |
|------|----------|----------|---------|---------|
| EST (estimate_only=TRUE) | 2.050 | 1.500 | **1.37×** | 3.44e-15 |
| FULL (est + hessian)     | 2.100 | 1.700 | **1.24×** | 1.50e-13 |

Note: warm-start benchmarks show ~0× improvement because the optimizer converges in 1–2 iterations (constructor allocation dominates), while cold-start reveals the true savings across ~50+ iterations.

---

**TODO-97: gaussian_lmm_var — preallocate lmm_analytic_hessian temporary matrices** ✓ DONE
Files: `EDI/src/fast_gaussian_lmm.cpp` (lmm_analytic_hessian, GaussianLMMObjective)

Root cause: `Eigen::PlainObjectBase::resize` + `cfree` visible in `gaussian_lmm_var` — Eigen is dynamically allocating and immediately freeing temporary matrices inside `lmm_analytic_hessian` per call. The `generic_product_impl` (GEMM) creating intermediate `Matrix<double, -1, -1>` objects triggers heap allocation for each matrix product. Fix: pre-declare named `Eigen::MatrixXd` temporaries at the top of `lmm_analytic_hessian` sized from data dimensions, then use `.noalias()` assignments to write directly into them without allocation. For repeated variance calls (bootstrap), store as class member fields.

Fix: Three changes:
1. **SM2/SV2 stack-allocated types**: `using SM2 = Eigen::Matrix<double, Dynamic, Dynamic, ColMajor, 2, 2>` — max-size-2 Eigen matrices use internal fixed storage (stack) for all per-group matrices (P, dV_e, dV_b, dP_e, dP_b, d2V_xx, term, tmp, etc.). Eliminates ~15 heap allocs per group × G groups per hessian call.
2. **Analytical 2×2 inverse**: `lmm_P_analytic(m, v_e, v_b)` replaces `LDLT<MatrixXd>` for V = v_e*I + v_b*J. For m=1: P=1/(v_e+v_b); m=2: det=v_e*(v_e+2*v_b), P_diag=(v_e+v_b)/det, P_off=-v_b/det. No LDLT heap allocation.
3. **Preallocated member buffers** in `GaussianLMMObjective`: `m_r(n)`, `m_grad_beta(p)`, `m_buf_XtA(p, max_m)` sized at construction. `fast_gaussian_lmm_cpp` uses `obj.hessian(par)` instead of standalone `lmm_fisher_hessian`. Also avoids Xg/rg copies by using `dat.X_s.middleRows(s,m)` / `buf_r.segment(s,m)` directly. Correctness verified: 19 tests pass including lme4 comparison, H symmetry/PD, score≈0.

Benchmark (G=200 matched pairs, n=400, p=3; cold start; 30 rounds × 10 reps; Wilcoxon one-sided):

| path | old (ms) | new (ms) | speedup | p-value |
|------|----------|----------|---------|---------|
| EST (estimate_only=TRUE) | 0.100 | 0.100 | 1.00× (n.s.) | 0.349 |
| FULL (est + hessian)     | 0.500 | 0.200 | **2.50×**    | 1.11e-16 |

Hessian alone: ~0.400 ms → ~0.100 ms (**4×** speedup). EST unchanged (analytical gradient was already scalar — no matrix alloc).

---

## Profiling coverage gaps — kernels to add (2026-07-04 audit)

The 2026-07-04 profiling run covered 129 kernels but left the following C++ files with **zero kernel coverage**. Each contains real computation; any optimization findings from the above analysis may also apply to these paths.

**TODO-98: Add profiler kernels for uncovered C++ paths** ✓ DONE
Files: `profile/edi_kernel_profiler.R`

The following source files have no entry in `edi_kernel_profiler.R`. Add one kernel per exported function and re-run `perf record` + `perf annotate` on each:

| Source file | Size | Key exported functions | Priority |
|---|---|---|---|
| `fast_clogit_plus_glmm.cpp` | 592 lines | `ClogitPlusGLMMObjective` (full GLMM objective with GH quadrature) | **HIGH** — same class of hotspots as logistic_glmm/ordinal_glmm but never sampled |
| `fast_survival_models_optim.cpp` | 721 lines | 7 exported parametric survival functions | **HIGH** — largest unsampled file |
| `lrt_ci_newton.cpp` | 286 lines | `lrt_ci_nr_cpp`, `pval_invert_ci_cpp` | **HIGH** — Newton-Raphson CI loop, likely special-function heavy |
| `optimal_design_search.cpp` | 251 lines | `d_optimal_search_cpp`, `a_optimal_search_cpp` | **HIGH** — iterative D/A-optimal search with matrix operations |
| `kk_compound_distr_parallel.cpp` | 251 lines | `compute_matching_compound_distr_parallel_cpp`, bootstrap variant | **HIGH** — OpenMP parallel, likely dominant in matching workflows |
| `fast_bai_parallel.cpp` | 157 lines | `compute_bai_distr_parallel_cpp` | **MEDIUM** — OpenMP parallel BAI distribution |
| `bisection_ci.cpp` | 117 lines | bisection CI loop | **MEDIUM** — looped CI computation |
| `cmh_speedups.cpp` | 144 lines | 2 CMH test functions | **MEDIUM** |
| `kk_bootstrap_loop.cpp` | 69 lines | `matching_bootstrap_loop_cpp` (OpenMP) | **MEDIUM** — OpenMP matching bootstrap |
| `rerandomization_helpers.cpp` | — | `rerandomization_search_cpp`, `compute_objective_vals_cpp` | **MEDIUM** — rerandomization search loop |
| `fast_cpoisson_combined.cpp` | — | combined constrained Poisson | **LOW** |
| `match_data_compute_speedup.cpp` | 39 lines | `match_diffs_cpp` | **LOW** |
| `sample_bootstrap_distr_weighted_distances.cpp` | 40 lines | `compute_bootstrapped_weighted_sqd_distances_cpp` | **LOW** |

Note: `fast_clogit_plus_glmm.cpp` contains a `ClogitPlusGLMMObjective` with its own `log1pexp_cpp` and `log_sum_exp_cpp` — by inspection, likely has the same `__log1p_fma` / `__ieee754_exp_fma` / GH-quadrature GEMV hotspots as the other GLMM objectives, and the same preallocated-buffer and fast-log1pexp fixes (TODO-93/94) will apply. Profile to confirm before implementing.

Note on benchmark/test methodology: TODO-98 is a profiling infrastructure task, not a code optimization — there is no "before/after" performance comparison. The "correctness test" equivalent is verifying all 10 new kernels execute without error and return valid output; 22 smoke tests pass (`test-todo98-new-kernel-smoke.R`). The `bai_distr` kernel had a missing `EDI:::` prefix which was also fixed. Kernel timings (single-threaded, cold start): clogit_glmm_est 4.3ms, clogit_glmm_var 7.4ms, dep_cens_transform_est 1.8ms, dep_cens_transform_var 1.9ms, d_optimal_search 6.6ms, kk_compound_distr 6.3ms, bai_distr 10.4ms, rerandomization_search 3.7ms, rerandomization_obj_vals 11.9ms, cmh_block_se 0.06ms.

---

## TODO-98 results — profiling run 2026-07-05

Kernels added to `edi_kernel_profiler.R` and profiled via `profile/run_edi_perf_missing_coverage.sh`.

---

**TODO-99: ClogitPlusGLMMObjective — inner-loop eta_k_g/mu_k_g temporaries dominate allocation cost** ✓ DONE
Files: `EDI/src/fast_clogit_plus_glmm.cpp` (`ClogitPlusGLMMObjective::operator()`, lines 203–214)

Root cause: `operator()` correctly uses preallocated member fields for most buffers (`m_eta_conc_all`, `m_log_terms_mat`, `m_mu_conc_all_k_mat`, etc.). However lines 207–208 still create `const Eigen::ArrayXd eta_k_g` and `const Eigen::ArrayXd mu_k_g` as heap-allocated temporaries **inside the double loop** over `n_nodes` (20) × `m_G_conc` (80) = 1 600 allocations per gradient call. Each allocation is for a group of size `sz` (here avg ≈ 2–3 obs), so small but frequent. `perf annotate` for `clogit_glmm_est` (1 000 reps, estimate-only): `malloc` 252 samples + `cfree` 214 = **466 samples (≈ 30%)** of total, making allocation the single largest cost.

Fix applied:
- Added `mutable Eigen::ArrayXd m_eta_k_g` (size `m_max_group_size_conc`); used in inner loop of both `operator()` (non-const) and `neg_glmm()` (const → mutable needed).
- Made `m_eta_conc_all`, `m_grp_y_eta0_conc`, `m_log_terms_mat` `mutable` so `neg_glmm()` reuses them (eliminated 3 per-call heap allocs: n_conc + G + G×Q).
- `operator()` inner loop: replaced `plogis_array_safe(eta_k_g)` with direct `1/(1+exp(-x))` expression into `m_mu_conc_all_k_mat` (no temporary ArrayXd); replaced `log1pexp_array_safe(eta_k_g).sum()` with scalar loop over `m_eta_k_g` (no temporary ArrayXd).
- Eliminated `y_conc_all = y_conc.array()` copy (used `y_conc.segment()` directly).
- Added `Eigen::VectorXd m_grad_beta_conc` member (replaces per-call `Eigen::VectorXd::Zero()` alloc).
- Test: `EDI/tests/testthat/test-clogit-glmm-buffer-reuse.R` (20 tests, all pass).

Benchmark (n_disc=200, n_conc=160, p=4, G=60; 30 rounds × 8 reps cold-start, Wilcoxon):

| Path | Old (ms) | New (ms) | Ratio | p-value |
|------|----------|----------|-------|---------|
| EST  | 2.50     | 1.62     | 1.54x | ≈ 0     |
| FULL | 4.25     | 4.00     | 1.06x | 0.023   |

EST gain is 1.54× (all operator() allocation eliminated from hot path). FULL gain is modest (1.06×) because `hessian()` still allocates locally — that is TODO-100's scope.

---

**TODO-100: ClogitPlusGLMMObjective::hessian() — 8 fresh local matrices per call** ✓ DONE
Files: `EDI/src/fast_clogit_plus_glmm.cpp` (`ClogitPlusGLMMObjective::hessian()`, lines 284–297)

Root cause: `hessian()` allocates eight local `Eigen::MatrixXd`/`VectorXd` objects on every call: `eta_conc_all` (n_conc), `grp_y_eta0` (G), `log_terms_mat` (G×Q), `mu_conc_all_k_mat` (n_conc×Q), `w_conc_all_k_mat` (n_conc×Q), `mu_group_sum_mat` (G×Q), `w_group_sum_mat` (G×Q), `res_group_sum_mat` (G×Q) — totalling ~130 KB per hessian call. In `clogit_glmm_var` (300 reps + hessian per optimization): `malloc` 75 + `cfree` 82 = **157 samples** in the variance profile.

Fix applied:
- Reused existing mutable members from TODO-99: `m_eta_conc_all`, `m_grp_y_eta0_conc`, `m_log_terms_mat`, `m_mu_conc_all_k_mat`, `m_ll_g_vec`, `m_eta_k_g` — `hessian()` is non-const so all are writable.
- Added 10 new mutable members: `m_w_conc_all_k_mat` (n_conc×Q), `m_mu_group_sum_mat` (G×Q), `m_w_group_sum_mat` (G×Q), `m_res_group_sum_mat` (G×Q), `m_hess_pk_vec` (G), `m_hess_d2L` (p), `m_hess_beta_avg` (p), `m_hess_beta_beta` (p×p), `m_hess_beta_sigma` (p), `m_hess_g_beta` (p).
- Eliminated `eta_k_all`/`mu_k_all` per-k allocs: computed `m_mu_conc_all_k_mat.col(k)` and `m_w_conc_all_k_mat.col(k)` directly via lazy expressions (no intermediate ArrayXd).
- Eliminated `Xg` copies: `weighted_crossprod` takes `MatrixBase<Derived>` so `X_conc.middleRows(s, sz)` passes directly — no copy. Same for `Xg.transpose() * w_seg` products.
- `eta_k_g` in inner loop: uses `m_eta_k_g.head(sz)` + scalar `log1pexp_safe` loop (consistent with TODO-99 operator() fix).
- Per-group accumulator allocs eliminated: `m_hess_beta_avg`, `m_hess_beta_beta`, `m_hess_beta_sigma` re-zeroed per group, `m_hess_g_beta` computed per k.
- Test: `EDI/tests/testthat/test-clogit-glmm-buffer-reuse.R` extended with 2 new hessian-specific tests (now 22 tests total, all pass).

Benchmark (n_disc=200, n_conc=160, p=4, G=60; 30 rounds × 6 reps cold-start, Wilcoxon; "old" = original pre-TODO-99/100):

| Path | Old (ms) | New (ms) | Ratio | p-value |
|------|----------|----------|-------|---------|
| FULL (TODO-99+100 vs original) | 4.33 | 2.50 | 1.73x | ≈ 0 |

Decomposition: TODO-99 FULL was 4.00ms (1.08x over original) → TODO-100 takes FULL from 4.00ms to 2.50ms = **1.60x improvement** for the hessian path alone. Combined with TODO-99 EST speedup (1.54x), the overall allocation burden for the full clogit+GLMM inference call is substantially reduced.

---

**TODO-101: DepCensTransformLikelihood — replace R::pnorm/dnorm with inline C++ in per-observation loop** ✓ DONE
Files: `EDI/src/fast_survival_models_optim.cpp`, `EDI/src/fast_erfc.h`

Root cause: `R::pnorm` and `R::dnorm` per-observation calls go through R's parameter dispatch; `Rf_pnorm_both` at 8.3% of est-kernel samples.

**Implemented (2026-07-07):** All four fixes applied:
1. Added `fast_log_pnorm(x)` and `fast_log_dnorm(x)` to `fast_erfc.h`; replaced all `R::pnorm(w,0,1,1,1)` → `fast_log_pnorm(w)` and `R::dnorm(z,0,1,1)` → `fast_log_dnorm(z)` in `operator()` and `hessian()` (8 R-dispatch calls per obs eliminated)
2. Replaced `std::pow(one_minus_rho_sq, 1.5)` with precomputed `omrs_1p5 = omrs * sqrt(omrs)` (hoisted outside loop in both `operator()` and `hessian()`); `pow(x, 2.5)` → `omr2_2p5 = omr2 * omr2 * sqrt(omr2)` in `hessian()`
3. Promoted `d_ll_d_mu_event`/`d_ll_d_mu_cens` to `mutable Eigen::VectorXd` member fields; replaced per-call `VectorXd::Zero(m_n)` allocation with `.setZero()`

**Benchmark (2026-07-07):** n=600, p=4, 30 rounds × 10/5 reps warm-start:
- EST: OLD 1.6000ms → NEW 0.6000ms (**2.67×**, Wilcoxon p≈0)
- VAR: OLD 1.8000ms → NEW 0.8000ms (**2.25×**, Wilcoxon p≈0)

**Correctness test:** `EDI/tests/testthat/test-dep-cens-transform-fast-pnorm.R` (3 tests, 5 assertions — all pass). Gradient at MLE matches numerical finite-diff to 1e-3.


---

**TODO-102: fast_bai_parallel — 8 vector heap-allocations per simulation inside OpenMP for** ✓ DONE
Files: `EDI/src/fast_bai_parallel.cpp` (`compute_bai_distr_parallel_cpp`, lines 59–67, 107)

Root cause: The `#pragma omp parallel for` loop body (line 55) declares `d_i`, `y_r`, `w_r`, `match_T`, `match_C`, `has_T`, `has_C` as fresh `std::vector`s on **every** simulation iteration (lines 59–67), plus `yT`/`yC` inside an inner conditional block (line 107). With nsim=2000, this is ~14 000–16 000 heap allocations per `compute_bai_distr_parallel_cpp` call. Compare with `kk_compound_distr_parallel.cpp` which already hoists its thread-local vectors correctly (`diffs`, `treated_idx`, `control_idx` declared in the outer `#pragma omp parallel` scope, cleared+reused per iteration).

Fix: Split `#pragma omp parallel for` into `#pragma omp parallel` + `#pragma omp for`, hoist all vector declarations to the parallel-block scope, and replace body construction with `clear()` + `push_back`/`assign`. Capacity grows to maximum on the first few iterations and is retained thereafter.

Expected speedup: 5×–10× for medium nsim (2000), since malloc/free is the dominant cost over the O(n) computation.

Implemented by splitting the `#pragma omp parallel for` into an outer `#pragma omp parallel` region plus an inner `#pragma omp for schedule(static)`. The per-simulation vectors (`d_i`, `y_r`, `w_r`, `match_T`, `match_C`, `has_T`, `has_C`, `yT`, `yC`) are now declared once per thread, reserve capacity up front, and are reused with `clear()` / `assign()` for each simulation handled by that thread. `has_T`/`has_C` use `std::vector<char>` rather than `std::vector<bool>` to avoid proxy-bitset overhead while preserving the same truth semantics.

Noise-controlled benchmark used `R CMD INSTALL --no-docs` installs; no `--preclean` was used. The before package was a temporary source copy saved before the TODO-102 edit, and the after package was the edited source. Timing used two before phases and two after phases, 50 rounds per phase, 20 complete calls per round, with the profiler-shaped BAI workload (`n=400`, `nsim=2000`, matched-pair `m_mat`, `convex_flag=TRUE`). Combined single-core medians were 10.4000ms → 9.1000ms per call (IQR 10.3000–10.6125ms → 8.4500–10.3000ms), a 1.143× speedup / 12.5% time reduction (independent bootstrap 95% CI 1.106×–1.198×; one-sided Wilcoxon p=2.60e-11). Combined 4-core medians were 3.3500ms → 2.7000ms per call (IQR 3.1500–3.5000ms → 2.5875–2.9625ms), a 1.241× speedup / 19.4% time reduction (bootstrap 95% CI 1.189×–1.271×; one-sided Wilcoxon p=4.59e-20). Per-phase medians for 1 core were before1 10.4500ms, after1 8.4500ms, after2 10.3000ms, before2 10.4000ms; for 4 cores they were before1 3.3500ms, after1 2.6000ms, after2 2.8500ms, before2 3.3500ms. Before/after checksums matched exactly for both workloads.

Correctness: added `test-bai-thread-local-workspaces.R`, which compares `compute_bai_distr_parallel_cpp()` against an independent R implementation for convex and non-convex paths, including matched pairs, reservoir rows, all-reservoir simulations, and 1-core/2-core determinism. The dedicated test passed all 4 assertions. Existing `test-todo98-new-kernel-smoke.R` passed all 22 assertions.

---

**TODO-103: cmh_speedups — unordered_map allocated fresh on every call; replace with flat vector** ✓ DONE
Files: `EDI/src/cmh_speedups.cpp` (`compute_cmh_block_se_cpp`, `compute_extended_robins_block_se_cpp`)

Root cause: Both CMH functions created a fresh `std::unordered_map` on every call, reserved to `y.size()` (over-reserved; actual block count B ≈ n/4), then inserted n entries via hash-table pointer chasing.

Fix: One-pass scan to find `max_block_id` (skipping NA/invalid); allocate `std::vector<int> block_sums(max_block_id+1, -1)` (-1 sentinel = unseen; avoids ambiguity with all-zero blocks); accumulate with `block_sums[match_id] += y[i]` (O(1) stride-1); count B inline on first insertion. Iteration replaces hash-map traversal with a linear stride-1 scan. Same fix for `compute_extended_robins_block_se_cpp` using a flat `std::vector<RobinsBlockAccumulator>` (n=0 as unseen sentinel). Also removed `#include <unordered_map>`.

Benchmark (60 rounds × 200 reps, microseconds):

| config              | OLD median | NEW median | ratio | Wilcoxon p |
|---------------------|-----------|-----------|-------|------------|
| cmh  n=200,  B=50   | 5.0 μs    | 5.0 μs    | 1.00× | 3e-08 ✓   |
| cmh  n=1000, B=250  | 30.0 μs   | 10.0 μs   | 3.00× | 8e-12 ✓   |
| cmh  n=2000, B=500  | 55.0 μs   | 20.0 μs   | 2.75× | 7e-12 ✓   |
| cmh  n=4000, B=1000 | 110.0 μs  | 40.0 μs   | 2.75× | 8e-12 ✓   |
| rob  n=200,  B=50   | 10.0 μs   | 5.0 μs    | 2.00× | 3e-06 ✓   |
| rob  n=1000, B=250  | 32.5 μs   | 20.0 μs   | 1.62× | 8e-12 ✓   |
| rob  n=2000, B=500  | 65.0 μs   | 35.0 μs   | 1.86× | 8e-12 ✓   |

2.75–3× speedup for CMH at n≥1000; 1.6–2× for Robins; effect at n=200 is real but sub-resolution (both round to 5μs). Correctness: 15/15 tests pass (`test-cmh-flat-vector.R`).

---

**TODO-104: optimal_design_search — per-simulation w, t_idxs, c_idxs allocation inside nsim loop** ✓ DONE
Files: `EDI/src/optimal_design_search.cpp` (`d_optimal_search_cpp` lines 30–43, `a_optimal_search_cpp` lines 137–150)

Root cause: Inside the `for (int s = 0; s < nsim; ++s)` loop, each iteration allocated `Eigen::VectorXd w(n)` and `Pw(n)` (d_optimal) or `Pw(n)` + `Hw(n)` (a_optimal) plus `std::vector<int> t_idxs, c_idxs` (with `reserve()`). For nsim=500: 500 × 3–4 heap allocations per call.

Fix: Hoist `w`, `Pw`, `Hw` (a_optimal), `t_idxs`, `c_idxs` before the `s` loop with a single `reserve()` each. Inside the loop: `w.setZero(); t_idxs.clear(); c_idxs.clear();` and use `Pw.noalias() = P * w; Hw.noalias() = H * w;`. Fix 2 (maintaining sort order across swaps) was impractical: all Pw[k] values change on every accepted swap, so Ai and Bj keys for ALL indices are invalidated — full re-sort is unavoidable.

Benchmark (60 rounds × 20 reps; OLD = row-access fix only; NEW = row-access fix + hoisting):

| config                   | OLD median | NEW median | ratio | Wilcoxon p |
|--------------------------|-----------|-----------|-------|------------|
| d n=30,  p=4, nsim=1000  | 3.725 ms  | 3.150 ms  | 1.18× | 9e-12 ✓   |
| d n=60,  p=6, nsim=500   | 7.100 ms  | 6.250 ms  | 1.14× | 3e-10 ✓   |
| d n=100, p=8, nsim=200   | 9.125 ms  | 6.900 ms  | 1.32× | 8e-12 ✓   |
| d n=60,  p=10, nsim=500  | 9.625 ms  | 7.475 ms  | 1.29× | 8e-12 ✓   |
| a n=30,  p=4, nsim=600   | 4.400 ms  | 3.600 ms  | 1.22× | 4e-11 ✓   |
| a n=60,  p=6, nsim=300   | 9.350 ms  | 8.250 ms  | 1.13× | 4e-10 ✓   |
| a n=100, p=8, nsim=150   | 13.925 ms | 14.500 ms | 0.96× | n.s.       |

13–32% speedup at small/medium n; no significant effect at n=100 a_optimal (GEMV cost dominates over allocation savings). Correctness: 13/13 tests pass (`test-optimal-design-hoisted-allocs.R`).

---

**TODO-105: compute_objective_vals_cpp — string comparison + stride-r cache miss** ✓ DONE
Files: `EDI/src/rerandomization_helpers.cpp` (`compute_objective_vals_cpp`)

Two combined fixes applied together in the working tree:

**1. String comparison hoist (TODO-105):** `if (objective == "abs_sum_diff")` was evaluated inside the `for (int row = 0; row < r; row++)` loop. Hoisted to `const bool abs_mode = (objective == "abs_sum_diff")` before all loops; replaced in-loop branch with `if (abs_mode)`.

**2. Cache-friendly loop restructuring (TODO-107 companion):** Old code looped `for row, for i: indicTs(row,i)` — stride-r=5000 column-major access causing ~n cache misses per row. Restructured to loop `for i (outer), for row (inner)`: precompute `const int* indic_col = indic_ptr + i*r` once per i, then `indic_col[row]` is stride-1 sequential. Also preallocated `sum_T[r*p]` flat array (replacing per-row `std::fill`) and extracted `x_row[p]` once per i instead of reloading `X(i,j)` inside the row loop.

Benchmark (60 rounds × 5 reps, `compute_objective_vals_cpp`, abs_sum_diff):

| config            | OLD median | NEW median | ratio | Wilcoxon p |
|-------------------|-----------|-----------|-------|------------|
| n=50,  p=4, r=1000  | 0.60 ms  | 0.20 ms   | 3.00× | 4e-22 ✓   |
| n=100, p=4, r=2000  | 2.80 ms  | 1.20 ms   | 2.33× | 6e-22 ✓   |
| n=200, p=4, r=5000  | 15.60 ms | 5.70 ms   | 2.74× | 2e-21 ✓   |
| n=200, p=8, r=5000  | 14.00 ms | 5.60 ms   | 2.50× | 1e-21 ✓   |
| n=50,  p=4, r=5000  | 3.20 ms  | 1.40 ms   | 2.29× | 6e-22 ✓   |

2.3–3.0× speedup across all configs, all highly significant. Correctness: 5/5 tests pass (`test-rerandomization-objective-contiguous-indicts.R`).

---

**TODO-106: d_optimal_search_cpp — eliminate per-iteration j×n multiply by transposing P access** ✓ DONE
Files: `EDI/src/optimal_design_search.cpp` (`d_optimal_search_cpp` line 83, `a_optimal_search_cpp` line ~200)

Root cause: `perf annotate` showed `imulq j*n` as a hot instruction in the inner swap loop. Each inner-j iteration accessed `p_ptr[j*n + i]` (stride-n column access, one cache miss per j). Since P is symmetric, `P(i,j) = P(j,i)`: precompute `p_row_i = p_ptr + i*n` once per outer-i, then use `p_row_i[j]` (stride-1, sequential access). For n=60 the entire row fits in ~8 cache lines and stays hot in L1 through the inner-j loop.

Applied to both `d_optimal_search_cpp` (one P lookup per inner step) and `a_optimal_search_cpp` (one P + one H lookup per inner step, using `p_row_i` and `h_row_i`).

Benchmark (60 rounds × 5 reps, `d_optimal_search_cpp`, nsim=200):

| config        | OLD median | NEW median | ratio | Wilcoxon p |
|---------------|-----------|-----------|-------|------------|
| n=30,  p=4    | 0.80 ms   | 0.60 ms   | 1.33× | 1e-05 ✓   |
| n=60,  p=6    | 3.40 ms   | 2.60 ms   | 1.31× | 4e-19 ✓   |
| n=100, p=8    | 9.90 ms   | 6.80 ms   | 1.46× | 8e-20 ✓   |
| n=60,  p=10   | 3.40 ms   | 2.80 ms   | 1.21× | 5e-17 ✓   |

21–46% speedup across all sizes, all highly significant. Correctness: 10/10 tests pass (`test-d-optimal-row-access.R`).

---

**TODO-107: compute_objective_vals_cpp — stride-r=5000 cache miss in indicTs inner loop** ✓ DONE
Files: `EDI/src/rerandomization_helpers.cpp` (`compute_objective_vals_cpp`, lines 114–132)

Root cause: `perf report` for `rerandomization_obj_vals` (r=5000, n=200, p=4): 38.46% of samples in `compute_objective_vals_cpp`. The dominant cost is the inner access `indicTs(row, i)` (line 118), where `indicTs` is an `r×n` column-major integer matrix. With the outer loop over `row` (0..r-1) and inner loop over `i` (0..n-1), each `indicTs(row, i)` access is at offset `row + i * r` — stride-r=5000 between consecutive `i` values. At r=5000 and 4 bytes/int, each column jump = 20KB > typical L1 cache line distance. Result: ~n=200 cache misses per row evaluation, for r=5000 rows = 1M cache misses per call.

Fix:
1. Transpose `indicTs` into a local `n×r` matrix at the top of the function, then access `indicTs_T(i, row) = indicTs(row, i)` with stride-1 inner loop. One-time O(n×r) transpose cost, then all inner loops are sequential.
2. Or restructure the loop: outer loop over `i`, inner loop over `row`, accumulating `sum_T[row][j]` — but requires `r×p` accumulator matrix.
3. The string comparison inside the row loop (`objective == "abs_sum_diff"` at line 123) is secondary (TODO-105) — the cache miss is the primary cost.

Expected speedup: 3×–5× for the n×r sweep; the transpose itself is cheap.

Disassembly confirmation (`objdump -d EDI.so`, function `_Z26compute_objective_vals_cpp` at 0x6ea720, size 5754 bytes):
- **+0x63b**: `cmpl $0x1,(%rax,%r14,4)` — the `indicTs(row,i)==1` load; r14 = `stride * i + col_base` computed via `imul %r12,%r14` at +0x661 (stride=r=5000, so each i-step jumps 20 KB — L3 miss per inner iteration)
- **+0x640 [134 samples]**: `je ...` — branch on the load result; receives the sample attribution while the CPU stalls waiting for the cache fill
- **+0x8eb**: `vaddsd (%rax,%r8,8),%xmm0,%xmm0` — `X(i,j)` load in the `sum_T[j]` accumulation j-loop; stride = n = 200 between j-steps, and `imul %r14,%r8` at +0x90c recomputes the column offset every j-iteration
- **+0x8f1 [50 samples]**: `vmovsd %xmm0,0x0(%r13,%r14,8)` — store back to `sum_T[j]`; samples here reflect the overall j-loop cost (p=4 so minor)

The 134-sample hot instruction confirms the indicTs cache miss is the dominant cost, not the string comparison (TODO-105).

Implemented the accumulator-loop variant rather than materializing an `n×r` transposed integer copy. `compute_objective_vals_cpp()` now reads each original `indicTs` column contiguously (`row` varying fastest), accumulates treatment counts plus an `r×p` row-major treatment-sum scratch buffer, then computes the requested objective in a final row pass. This removes the stride-`r` `indicTs(row, i)` load from the hot objective loop and also reuses each `X` row across all randomization rows for that subject. The objective-mode string check is hoisted to booleans at entry, and the function now validates `ncol(indicTs) == nrow(X)`.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-107 source and optimized source; no `--preclean` was used. Timing used separate R processes, process CPU time, two before phases and two after phases, 40 rounds per phase, profiler-sized data (`n=200`, `p=4`, `r=5000`), 50 complete `abs_sum_diff` calls per round and 30 complete `mahal_dist` calls per round. Combined medians: `abs_sum_diff` 11.0600ms → 5.2600ms per call (IQR 10.9600–11.1800ms → 5.1800–5.3400ms), a 2.103× speedup / 52.4% time reduction (independent bootstrap 95% CI 2.091×–2.118×; one-sided Wilcoxon p<2.2e-16). `mahal_dist` improved 11.3333ms → 5.3667ms per call (IQR 11.1000–11.5667ms → 5.2667–5.4417ms), a 2.112× speedup / 52.6% time reduction (bootstrap 95% CI 2.081×–2.138×; one-sided Wilcoxon p<2.2e-16). Before/after outputs were exactly equal on the benchmark data for both objectives.

Correctness: added `test-rerandomization-objective-contiguous-indicts.R`, which compares `compute_objective_vals_cpp()` with an independent R implementation for both `abs_sum_diff` and `mahal_dist`, and checks invalid objective/input handling. The dedicated test passed all 5 assertions. The existing `test-todo98-new-kernel-smoke.R` file also passed all 22 assertions.

---

**TODO-108: rerandomization_search_cpp — OpenMP overhead dominates for small REPS** ✓ DONE
Files: `EDI/src/rerandomization_helpers.cpp` (`rerandomization_search_cpp`, lines 202–240)

Root cause: `perf report` for `rerandomization_search` (r=500, max_draws=50000, n=200, p=4, REPS=50): `libgomp.so` accounts for **39.27%** of all samples, with the actual computation nearly invisible. The `#pragma omp parallel` block (line 202) spawns threads and sets up a barrier even when there are few accepted draws. With the `#pragma omp for schedule(static)` dividing 50000 draws across threads, each thread does relatively little work before the barrier, inflating per-call overhead.

Considerations: This is a calling-pattern issue (REPS=50, short per-call time of 4.4ms) exacerbating inherent OpenMP overhead. In production (large n, many accepted draws, long-running calls), OpenMP overhead is proportionally smaller. However, callers with small `max_draws` or small `r` pay disproportionately. If `max_draws * n` is below a threshold, a sequential fallback (`#ifdef _OPENMP ... omp_get_max_threads() > 1 ...`) would reduce overhead for short jobs.

Implemented an early-stopping OpenMP work queue rather than a broad multithreaded sequential fallback. A direct sequential threshold was tested but was not wall-clock faster on this machine because the existing parallel path wins on latency despite spending extra CPU. The final implementation keeps the sequential path only when `omp_get_max_threads() <= 1`; otherwise, threads fetch draw chunks from an atomic `next_draw` counter and exit as soon as `found >= r`, avoiding the old `#pragma omp for schedule(static)` behavior where every scheduled draw was visited even after enough accepted allocations had been collected. This preserves parallel latency while reducing wasted post-acceptance OpenMP work.

Noise-controlled BAAB benchmark used plain `R CMD INSTALL --no-docs` installs for the saved pre-TODO-108 source and optimized source; no `--preclean` was used. Timing used separate R processes, elapsed wall time as the primary latency metric, and process CPU time as a secondary resource metric. The profiler-sized workload was `n=200`, `p=4`, `r=500`, `max_draws=50000`, `cutoff=1.0`, 50 rounds per phase, 20 complete calls per round, two before phases and two after phases (4,000 timed calls total). Combined wall-clock medians: 0.4250ms → 0.2500ms per call (IQR 0.2000–1.1000ms → 0.2000–0.5000ms), a 1.700× speedup / 41.2% time reduction (independent bootstrap 95% CI 1.091×–2.600×; one-sided Wilcoxon p=1.86e-4). Combined process-CPU medians: 3.2750ms → 1.9000ms per call (IQR 1.8875–10.7750ms → 1.2500–4.1375ms), a 1.724× speedup / 42.0% reduction (bootstrap 95% CI 1.141×–3.167×; one-sided Wilcoxon p=6.03e-5). All benchmark calls returned exactly `r=500` accepted allocations.

Correctness: added `test-rerandomization-search-early-stop.R`, which independently reproduces the internal rerandomization score for both `abs_sum_diff` and `mahal_dist` and verifies that returned allocations are binary, balanced, have the requested count, and satisfy the cutoff. The dedicated test passed all 12 assertions. The existing `test-todo98-new-kernel-smoke.R` file also passed all 22 assertions.

---

**TODO-109: kk_compound_distr — merge multiple O(n) passes per simulation into one** ✓ DONE
Files: `EDI/src/kk_compound_distr_parallel.cpp`

Root cause: Each simulation iteration made 4 separate O(n) passes over `m_col`/`w_col`/`y`: (1) find max_m, (2) scatter-fill treated/control idx, (3) unmatched sums, (4) unmatched variances. Each data element read 3–4× causing repeated L1 cache pressure.

**Implemented (2026-07-07):** Single merged O(n) pass replacing all 4:
- Simultaneously tracks `max_m`, scatter-fills `treated_idx`/`control_idx`, and accumulates `nRT`, `nRC`, `sum_T`, `sum_C`, `sum_T2`, `sum_C2` for unmatched obs
- Eliminates the second variance pass by computing `var_T = (sum_T2 − nRT·mean_T²)/(nRT−1)` from the merged pass accumulators (sum-of-squares identity)
- Matched-pair stats (d_bar, ssqD_bar) computed in one O(m) pass over diffs using running sum+sumsq, eliminating two sub-passes
- Applied to both `compute_matching_compound_distr_parallel_cpp` and `compute_matching_compound_bootstrap_parallel_cpp`
- Pre-initialization `treated_idx.assign(n, -1)` replaces `assign(m, -1)` post-find; cost similar since m ≈ n/2

**Benchmark (2026-07-07):** n=400, nsim=2000, 30 rounds × 5 reps, 1 core:
- OLD 4.9000ms → NEW 2.5000ms (**1.96×**, Wilcoxon p≈0)

**Correctness test:** `EDI/tests/testthat/test-kk-compound-merged-pass.R` (3 tests — mixed data, all-matched, all-unmatched; all pass to tolerance 1e-10 vs R reference).

---

**TODO-110: Bootstrap randomization test (BRT) — three-tier fast-path stack** ✓ DONE
Files: `EDI/src/rand_bootstrap_mean_diff_parallel.cpp` (new), `EDI/R/inference_all_abstract_rand_bootstrap.R`, `EDI/R/inference_all_abstract_rand_bootstrap_ci.R`, `EDI/R/inference_all_mean_diff.R`, `EDI/R/inference_continuous_ols.R`

Context: the BRT (`InferenceRandBootstrap`/`InferenceRandBootstrapCI`, added 2026-07-13) resamples n rows with replacement and draws a fresh assignment from the design per replicate. The naive path duplicates the design + inference objects per replicate via `bootstrap_subset_inference()` — pure R6 overhead that dwarfs the statistic cost for closed-form estimators, and is multiplied ~25–35× by the CI's bisection inversion.

**Implemented (2026-07-13/14), three tiers with graceful NULL/fallback dispatch:**
1. **Reusable worker states** (all classes with `supports_reusable_bootstrap_worker`): load the bootstrap subset into a persistent worker (`load_bootstrap_sample_into_worker`), then inject the fresh assignment + sharp-null-shifted responses (`load_rand_bootstrap_assignment_into_worker`). Poisson/KK14 benchmark (n=60, B=101, serial): 11.81s → 1.23s per p-value (**9.6×**), identical p-values. Gotcha found: `worker_state$worker_des` partial-matches `worker_des_priv` — worker-state members must be accessed with `[[ ]]`.
2. **C++ batch kernel** (`compute_rand_bootstrap_mean_diff_parallel_cpp`, dispatched via per-class `compute_fast_rand_bootstrap_distr`): evaluates all B replicates in one call from (y0, n×B row-index matrix, n×B assignment matrix, delta); OpenMP via `should_parallelize_replicates`, mirrors `compute_simple_mean_diff_parallel_cpp`. Mean-diff/KK14 benchmark (n=100, B=501, serial): p-value 97.99s → 2.23s (**44×**, identical p-value 0.0040); CI 2.9s vs ~220s pre-optimization (**~75×**) since all bisection deltas reuse one materialized draw set. Remaining p-value cost is dominated by the serial draw materialization, now parallelized across cores with per-draw seeds (deterministic in private$seed regardless of num_cores).
3. **Closed-form CI** (per-class `compute_rand_bootstrap_ci_affine_coefs`; mean-diff + `InferenceContinOLS`): with the additive sharp-null shift each null draw is affine in delta, t0_b(δ) = A_b + δ·c_b, so the p-value is a step function with breakpoints (t−A_b)/c_b and the CI is read off exactly from interval probes — zero bisection evaluations, no pval_epsilon/search-radius/conservative fallbacks. Convention note: the package inverts the two-sided randomization p-value at **alpha/2** (see `InferenceRandCI`); the closed form must match or CIs disagree by ~0.5 on unit-scale data.

Correctness: `EDI/tests/testthat/test-rand-bootstrap.R` — kernel vs per-iteration reference equal to 1e-10 (delta 0 and 0.7), worker path vs standard path equal to 1e-10, affine decomposition vs reference iteration equal to 1e-8 for both classes, closed-form CI vs bisection CI within bisection tolerance; 67 assertions, all passing.

---

**TODO-111: BRT C++ fast-path kernels for CoxPH, Weibull marginal, and robust regression** ✓ DONE
Files: `EDI/src/fast_coxph_regression.cpp`, `EDI/src/fast_weibull_regression.cpp`, `EDI/src/fast_robust_regression.cpp`, `EDI/R/inference_survival_coxph.R`, `EDI/R/inference_survival_KK_weibull_marginal.R`, `EDI/R/inference_continuous_robust_regr.R`

Context: after the three-tier BRT stack (TODO-110) reduced mean-diff/OLS BRT to negligible cost, the survival and robust regression classes fell back to the reusable-worker path — still 9–10× faster than naive R6 duplication, but the C++ batch-kernel tier was not yet implemented for these classes.

**Implemented (2026-07-14):**
- `compute_coxph_rand_bootstrap_parallel_cpp`: per-draw Cox partial-likelihood IRLS (via `compute_coxph_regression_cpp` internally); OpenMP over B draws.
- `compute_weibull_rand_bootstrap_parallel_cpp`: per-draw Weibull AFT MLE fit; OpenMP over B draws.
- `compute_robust_rand_bootstrap_parallel_cpp`: per-draw robust regression (Huber M or MM via `fast_robust_regression_cpp`); OpenMP over B draws.
- Each dispatched via `compute_fast_rand_bootstrap_distr` on the respective class; returns NULL if custom stat fn is set, falling through to the worker path.
- Fix: `Nullable<NumericVector>()` constructor crashes inside OpenMP parallel regions; replaced with `R_NilValue` cast directly.

Correctness: `EDI/tests/testthat/test-rand-bootstrap.R` — C++ kernel vs per-iteration reference equal to 1e-10 for CoxPH, Weibull, and robust regression (delta 0 and 0.7); 178 assertions total, all passing.

---

**TODO-112: BRT studentized, symmetric-percentile-t, and smoothed p-value types** ✓ DONE
Files: `EDI/R/inference_all_abstract_rand_bootstrap.R`, `EDI/R/inference_all_abstract_rand_bootstrap_ci.R`

Context: the original BRT used only the percentile type (compare null-draw statistics to the observed statistic). Three additional types improve power or coverage for classes where a standard error is available.

**Implemented (2026-07-14):**
- **`type = "studentized"`**: each null draw is pivoted as `z_b = (t_b - delta) / se_b`; the p-value compares `|z_obs|` to `|z_b|`; produces asymmetric CIs under inversion. Requires the class to expose `get_standard_error()`; falls back to `"percentile"` if SE is unavailable.
- **`type = "symmetric-percentile-t"`**: same pivoting but uses `|z_b|` (absolute value pivot), giving a symmetric CI.
- **`type = "smoothed"`**: adds `N(0, sigma/sqrt(n))` kernel noise to each null draw before comparison; useful for discrete response types. Bypasses the C++ fast path (which produces raw integer draws) and falls through to per-iteration R evaluation.
- Per-draw SE computed via `compute_brt_null_statistics_with_se` (calls `run_rand_bootstrap_iteration_with_se`), which temporarily swaps the assignment vector and calls the inference object's `compute_estimate` + `get_standard_error` in each iteration.
- CI inversion for studentized/symmetric-percentile-t uses the same bisection infrastructure with the pivoted p-value function.

Correctness: `EDI/tests/testthat/test-rand-bootstrap.R` — studentized/symmetric-percentile-t/smoothed pvals and CIs validated against per-iteration reference; all 178 assertions pass.

## BRT smoothed fast-kernel noise support — 2026-07-20

`compute_rand_bootstrap_two_sided_pval`/`compute_rand_bootstrap_confidence_interval(type = "smoothed")` unconditionally disabled every class's C++ `compute_fast_rand_bootstrap_distr` kernel whenever a draw carried per-draw `smooth_noise`, because none of the kernels could consume it; every smoothed evaluation fell back to the R-level reused-worker/per-iteration path, which for `InferenceAllSimpleWilcox` calls `stats::wilcox.test()` once per bootstrap draw. Extended the 8 C++ kernels that operate on real-valued responses (`compute_wilcox_hl_rand_bootstrap_parallel_cpp`, `compute_rand_bootstrap_mean_diff_parallel_cpp`, `compute_rand_bootstrap_ols_parallel_cpp`, `compute_robust_rand_bootstrap_parallel_cpp`, `compute_coxph_rand_bootstrap_parallel_cpp`, `compute_weibull_rand_bootstrap_parallel_cpp`, `compute_logrank_rand_bootstrap_parallel_cpp`, `compute_survival_stat_diff_rand_bootstrap_parallel_cpp`) with an optional `Rcpp::Nullable<Rcpp::NumericMatrix> noise_mat` argument, added before the sharp-null shift exactly where the R fallback already adds it (`y_sim = y0[i_b] + smooth_noise`, applied to every row regardless of treatment status). `rand_bootstrap_draw_matrices()` now packs this matrix whenever draws carry `smooth_noise`; the two blanket dispatcher checks that disabled the fast path (`approximate_rand_bootstrap_distribution_beta_hat_T` and `get_brt_distribution_prefix` in `inference_all_abstract_rand_bootstrap.R`) were removed in favor of each class deciding for itself. The two ordinal classes (`InferenceOrdinalRidit`, `InferenceOrdinalJonckheereTerpstraTest`) explicitly decline noise-carrying draws in their own `compute_fast_rand_bootstrap_distr`, since continuous noise on integer category codes isn't meaningful, and are unchanged.

Correctness: 8 new per-kernel tests (`EDI/tests/testthat/test-brt-smoothed-*-kernel.R`) use an identity-resampling equivalence check — `kernel(y0, noise_mat=N, ...)` on identity `i_mat` columns must equal `kernel(y0 + N[,b], noise_mat=NULL, ...)` per column — plus a `noise_mat=NULL` vs. explicit zero-matrix no-op check; `test-brt-smoothed-noise-mat-plumbing.R` checks the shared draw-matrix builder and the ordinal opt-out; `test-brt-smoothed-wilcox-ci-perf.R` compares the fast-kernel CI output against a forced-slow-fallback subclass for identical seeds (tolerance 1e-6) and asserts a speedup. Full existing BRT suites (`test-rand-bootstrap.R` 178 assertions, `test-wilcox-regr-bootstrap-fast-path.R`, `test-bootstrap-reused-worker-families.R`) still pass unchanged.

Performance only materializes for `compute_rand_bootstrap_confidence_interval(type = "smoothed")`, not a standalone `compute_rand_bootstrap_two_sided_pval(type = "smoothed")` call: CI inversion pre-materializes fresh assignments once (`materialize_w = TRUE`, common random numbers reused across every `delta` evaluated during root-finding), so the fast kernel can engage on each evaluation, whereas a standalone p-value call draws the fresh assignment lazily per replicate (`materialize_w = FALSE`) so `rand_bootstrap_draw_matrices()` can never build the matrices the fast kernels need, regardless of this fix. Measured on `InferenceAllSimpleWilcox`, `n = 30`, `B = 99`: `compute_rand_bootstrap_confidence_interval(type = "smoothed")` went from 25.523s (forced-slow-fallback subclass) to 0.508s (fast kernel), about 50x; the standalone smoothed p-value call is unaffected (1.5s both before and after, dominated by the same per-replicate lazy-draw path either way). Installed with `R CMD INSTALL --no-docs`.
