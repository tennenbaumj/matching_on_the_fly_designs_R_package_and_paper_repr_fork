# Path Audit Hardening Report

Generated: 2026-07-20

## Summary

The current `path_audits.html` intentionally uses two successful-run colors:

- Dark green: the method is expected to return numeric output under the comprehensive-test contract.
- Light green: the method must be attempted, but may legitimately return explicit non-estimable output because the estimate, standard error, resampling distribution, or interval inversion target is undefined or numerically degenerate.

`status=error` is different. It indicates a harness or package failure and is not acceptable in a bug-free package/test harness. Those rows should continue to be monitored and fixed.

The current rendered audit is very conservative:

```text
dark green:   22
light green: 2551
SLOW:        497
N/A:         890
```

This does not mean 2551 paths have demonstrated failures. It means the renderer currently marks whole method families as potentially non-estimable unless we have a per-class/per-method guarantee.

## Observed Non-Estimable Mechanisms

Current result CSVs show recognized non-estimable output concentrated in these mechanisms:

```text
bootstrap_standard_error_ci_unavailable
bootstrap_too_few_finite_standard_errors
bootstrap_bca_adjustment_on_boundary
bootstrap_bca_jackknife_unavailable
bayesian_bootstrap_* standard-error / BCa boundary failures
jackknife_nonfinite_replicate_estimates
randomization_too_few_finite_estimates
score_test_unavailable
rand_bootstrap_ci_estimate_unavailable
```

These are mostly real finite-sample or numerical limitations, not package errors, provided they are explicitly typed and recorded as non-estimable.

## What Can Become Dark Green

### Good Candidates

Bootstrap percentile/basic/symmetric paths are the strongest candidates for promotion. They mainly require a finite bootstrap estimate distribution, not finite per-replicate standard errors or jackknife accelerations.

Likelihood-based model paths may also be promotable class by class after stronger fallback fitting, rank reduction, separation detection, and bracketing improvements.

Simple closed-form model-based Wald paths can remain dark green when they have a finite-sample numeric guarantee under the comprehensive-test contract.

## Actionable To-Dos

### Audit Refinement

- [x] Add per-method audit metadata instead of coloring entire families at once.
  Suggested fields:
  - `always_numeric_methods`
  - `maybe_nonestimable_methods`
  - `slow_methods`
  - `unsupported_methods`

- [ ] Promote cells to dark green only when both conditions hold:
  - the method has a defensible numeric guarantee or statistically valid fallback
  - comprehensive tests show no recognized non-estimable outputs for that method/class/design coverage

- [ ] Generate an audit reconciliation table from result CSVs:
  - expected attempted paths from `path_audits_source.R`
  - actually called paths from comprehensive-test logs/CSVs
  - paths returning numeric output
  - paths returning explicit non-estimable output
  - paths returning `status=error`

- [ ] Add a third internal audit state, even if not displayed separately yet:
  - `numeric_observed`
  - `nonestimable_observed`
  - `error_observed`

- [ ] Split bootstrap family coloring by method:
  - percentile/basic/symmetric p-value and CI candidates can be darker after validation
  - studentized and BCa should remain light unless class-specific guarantees exist

### Package Hardening

- [ ] Add adaptive bootstrap estimate collection for percentile/basic/symmetric methods:
  - continue drawing until a minimum finite count is reached
  - cap attempts to avoid infinite loops
  - record typed non-estimable output if the finite-count target cannot be met

- [ ] Add preflight checks before bootstrap/refit loops:
  - both treatment arms present
  - enough events/categories/nonzero responses
  - no structurally empty model component
  - rank is sufficient after design-matrix hardening

- [ ] Add class-specific robust bootstrap refit fallbacks:
  - treatment-only fallback where covariates are singular
  - QR column dropping
  - separation-aware fallback for binary/ordinal models
  - event/category preserving bootstrap where statistically appropriate

- [ ] Harden likelihood score/LR/gradient paths:
  - improve cold starts
  - add ridge/penalized fallback where appropriate
  - detect separation and boundary estimates early
  - return typed non-estimable output instead of untyped errors

- [ ] Harden randomization paths:
  - pre-check observed statistic availability
  - adaptively draw more assignments when too few finite statistics are available
  - preserve required strata/block/event/category structure where valid
  - return typed non-estimable output when valid assignment draws cannot identify the statistic

- [ ] Harden jackknife paths:
  - prefer delete-block/delete-cluster where delete-one breaks design structure
  - detect when a deletion creates an empty arm, no events, no categories, or singular fit
  - keep typed non-estimable output for irreducible cases

### Harness Invariants

- `status=error` is a failing condition requiring investigation.

- Explicit non-estimable output is acceptable only when:
  - the package sets `is_nonestimable()`
  - the reason and stage are typed
  - the harness recognizer permits that method family

- Generic `NA`, `NaN`, or `Inf` should not be accepted unless paired with explicit non-estimable state or a narrowly documented exception.

- Todo: add summary counts to the comprehensive-test output:
  - numeric ok
  - explicit non-estimable ok
  - skipped/SLOW
  - unsupported/N/A
  - error

## Suggested Priority Order

1. Refine the audit metadata so light green is not family-wide by default.
2. Add result-driven reconciliation from comprehensive-test CSVs.
3. Harden percentile/basic bootstrap paths first.
4. Harden likelihood score/LR/gradient paths class by class.
5. Keep studentized, BCa, jackknife, randomization, and BRT light green unless a class-specific guarantee is established.


# cd /home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr && Rscript package_tests/comprehensive_tests.R 40 1 proportion  2>&1 | tee package_tests/comprehensive_tests_proportion_$(date +%Y%m%d_%H%M%S).log