# EDI inference class path audit
# Run html_from_audit(audit_classes, "path_audits.html") to regenerate the HTML table.
# testing_types: "full"=c(wald,score,grad,lr); "wald"=wald only; "wald_lr"; "wald_score_grad"; "wald_score_lr"; "none"=character(0)
# rand/rci applicable by response type: rand→{cont,surv,prop,incid,count,ord}; rci→{cont,prop,surv} (count/ord always skip_ci_rand; incid no rci)
# pboot: TRUE=supports parametric bootstrap (is(InferenceParamBootstrap) && supports_lik_ratio_param_bootstrap()=TRUE);
#        FALSE=is ParamBootstrap but returns FALSE; NA=not InferenceParamBootstrap at all
# Bayesian bootstrap has 6 distinct flavors (confirmed against comprehensive_tests_*.log "Calling ...()" labels
# across all response types), each an explicit "compute_bayesian_bootstrap_{two_sided_pval,confidence_interval}[_type]"
# path: percentile (pval+ci, default/untyped call), symmetric (pval only), basic (ci only), wald (pval+ci),
# bca (pval+ci), studentized (pval+ci, aka "bootstrap-t"). skip_bbt/skip_bbt_ci still apply uniformly to all 6 --
# no per-flavor skip data is tracked yet.
# Optional per-row method overrides: always_numeric_methods, maybe_nonestimable_methods, slow_methods,
# unsupported_methods. Values should match the comprehensive_tests.R function_run labels used below.

audit_classes = list(

  # ── GLOBAL ──────────────────────────────────────────────────────────────────
  list(name="InferenceAllSimpleMeanDiff",           section="Global",     resp="all",   kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, slow_methods=c("compute_bootstrap_two_sided_pval_studentized"), notes="rci for cont only (non-cont AllSimpleMeanDiff skipped); bootstrap studentized pval avg 178.8s / max 2001.5s at n=12 slow"),
  list(name="InferenceAllSimpleMeanDiffPooledVar",  section="Global",     resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, notes="inherits AllSimpleMeanDiff; tested cont only"),
  list(name="InferenceAllSimpleWilcox",             section="Global",     resp="all",   kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=NA,    jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, notes="boot enabled; Bayesian bootstrap structurally unsupported (fractional weights undefined for rank statistic); BRT CI smoothed avg 94s slow"),

  # ── CONTINUOUS KK ────────────────────────────────────────────────────────────
  list(name="InferenceContinKKGLMM",               section="Continuous", resp="cont",  kk=TRUE,  types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_jack_slow=TRUE, skip_rand=TRUE,  rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp",  pboot=TRUE,  notes="InferenceParamBootstrap; all resampling too slow; jackknife estimate avg 54s / max 102s (skip_jack_slow)"),
  list(name="InferenceContinKKOLSOneLik",          section="Continuous", resp="cont",  kk=TRUE,  types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  bartlett_exact=TRUE, notes="KKPassThroughCompound → ParamBootstrap; explicit pboot=TRUE; exact Bartlett implemented (holds sigma2 fixed, LR reduces exactly to the classical F(1,n-p) statistic; verified against lm())"),
  list(name="InferenceContinKKRobustRegrOneLik",   section="Continuous", resp="cont",  kk=TRUE,  types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, slow_methods=c("compute_rand_confidence_interval(custom)"), notes="KKPassThroughCompoundNoParamBootstrap; wald only; custom rand CI avg 336.6s / max 1994.8s at n=6 slow"),
  list(name="InferenceContinKKQuantileRegrOneLik", section="Continuous", resp="cont",  kk=TRUE,  types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA,    slow_methods=c("compute_bootstrap_confidence_interval"), notes="AbstractQuantileRandCI → KKPassThroughCompoundNoParamBootstrap; bootstrap CI mean 109.0s / max 9434.3s at n=98 slow"),

  # ── CONTINUOUS non-KK ────────────────────────────────────────────────────────
  list(name="InferenceContinLin",                  section="Continuous", resp="cont",  kk=FALSE, types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceContinOLS",                  section="Continuous", resp="cont",  kk=FALSE, types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceContinRobustRegr",           section="Continuous", resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp",  pboot=NA,    skip_brt="ci", notes="InferenceAsymp"),
  list(name="InferenceContinQuantileRegr",         section="Continuous", resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA,    notes="InferenceAsymp"),

  # ── INCIDENCE KK ─────────────────────────────────────────────────────────────
  list(name="InferenceIncidKKCondLogitOneLik",          section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap directly; explicit TRUE"),
  list(name="InferenceIncidKKCondLogitPlusGLMMOneLik",  section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, skip_bbt_ci=TRUE, slow_methods=c("compute_rand_two_sided_pval(delta=0.5)", "compute_bayesian_bootstrap_two_sided_pval", "compute_bayesian_bootstrap_two_sided_pval_symmetric"), notes="AbstractKKCondLogitPlusGLMM → InferenceParamBootstrap; supports_lik_ratio_param_bootstrap=FALSE (no simulate_under_lik_null); rand delta=0.5 pval avg 205s; Bayesian bootstrap pval avg 32.4s / max 68.4s at n=32 slow; symmetric Bayesian bootstrap pval avg 30.6s / max 54.9s at n=18 slow"),
  list(name="InferenceIncidKKGEE",                     section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidKKNewcombeRiskDiff",         section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="KKPassThroughCompoundNoParamBootstrap"),
  list(name="InferenceIncidKKGCompRiskDiff",            section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, notes="AbstractKKMarginalIncid → ParamBootstrap; supports_lik_ratio_param_bootstrap()=FALSE (no simulate_under_lik_null)"),
  list(name="InferenceIncidKKGCompRiskRatio",           section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, slow_methods=c("compute_bootstrap_two_sided_pval_symmetric"), notes="AbstractKKMarginalIncid → ParamBootstrap; supports_lik_ratio_param_bootstrap()=FALSE (no simulate_under_lik_null); bootstrap symmetric pval max 10456.5s slow"),
  list(name="InferenceIncidKKModifiedPoisson",          section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AbstractKKModifiedPoisson → AbstractKKMarginalIncid → ParamBootstrap; simulate_under_lik_null: Poisson draw"),

  # ── INCIDENCE non-KK ──────────────────────────────────────────────────────────
  list(name="InferenceIncidExactFisher",               section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact → InferenceJackknife; runs exact pval+CI only"),
  list(name="InferenceIncidExactBinomial",             section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact; KK14/FixedBinaryMatch"),
  list(name="InferenceIncidenceExactZhang",            section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="",  skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact; exact Zhang disables inherited bootstrap/Bayesian bootstrap because generic resampling weights do not define a valid exact Zhang estimator"),
  list(name="InferenceIncidWald",                     section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="inherits AllSimpleMeanDiff (ParamBootstrap); explicit FALSE"),
  list(name="InferenceIncidLogRegr",                  section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache; approximate Bartlett-corrected LR (pval+CI) implemented via generic InferenceParamBootstrap Monte-Carlo factor"),
  list(name="InferenceIncidProbitRegr",               section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferenceIncidMiettinenNurminenRiskDiff",section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidNewcombeRiskDiff",         section="Incidence", resp="incid", kk=FALSE, types="none",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp; get_supported_testing_types_impl=character(0); no direct tests"),
  list(name="InferenceIncidRiskDiff",                 section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    slow_methods=c("compute_bootstrap_confidence_interval", "compute_bootstrap_confidence_interval_studentized"), notes="AsympLikStdModCacheNoParamBootstrap; supports_likelihood_tests=FALSE; bootstrap CI/studentized CI avg 261.1s / max 2014.1s at n=8 slow"),
  list(name="InferenceIncidGCompRiskDiff",            section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidGCompRiskRatio",           section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidModifiedPoisson",          section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="ParamBootstrap; explicit FALSE"),
  list(name="InferenceIncidLogBinomial",              section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferenceIncidBinomialIdentityRiskDiff", section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),

  # ── PROPORTION ───────────────────────────────────────────────────────────────
  list(name="InferencePropKKGEE",                   section="Proportion", resp="prop",  kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=NA,    slow_methods=c("compute_rand_confidence_interval", "compute_rand_bootstrap_confidence_interval", "compute_rand_bootstrap_confidence_interval_studentized", "compute_rand_bootstrap_confidence_interval_symmetric-percentile-t", "compute_rand_bootstrap_confidence_interval_smoothed"), notes="InferenceAsymp; rand CI avg 30.1s at n=18; BRT CI avg 36.9-42.6s at n=18 slow"),
  list(name="InferencePropKKGLMM",                  section="Proportion", resp="prop",  kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=TRUE,  notes="AbstractKKCondLogitPlusGLMM → InferenceParamBootstrap; simulate_under_lik_null added"),
  list(name="InferencePropKKQuantileRegrOneLik",    section="Proportion", resp="prop",  kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=NA,    notes="AbstractKKQuantileRegrOneLik → KKPassThroughCompoundNoParamBootstrap; BRT CI studentized avg 34s + smoothed avg 33s slow"),
  list(name="InferencePropQuantileRegr",            section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp", pboot=NA,    slow_methods=c("compute_rand_confidence_interval"), notes="InferenceAsymp; logit-scale quantile regression, mirrors InferenceContinQuantileRegr; rand CI inversion is too slow/pathological with covariate adjustment"),
  list(name="InferencePropBetaRegr",                section="Proportion", resp="prop",  kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=TRUE,  slow_methods=c("compute_rand_confidence_interval", "compute_rand_bootstrap_confidence_interval", "compute_rand_bootstrap_confidence_interval_studentized", "compute_rand_bootstrap_confidence_interval_symmetric-percentile-t", "compute_rand_bootstrap_confidence_interval_smoothed"), notes="AsympLikStdModCache; rand CI avg 116.6s / p80 46.1s at n=22 slow; selected BRT CI paths avg >30s slow"),
  list(name="InferencePropZeroOneInflatedBetaRegr", section="Proportion", resp="prop",  kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp", pboot=TRUE,  notes="AsympLikStdModCache; pboot=TRUE (separate from bootstrap)"),
  list(name="InferencePropGCompMeanDiff",           section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=TRUE,  rand_resp="csp", skip_rpv=TRUE,  skip_rci=TRUE,  rci_resp="cp", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferencePropFractionalLogit",         section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp", pboot=NA,    notes="AsympLikStdModCacheNoParamBootstrap; supports_likelihood_tests=FALSE"),

  # ── COUNT ────────────────────────────────────────────────────────────────────
  list(name="InferenceCountPoissonKKGEE",           section="Count",      resp="count", kk=TRUE,  types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA,    notes="InferenceAsymp; rand/rci N/A for count"),
  list(name="InferenceCountKKGLMM",                 section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE, skip_jack_slow=TRUE, skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; supports_lik_ratio_param_bootstrap=TRUE when use_rcpp; jackknife hard-excluded in comprehensive_tests.R (supports_jackknife name exclusion), displayed as SLOW"),
  list(name="InferenceCountKKHurdlePoissonOneLik",  section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; simulate_under_lik_null: ZTP+GLMM matched + Poisson reservoir; in skip_ci_rand; rand N/A count"),
  list(name="InferenceCountKKCondPoissonOneLik",    section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; simulate_under_lik_null: cond-Poisson pairs Binomial + reservoir Poisson; in skip_ci_rand; rand N/A count"),
  list(name="InferenceCountPoisson",                section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="CountLikelihood → AsympLikStdModCache"),
  list(name="InferenceCountRobustPoisson",          section="Count",      resp="count", kk=FALSE, types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA, notes="CountCompositeLikelihood; explicit FALSE"),
  list(name="InferenceCountQuasiPoisson",           section="Count",      resp="count", kk=FALSE, types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA, notes="CountCompositeLikelihood; explicit FALSE"),
  list(name="InferenceCountNegBin",                 section="Count",      resp="count", kk=FALSE, types="full",skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE, notes="InferenceCountLikelihood → InferenceParamBootstrap; simulate_under_lik_null: rnbinom draw"),
  list(name="InferenceCountZeroInflatedPoisson",    section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; use_rcpp; 'Zero-Inflated Poisson' in pboot list; in skip_ci_rand; raw-LR bootstrap + Bartlett approx carve-out removed after 900-replicate calibration stress test found no miscalibration"),
  list(name="InferenceCountZeroInflatedNegBin",     section="Count",      resp="count", kk=FALSE, types="full",  skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; use_rcpp; 'Zero-Inflated Negative Binomial' in pboot list"),
  list(name="InferenceCountHurdlePoisson",          section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; 'Hurdle Poisson' in pboot list; in skip_ci_rand; raw-LR bootstrap + Bartlett approx carve-out removed after 900-replicate calibration stress test found no miscalibration"),
  list(name="InferenceCountHurdleNegBin",           section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE, notes="InferenceCountLikelihood → InferenceParamBootstrap; simulate_under_lik_null implemented; in skip_ci_rand"),

  # ── SURVIVAL ─────────────────────────────────────────────────────────────────
  list(name="InferenceSurvivalGehanWilcox",            section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    skip_brt="ci", slow_methods=c("compute_bootstrap_confidence_interval_studentized"), notes="InferenceAsymp; all BRT CI types avg 33-39s slow; bootstrap studentized CI max 10443.4s slow"),
  list(name="InferenceSurvivalLogRank",                section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalRestrictedMeanDiff",     section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalKMDiff",                 section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalWeibullRegr",            section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  skip_brt="ci", slow_methods=c("compute_rand_confidence_interval"), notes="AsympLikStdModCache → ParamBootstrap; all BRT CI types avg 42-518s slow; rand CI avg >30s slow"),
  list(name="InferenceSurvivalDepCensTransformRegr",   section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  slow_methods=c("compute_rand_bootstrap_two_sided_pval_smoothed"), notes="AsympLikStdModCache; simulate_under_lik_null added; BRT smoothed pval avg >30s slow"),
  list(name="InferenceSurvivalCoxPHRegr",              section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AsympLikStdModCache; pboot=use_rcpp (default)"),
  list(name="InferenceSurvivalStratCoxPHRegr",         section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="InferenceParamBootstrap directly; pboot=use_rcpp"),
  list(name="InferenceSurvivalKKClaytonCopulaOneLik",  section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE, skip_jack_slow=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=TRUE,  skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  skip_pboot_ci=TRUE, slow_methods=c("compute_rand_confidence_interval", "compute_rand_confidence_interval(custom)", "compute_lik_ratio_confidence_interval"), notes="InferenceParamBootstrap; pboot=TRUE even though skip_boot=TRUE (separate test family); jackknife structurally supported but skip_jack_slow=TRUE; rand CI avg 507.7s / max 2094.4s at n=6; custom rand CI avg 41.9s / max 1993.3s at n=53; lik-ratio CI avg 39.1s / max 234.2s at n=6 slow"),
  list(name="InferenceSurvivalKKLWACoxPHOneLik",       section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AbstractKKLWACoxOneLik → ParamBootstrap; explicit TRUE"),
  list(name="InferenceSurvivalKKStratCoxPHOneLik",     section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="InferenceParamBootstrap; explicit TRUE"),
  list(name="InferenceSurvivalKKWeibullFrailtyOneLik", section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  slow_methods=c("compute_rand_confidence_interval", "compute_score_confidence_interval"), notes="AbstractKKWeibullFrailtyOneLik → ParamBootstrap; pboot=use_rcpp; skip_boot=TRUE; rand CI avg 45.7s / max 171.5s at n=6; score CI avg 50.1s / max 324.8s at n=13 slow"),
  list(name="InferenceSurvivalKKWeibullMarginal",      section="Survival",  resp="surv", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),

  # ── ORDINAL ───────────────────────────────────────────────────────────────────
  list(name="InferenceOrdinalKKGEE",                    section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp; rand/rci N/A for ordinal"),
  list(name="InferenceOrdinalKKGLMM",                   section="Ordinal", resp="ord", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; simulate_under_lik_null added; supports=isTRUE(use_rcpp); BRT pval smoothed avg 317s slow"),
  list(name="InferenceOrdinalKKCLMM",                   section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="AbstractKKOrdinalCLMM → AsympLik; supports_likelihood_tests=FALSE"),
  list(name="InferenceOrdinalKKCLMMProbit",             section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="AbstractKKOrdinalCLMM → AsympLik"),
  list(name="InferenceOrdinalKKCLMMCauchit",            section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="AbstractKKOrdinalCLMM → AsympLik"),
  list(name="InferenceOrdinalKKCLMMCloglog",            section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="AbstractKKOrdinalCLMM → AsympLik"),
  list(name="InferenceOrdinalKKCondAdjCatLogitRegr",    section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsympLik; supports_likelihood_tests=FALSE"),
  list(name="InferenceOrdinalPairedSignTest",           section="Ordinal", resp="ord", kk=TRUE,  types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsympLik; supports_likelihood_tests=FALSE"),
  list(name="InferenceOrdinalAdjCatLogitRegr",          section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  slow_methods=c("compute_rand_bootstrap_two_sided_pval_smoothed"), notes="AsympLikStdModCache; simulate_under_lik_null added; BRT smoothed pval avg 499s / max 2694s slow"),
  list(name="InferenceOrdinalStereotypeLogitRegr",      section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  slow_methods=c("compute_rand_bootstrap_two_sided_pval_smoothed"), notes="AsympLikStdModCache; simulate_under_lik_null added; BRT smoothed pval avg >30s slow"),
  list(name="InferenceOrdinalContRatioRegr",            section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache; simulate_under_lik_null added; BRT pval smoothed avg 66s slow"),
  list(name="InferenceOrdinalCloglogRegr",              section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE even with skip_boot"),
  list(name="InferenceOrdinalCauchitRegr",              section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE"),
  list(name="InferenceOrdinalOrderedProbitRegr",        section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE"),
  list(name="InferenceOrdinalPropOddsRegr",             section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceOrdinalPartialProportionalOddsRegr",section="Ordinal",resp="ord",kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceOrdinalGCompMeanDiff",            section="Ordinal", resp="ord", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceOrdinalJonckheereTerpstraTest",   section="Ordinal", resp="ord", kk=FALSE, types="none", skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=TRUE,  skip_bbt=TRUE,  jack=FALSE, skip_rand=TRUE,  rand_resp="",        skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, notes="InferenceAsymp; SPECIAL: only exact pval + estimate called"),
  list(name="InferenceOrdinalRidit",                   section="Ordinal", resp="ord", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp")
)

# ── HTML renderer ─────────────────────────────────────────────────────────────
load_nonestimability_stats = function(result_dir = "package_tests") {
  defaults_file = file.path(result_dir, "path_audits_nonestimability_defaults.csv")
  stats_env_from_rates = function(rates) {
    stats_env = new.env(parent = emptyenv())
    required = c("base_class", "function_run", "n", "nonestimable", "rate")
    if (is.null(rates) || !all(required %in% names(rates))) return(stats_env)
    for (i in seq_len(nrow(rates))) {
      assign(
        paste(rates$base_class[i], rates$function_run[i], sep = "||"),
        list(
          n = as.integer(rates$n[i]),
          nonestimable = as.integer(rates$nonestimable[i]),
          rate = as.numeric(rates$rate[i])
        ),
        envir = stats_env
      )
    }
    stats_env
  }
  read_default_rates = function() {
    if (!file.exists(defaults_file)) return(NULL)
    tryCatch(utils::read.csv(defaults_file, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  load_default_stats = function() {
    stats_env_from_rates(read_default_rates())
  }
  maybe_update_default_stats = function(rates) {
    if (is.null(rates) || !all(c("n", "base_class", "function_run", "nonestimable", "rate") %in% names(rates))) {
      return(invisible(FALSE))
    }
    default_rates = read_default_rates()
    current_total = sum(as.numeric(rates$n), na.rm = TRUE)
    default_total = if (is.null(default_rates) || !"n" %in% names(default_rates)) {
      -Inf
    } else {
      sum(as.numeric(default_rates$n), na.rm = TRUE)
    }
    if (is.finite(current_total) && current_total > default_total) {
      utils::write.csv(rates, defaults_file, row.names = FALSE)
      return(invisible(TRUE))
    }
    invisible(FALSE)
  }
  files = list.files(result_dir, pattern = "^comprehensive_tests_results.*\\.csv$", full.names = TRUE)
  files = files[!grepl("_filtered_", basename(files), fixed = TRUE)]
  if (!length(files) || !requireNamespace("data.table", quietly = TRUE)) {
    return(load_default_stats())
  }
  pieces = lapply(files, function(f) {
    tryCatch(data.table::fread(f, showProgress = FALSE, fill = TRUE), error = function(e) NULL)
  })
  pieces = Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(load_default_stats())
  dt = data.table::rbindlist(pieces, fill = TRUE)
  required = c("inference_class", "function_run", "error_message")
  if (!all(required %in% names(dt))) return(load_default_stats())
  dt[, base_class := sub(" [\\(\\[].*$", "", inference_class)]
  dt[, explicit_nonestimable := grepl("^Explicitly non-estimable", error_message)]
  rates = dt[, .(
    n = .N,
    nonestimable = sum(explicit_nonestimable, na.rm = TRUE)
  ), by = .(base_class, function_run)]
  rates[, rate := data.table::fifelse(n > 0L, nonestimable / n, NA_real_)]
  data.table::setorder(rates, base_class, function_run)
  maybe_update_default_stats(rates)
  stats_env_from_rates(rates)
}

html_from_audit = function(classes, outfile = "path_audits.html") {
  GREEN = "#2e7d32"; LIGHT_GREEN = "#a5d6a7"; PALE_GREEN = "#dff3df"; LIME = "#dce775"; YELLOW = "#fff176"; ORANGE_LIGHT = "#ffcc80"; ORANGE = "#ff9800"; ORANGE_DARK = "#ef6c00"; LIGHT_RED = "#ffcdd2"; DARK_GREY = "#333333"; GREY = "#e0e0e0"
  `%||%` = function(x, y) if (is.null(x)) y else x
  nonestimability_stats = load_nonestimability_stats(dirname(outfile))
  html_escape = function(x) {
    x = gsub("&", "&amp;", x, fixed = TRUE)
    x = gsub("<", "&lt;", x, fixed = TRUE)
    x = gsub(">", "&gt;", x, fixed = TRUE)
    x = gsub('"', "&quot;", x, fixed = TRUE)
    x
  }
  low_estimability_rows = function() {
    audited_names = vapply(classes, function(r) r$name, character(1))
    keys = ls(nonestimability_stats, all.names = TRUE)
    keys = keys[grepl("\\|\\|", keys)]
    rows = lapply(keys, function(key) {
      info = get(key, envir = nonestimability_stats, inherits = FALSE)
      parts = strsplit(key, "\\|\\|")[[1]]
      if (length(parts) != 2L || !parts[1] %in% audited_names) return(NULL)
      estimable_rate = 1 - info$rate
      if (!is.finite(estimable_rate) || estimable_rate >= 0.01) return(NULL)
      data.frame(
        base_class = parts[1],
        function_run = parts[2],
        n = as.integer(info$n),
        nonestimable = as.integer(info$nonestimable),
        estimable_rate = as.numeric(estimable_rate),
        stringsAsFactors = FALSE
      )
    })
    rows = Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(base_class = character(), function_run = character(), n = integer(), nonestimable = integer(), estimable_rate = numeric(), stringsAsFactors = FALSE))
    }
    rows = do.call(rbind, rows)
    rows[order(rows$base_class, rows$function_run), , drop = FALSE]
  }
  low_estimability_html = function() {
    rows = low_estimability_rows()
    if (!nrow(rows)) {
      return('<h2>&lt;1% Estimable Paths</h2><p>No audited paths are currently below 1% observed estimability.</p>')
    }
    row_html = vapply(seq_len(nrow(rows)), function(i) {
      sprintf(
        '<tr><td style="font-family:monospace;white-space:nowrap">%s</td><td style="font-family:monospace;white-space:nowrap">%s</td><td style="text-align:right">%d</td><td style="text-align:right">%d</td><td style="text-align:right">%.2f%%</td></tr>',
        html_escape(rows$base_class[i]), html_escape(rows$function_run[i]), rows$n[i], rows$nonestimable[i], 100 * rows$estimable_rate[i]
      )
    }, character(1))
    paste0(
      '<h2>&lt;1% Estimable Paths</h2>',
      '<p style="color:#555">Audited class/function paths with observed estimability below 1%, using the same current/default comprehensive results snapshot as the color scale.</p>',
      '<table class="low-estimability"><thead><tr><th>Inference Type</th><th>Function</th><th>n</th><th>Non-Estimable</th><th>Estimable</th></tr></thead><tbody>',
      paste(row_html, collapse = "\n"),
      '</tbody></table>'
    )
  }
  cell = function(color, s, title = NULL, text_color = NULL) {
    title_attr = if (is.null(title)) "" else sprintf(' title="%s"', title)
    style = sprintf("background:%s;text-align:center;white-space:nowrap", color)
    if (!is.null(text_color)) style = paste0(style, ";color:", text_color)
    sprintf('<td style="%s"%s>%s</td>', style, title_attr, s)
  }
  status_cell = function(status, info = NULL) {
    rate_title = function(default_title) {
      if (is.null(info)) return(default_title)
      estimable_rate = 1 - info$rate
      sprintf(
        "%s; observed approx estimable rate %.1f%%; explicit non-estimable rate %.1f%% (%d/%d rows)",
        default_title, 100 * estimable_rate, 100 * info$rate, as.integer(info$nonestimable), as.integer(info$n)
      )
    }
    switch(status,
      always_numeric = cell(GREEN,       "✓",        "Theoretically guaranteed by the comprehensive-test contract", text_color = "#fff"),
      maybe_unknown  = cell(PALE_GREEN,  "unknown",  "Attempted, empirical estimability unknown; no result rows found"),
      maybe_0        = cell(LIGHT_GREEN, "100%",     rate_title("Attempted, 100% estimable observed")),
      maybe_low      = cell(PALE_GREEN,  "[95-100)%",rate_title("Attempted, [95-100)% estimable observed")),
      maybe_mid      = cell(LIME,        "[75-95)%", rate_title("Attempted, [75-95)% estimable observed")),
      maybe_high     = cell(YELLOW,      "[25-75)%", rate_title("Attempted, [25-75)% estimable observed")),
      maybe_5_25     = cell(ORANGE_LIGHT,"[5-25)%",  rate_title("Attempted, [5-25)% estimable observed")),
      maybe_1_5      = cell(ORANGE,      "[1-5)%",   rate_title("Attempted, [1-5)% estimable observed")),
      maybe_lt_1     = cell(ORANGE_DARK, "(0-1)%",   rate_title("Attempted, (0-1)% estimable observed"), text_color = "#fff"),
      maybe_zero     = cell(ORANGE_DARK, "0%",       rate_title("Attempted, 0% estimable observed"), text_color = "#fff"),
      slow           = cell(LIGHT_RED,   "SLOW", "Skipped or lightly tested in comprehensive tests because this path is too slow"),
      unsupported    = cell(DARK_GREY,   "NTS",  "Not theoretically supported by this model", text_color = "#fff"),
      not_implemented= cell(GREY,        "NI",   "Not implemented yet"),
      stop("Unknown audit cell status: ", status)
    )
  }
  ok    = function() status_cell("always_numeric")
  maybe = function() status_cell("maybe_unknown")
  bad   = function() status_cell("slow")
  na    = function() status_cell("unsupported")
  ni    = function() status_cell("not_implemented")

  row_methods = function(r, field) {
    as.character(r[[field]] %||% character())
  }
  method_status = function(r, method_id, default_status) {
    if (method_id %in% row_methods(r, "unsupported_methods")) return("unsupported")
    if (method_id %in% row_methods(r, "slow_methods")) return("slow")
    if (method_id %in% row_methods(r, "always_numeric_methods")) return("always_numeric")
    if (method_id %in% row_methods(r, "maybe_nonestimable_methods")) return("maybe")
    default_status
  }
  method_info = function(r, method_id) {
    key = paste(r$name, method_id, sep = "||")
    if (exists(key, envir = nonestimability_stats, inherits = FALSE)) {
      return(get(key, envir = nonestimability_stats, inherits = FALSE))
    }
    NULL
  }
  maybe_status = function(info) {
    if (is.null(info) || !is.finite(info$rate)) return("maybe_unknown")
    if (info$rate <= 0) return("maybe_0")
    if (info$rate <= 0.05) return("maybe_low")
    if (info$rate <= 0.25) return("maybe_mid")
    if (info$rate <= 0.75) return("maybe_high")
    if (info$rate < 0.95) return("maybe_5_25")
    if (info$rate < 0.99) return("maybe_1_5")
    if (info$rate < 1) return("maybe_lt_1")
    "maybe_zero"
  }
  method_cell = function(r, method_id, default_status) {
    status = method_status(r, method_id, default_status)
    info = if (identical(status, "maybe")) method_info(r, method_id) else NULL
    if (identical(status, "maybe")) status = maybe_status(info)
    status_cell(status, info)
  }

  stable_model_based_numeric = function(r, ttype) {
    r$name %in% c(
      "InferenceAllSimpleMeanDiff",
      "InferenceAllSimpleMeanDiffPooledVar",
      "InferenceAllSimpleWilcox",
      "InferenceIncidWald",
      "InferenceIncidMiettinenNurminenRiskDiff",
      "InferenceIncidNewcombeRiskDiff",
      "InferenceOrdinalJonckheereTerpstraTest",
      "InferenceOrdinalRidit"
    ) && ttype == "wald"
  }

  # ── Model estimate ("est", nested under Asymptotic > Model-Based) ───────────
  # compute_estimate() is defined on every audited class, including the exact-only
  # classes (InferenceIncidExactFisher/Binomial, InferenceIncidenceExactZhang) and
  # InferenceOrdinalJonckheereTerpstraTest -- verified against EDI/R source. Never
  # unsupported/slow; only always_numeric (stable closed-form classes) vs maybe --
  # same color logic as the model-based wald/score/lr/grad pval/ci cells below.
  cell_estimate = function(r) {
    default = if (stable_model_based_numeric(r, "wald")) "always_numeric" else "maybe"
    method_cell(r, "compute_estimate", default)
  }

  # ── Asymptotic ──────────────────────────────────────────────────────────────
  type_ok = function(r, ttype) {
    if (isTRUE(r$skip_asymp)) return(FALSE)
    switch(r$types,
      full            = TRUE,
      wald            = ttype == "wald",
      wald_lr         = ttype %in% c("wald","lr"),
      wald_score_grad = ttype %in% c("wald","score","grad"),
      wald_score_lr   = ttype %in% c("wald","score","lr"),
      none            = FALSE,
      FALSE)
  }
  method_id_for_model_based = function(ttype, kind) {
    prefix = switch(ttype,
      wald = "compute_wald",
      score = "compute_score",
      lr = "compute_lik_ratio",
      grad = "compute_gradient")
    suffix = if (kind == "pval") "two_sided_pval" else "confidence_interval"
    paste(prefix, suffix, sep = "_")
  }
  cell_ap = function(r, ttype) {
    method_id = method_id_for_model_based(ttype, "pval")
    if (!type_ok(r, ttype)) return(method_cell(r, method_id, "unsupported"))
    default = if (stable_model_based_numeric(r, ttype)) "always_numeric" else "maybe"
    method_cell(r, method_id, default)
  }
  cell_ac = function(r, ttype) {
    method_id = method_id_for_model_based(ttype, "ci")
    if (!type_ok(r, ttype)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    default = if (stable_model_based_numeric(r, ttype)) "always_numeric" else "maybe"
    method_cell(r, method_id, default)
  }
  # Approximate Bartlett-corrected likrat ("LR-Bart-app", nested next to "LR" under
  # Model-Based): base-class plumbing lives in InferenceAsympLik
  # (compute_lik_ratio_bartlett_approx_two_sided_pval / _confidence_interval); InferenceParamBootstrap
  # overrides supports_bartlett_likelihood_ratio_approx() to isTRUE(supports_lik_ratio_param_bootstrap())
  # -- i.e. Monte-Carlo Bartlett support automatically follows parametric-bootstrap LR support
  # (the same "pboot" field already tracked per row), reusing the simulate_under_lik_null()
  # machinery for free. Classes without any likelihood/partial-likelihood LR test at all (same
  # type_ok("lr") gate used by the "LR" column itself) are structurally NTS here, same as "LR".
  # Among the remainder, rows follow pboot exactly: pboot=TRUE -> maybe, pboot=FALSE/NA -> NI.
  cell_likrat_bart_p = function(r) {
    method_id = "compute_lik_ratio_bartlett_two_sided_pval"
    if (!type_ok(r, "lr")) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$bartlett_exact)) return(method_cell(r, method_id, "not_implemented"))
    if (!isTRUE(r$pboot)) return(method_cell(r, method_id, "not_implemented"))
    method_cell(r, method_id, "maybe")
  }
  cell_likrat_bart_c = function(r) {
    method_id = "compute_lik_ratio_bartlett_confidence_interval"
    if (!type_ok(r, "lr")) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$bartlett_exact)) return(method_cell(r, method_id, "not_implemented"))
    if (!isTRUE(r$pboot)) return(method_cell(r, method_id, "not_implemented"))
    if (!isTRUE(r$run_pboot_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  # "LR-Bart-ex": exact-factor (closed-form analytic) Bartlett correction, nested
  # next to "LR-Bart-app". Implemented so far only for InferenceContinKKOLSOneLik
  # (bartlett_exact=TRUE): holding sigma2 fixed makes the package's own LR statistic
  # algebraically identical to the classical partial F(1,n-p) statistic, an exact
  # finite-sample pivot (verified against base R's lm()), not a Cordeiro-style
  # tensor derivation. Every other family remains NI -- see
  # package_metadata/likrat_correction_bartlett.md's practical-derivation-risk table.
  cell_likrat_bart_ex_p = function(r) {
    method_id = "compute_lik_ratio_bartlett_two_sided_pval"
    if (!type_ok(r, "lr")) return(method_cell(r, method_id, "unsupported"))
    if (!isTRUE(r$bartlett_exact)) return(method_cell(r, method_id, "not_implemented"))
    method_cell(r, method_id, "maybe")
  }
  cell_likrat_bart_ex_c = function(r) {
    method_id = "compute_lik_ratio_bartlett_confidence_interval"
    if (!type_ok(r, "lr")) return(method_cell(r, method_id, "unsupported"))
    if (!isTRUE(r$bartlett_exact)) return(method_cell(r, method_id, "not_implemented"))
    if (!isTRUE(r$run_pboot_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # "other": class-specific asymptotic pval methods with non-generic names, not
  # reachable through the wald/score/lr/grad dispatch above. No CI counterpart exists.
  is_log_rank_other = function(r) {
    r$name %in% c("InferenceSurvivalLogRank", "InferenceSurvivalKMDiff")
  }
  cell_am_other_p = function(r) {
    method_id = "compute_asymp_log_rank_two_sided_pval_for_treatment_effect"
    if (!is_log_rank_other(r) || isTRUE(r$skip_asymp)) return(method_cell(r, method_id, "unsupported"))
    method_cell(r, method_id, "maybe")
  }

  # ── Nonparam bootstrap (classical) ─────────────────────────────────────────
  cell_bp = function(r, method_id) {
    if (isTRUE(r$skip_boot)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_bc = function(r, method_id) {
    if (isTRUE(r$skip_boot)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_boot_ci)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Bayesian bootstrap ──────────────────────────────────────────────────────
  cell_bbt_p = function(r, method_id) {
    if (is.na(r$skip_bbt)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_bbt)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_bbt_c = function(r, method_id) {
    if (is.na(r$skip_bbt)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_bbt)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_bbt_ci)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Studentized bootstrap (separate skip from general boot) ──────────────────
  cell_bst_p = function(r, method_id) {
    if (isTRUE(r$skip_boot)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_stud)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_bst_c = function(r, method_id) {
    if (isTRUE(r$skip_boot)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_stud)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Parametric bootstrap ────────────────────────────────────────────────────
  # pboot=FALSE means the class inherits InferenceParamBootstrap but
  # supports_lik_ratio_param_bootstrap()=FALSE because simulate_under_lik_null()
  # was never implemented for it -- structurally not implemented (NI), not slow.
  cell_pb_p = function(r, method_id) {
    if (is.na(r$pboot)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$pboot)) method_cell(r, method_id, "maybe") else method_cell(r, method_id, "not_implemented")
  }
  cell_pb_c = function(r, method_id) {
    if (is.na(r$pboot)) return(method_cell(r, method_id, "unsupported"))
    if (!isTRUE(r$pboot)) return(method_cell(r, method_id, "not_implemented"))
    if (!isTRUE(r$run_pboot_ci)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_pboot_ci)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Parametric bootstrap estimate (blank subheader, nested left of "LR") ────
  # compute_param_bootstrap_estimate()/compute_param_bootstrap_pval()/
  # compute_param_bootstrap_confidence_interval() (InferenceParamBootstrap):
  # single-level parametric-bootstrap bias correction of the point estimate
  # (2*theta_hat - mean(theta_star)) plus its "basic"/reflected-quantile pval
  # and CI companions -- all three share the same simulate_under_lik_null()
  # contract as the "LR" pair and are gated by the same pboot field (they
  # default to isTRUE(supports_lik_ratio_param_bootstrap())). Per explicit
  # instruction: is.na(pboot) (not InferenceParamBootstrap at all) -> NTS;
  # everything else (pboot TRUE or FALSE) is left blank pending audit.
  cell_pbest_blank = function() '<td></td>'
  cell_pbest_estimate = function(r) {
    method_id = "compute_param_bootstrap_estimate"
    if (is.na(r$pboot)) return(method_cell(r, method_id, "unsupported"))
    cell_pbest_blank()
  }
  cell_pbest_p = function(r) {
    method_id = "compute_param_bootstrap_pval"
    if (is.na(r$pboot)) return(method_cell(r, method_id, "unsupported"))
    cell_pbest_blank()
  }
  cell_pbest_c = function(r) {
    method_id = "compute_param_bootstrap_confidence_interval"
    if (is.na(r$pboot)) return(method_cell(r, method_id, "unsupported"))
    cell_pbest_blank()
  }

  # ── Jackknife ───────────────────────────────────────────────────────────────
  # jack=FALSE means truly unsupported (no such row currently). jack=TRUE + skip_jack_slow=TRUE
  # means skipped in comprehensive tests -- either for runtime (InferenceSurvivalKKClaytonCopulaOneLik,
  # InferenceContinKKGLMM) or via a hard-coded name exclusion (InferenceCountKKGLMM).
  cell_jack_estimate = function(r) {
    method_id = "compute_jackknife_estimate"
    if (!isTRUE(r$jack)) return(method_cell(r, method_id, "unsupported"))
    if (identical(r$resp, "count")) return(method_cell(r, method_id, "not_implemented"))
    if (isTRUE(r$skip_jack_slow)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_jack_p = function(r, method_id) {
    if (!isTRUE(r$jack)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_jack_slow)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_jack_c = function(r, method_id) {
    if (!isTRUE(r$jack)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_jack_slow)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Randomization ───────────────────────────────────────────────────────────
  cell_rand_p = function(r, method_id) {
    if (!nchar(r$rand_resp)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_rand_c = function(r, method_id) {
    if (!nchar(r$rci_resp)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_rci) || isTRUE(r$skip_rand_ci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_rand_custom_p = function(r, method_id) {
    if (!nchar(r$rand_resp)) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_rand_custom_c = function(r, method_id) {
    if (!r$resp %in% c("cont", "prop", "surv")) return(method_cell(r, method_id, "unsupported"))
    if (isTRUE(r$skip_rci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Bootstrap-randomization (BRT) — all types share same per-class logic ───
  # skip_brt="ci" means pval runs but CI is too slow (e.g. ContinRobustRegr)
  cell_brt_p = function(r, method_id) {
    if (!nchar(r$rand_resp)) return(method_cell(r, method_id, "unsupported"))
    if (!isTRUE(r$run_brt)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(method_cell(r, method_id, "slow"))
    if (!identical(r$skip_brt, "ci") && isTRUE(r$skip_boot)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }
  cell_brt_c = function(r, method_id) {
    if (!nchar(r$rci_resp)) return(method_cell(r, method_id, "unsupported"))
    if (!isTRUE(r$run_brt)) return(method_cell(r, method_id, "slow"))
    if (identical(r$skip_brt, "ci")) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_boot) || isTRUE(r$skip_rand)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_rand_ci)) return(method_cell(r, method_id, "slow"))
    if (isTRUE(r$skip_ci) || isTRUE(r$skip_rci)) return(method_cell(r, method_id, "slow"))
    method_cell(r, method_id, "maybe")
  }

  # ── Exact inference (Fisher/Binomial/Zhang: pval+CI; Jonckheer: pval only) ───
  cell_exact_p = function(r) {
    method_id = "compute_exact_two_sided_pval_for_treatment_effect"
    if (isTRUE(r$exact_p)) method_cell(r, method_id, "always_numeric") else method_cell(r, method_id, "unsupported")
  }
  cell_exact_c = function(r) {
    method_id = "compute_exact_confidence_interval"
    if (isTRUE(r$exact_c)) method_cell(r, method_id, "always_numeric") else method_cell(r, method_id, "unsupported")
  }

  # ── Four frozen header rows ──────────────────────────────────────────────────
  # Data cols: 14 (model-based, incl. 1 model "est", 1 log-rank "other", and the
  #            2-col "LR-Bart-app" + 2-col "LR-Bart-ex" plumbing-only pairs next to "LR")
  #            + 2 (exact-other) + 5 (exact-rand) + 18 (npboot, incl. 10 for the 6 bayes flavors)
  #            + 5 (pboot: 3-col parametric-bootstrap-estimate [blank subheader, est/pval/ci]
  #                + 2-col "LR" pval/ci) + 3 (jack, incl. 1 jackknife "est") + 9 (brt) = 56
  NCOL = 57  # 1 (class) + 56 (data)
  hdr = paste0(
    # Row 0: meta-categories — Asymptotic (49) | Exact (7)
    '<tr class="hdr0">',
      '<th rowspan="4" style="text-align:left">Inference Type</th>',
      '<th colspan="49">Asymptotic</th>',
      '<th colspan="7">Exact</th>',
    '</tr>',
    # Row 1: categories
    '<tr class="hdr1">',
      # Under Asymptotic (44)
      '<th colspan="14">Model-Based</th>',
      '<th colspan="18">Nonparam Boot</th>',
      '<th colspan="5">Param Boot</th>',
      '<th colspan="3">Jackknife</th>',
      '<th colspan="9">Boot-Rand</th>',
      # Under Exact (7)
      '<th colspan="2">other</th>',
      '<th colspan="5">Randomization</th>',
    '</tr>',
    # Row 2: test types
    '<tr class="hdr2">',
      # Model-Based: model "est" (colspan=1), then 4 types × colspan=2, plus "LR-Bart-app"
      # and "LR-Bart-ex" (colspan=2 each, plumbing only, always NI) and "other" (pval-only, colspan=1)
      '<th colspan="1">est</th>',
      '<th colspan="2">wald</th>',
      '<th colspan="2">score</th>',
      '<th colspan="2">LR</th>',
      '<th colspan="2">LR-Bart-app</th>',
      '<th colspan="2">LR-Bart-ex</th>',
      '<th colspan="2">grad</th>',
      '<th colspan="1">other</th>',
      # Nonparam Boot: pctile(2) symm(1) basic(1) bca(2) stud(2) [classical, 8]
      #                bayes-pctile(2) bayes-symm(1) bayes-basic(1) bayes-wald(2) bayes-bca(2) bayes-stud(2) [bayes, 10]
      '<th colspan="2">pctile</th>',
      '<th colspan="1">symm</th>',
      '<th colspan="1">basic</th>',
      '<th colspan="2">bca</th>',
      '<th colspan="2">stud</th>',
      '<th colspan="2">bayes-pctile</th>',
      '<th colspan="1">bayes-symm</th>',
      '<th colspan="1">bayes-basic</th>',
      '<th colspan="2">bayes-wald</th>',
      '<th colspan="2">bayes-bca</th>',
      '<th colspan="2">bayes-stud</th>',
      # Param Boot: blank subheader (parametric-bootstrap-estimate: est/pval/ci) then "LR"
      '<th colspan="3"></th>',
      '<th colspan="2">LR</th>',
      # Jackknife (no named sub-type)
      '<th colspan="3"></th>',
      # Boot-Rand: 4 types
      '<th colspan="3">pctile</th>',
      '<th colspan="2">stud</th>',
      '<th colspan="2">symm-t</th>',
      '<th colspan="2">smooth</th>',
      # Exact other (no sub-type)
      '<th colspan="2"></th>',
      # Randomization: vanilla(3) + custom(2)
      '<th colspan="3">vanilla</th>',
      '<th colspan="2">custom</th>',
    '</tr>',
    # Row 3: pval / ci leaves
    '<tr class="hdr3">',
      # Model-Based (14)
      '<th>est</th>',
      '<th>pval</th><th>ci</th>',   # wald
      '<th>pval</th><th>ci</th>',   # score
      '<th>pval</th><th>ci</th>',   # likrat
      '<th>pval</th><th>ci</th>',   # LR-Bart-app (plumbing only, always NI)
      '<th>pval</th><th>ci</th>',   # LR-Bart-ex (plumbing only, always NI)
      '<th>pval</th><th>ci</th>',   # grad
      '<th>pval</th>',              # other (log-rank, pval only)
      # Nonparam Boot (18)
      '<th>pval</th><th>ci</th>',   # pctile
      '<th>pval</th>',               # symm (pval only)
      '<th>ci</th>',                 # basic (ci only)
      '<th>pval</th><th>ci</th>',   # bca
      '<th>pval</th><th>ci</th>',   # stud
      '<th>pval</th><th>ci</th>',   # bayes-pctile
      '<th>pval</th>',               # bayes-symm (pval only)
      '<th>ci</th>',                 # bayes-basic (ci only)
      '<th>pval</th><th>ci</th>',   # bayes-wald
      '<th>pval</th><th>ci</th>',   # bayes-bca
      '<th>pval</th><th>ci</th>',   # bayes-stud
      # Param Boot (5): parametric-bootstrap-estimate (est/pval/ci) then LR (pval/ci)
      '<th>est</th><th>pval</th><th>ci</th>',
      '<th>pval</th><th>ci</th>',
      # Jackknife (3)
      '<th>est</th><th>pval</th><th>ci</th>',
      # Boot-Rand (9)
      '<th>pval</th><th>pval(&delta;)</th><th>ci</th>',   # pctile
      '<th>pval</th><th>ci</th>',                          # stud
      '<th>pval</th><th>ci</th>',                          # symm-t
      '<th>pval</th><th>ci</th>',                          # smooth
      # Exact other (2)
      '<th>pval</th><th>ci</th>',
      # Randomization (5)
      '<th>pval</th><th>pval(&delta;)</th><th>ci</th>',   # vanilla
      '<th>pval</th><th>ci</th>',                          # custom
    '</tr>'
  )

  sections = unique(vapply(classes, `[[`, "", "section"))
  body = ""
  for (sec in sections) {
    body = paste0(body, sprintf(
      '<tr><td colspan="%d" style="background:#263238;color:white;padding:4px 8px;font-weight:bold">%s</td></tr>\n',
      NCOL, sec))
    for (r in classes[vapply(classes, function(x) x$section == sec, logical(1))]) {
      nm = sub("^Inference", "", r$name)
      body = paste0(body, "<tr>",
        sprintf('<td style="font-family:monospace;padding:2px 8px;white-space:nowrap">%s</td>', nm),
        # Model-Based (14): model "est" first, then wald/score/lr/LR-Bart-app/LR-Bart-ex/grad pval+ci, then "other"
        cell_estimate(r),
        cell_ap(r,"wald"),  cell_ac(r,"wald"),
        cell_ap(r,"score"), cell_ac(r,"score"),
        cell_ap(r,"lr"),    cell_ac(r,"lr"),
        cell_likrat_bart_p(r), cell_likrat_bart_c(r),
        cell_likrat_bart_ex_p(r), cell_likrat_bart_ex_c(r),
        cell_ap(r,"grad"),  cell_ac(r,"grad"),
        cell_am_other_p(r),
        # Nonparam Boot (18)
        cell_bp(r, "compute_bootstrap_two_sided_pval"), cell_bc(r, "compute_bootstrap_confidence_interval"),        # pctile p + ci
        cell_bp(r, "compute_bootstrap_two_sided_pval_symmetric"),                                                   # symm p (no ci type)
        cell_bc(r, "compute_bootstrap_confidence_interval_basic"),                                                   # basic ci (no p type)
        cell_bp(r, "compute_bootstrap_two_sided_pval_bca"), cell_bc(r, "compute_bootstrap_confidence_interval_bca"), # bca p + ci
        cell_bst_p(r, "compute_bootstrap_two_sided_pval_studentized"), cell_bst_c(r, "compute_bootstrap_confidence_interval_studentized"), # stud p + ci
        cell_bbt_p(r, "compute_bayesian_bootstrap_two_sided_pval"), cell_bbt_c(r, "compute_bayesian_bootstrap_confidence_interval"),        # bayes-pctile p + ci
        cell_bbt_p(r, "compute_bayesian_bootstrap_two_sided_pval_symmetric"),                                                                  # bayes-symm p (no ci type)
        cell_bbt_c(r, "compute_bayesian_bootstrap_confidence_interval_basic"),                                                                  # bayes-basic ci (no p type)
        cell_bbt_p(r, "compute_bayesian_bootstrap_two_sided_pval_wald"), cell_bbt_c(r, "compute_bayesian_bootstrap_confidence_interval_wald"), # bayes-wald p + ci
        cell_bbt_p(r, "compute_bayesian_bootstrap_two_sided_pval_bca"), cell_bbt_c(r, "compute_bayesian_bootstrap_confidence_interval_bca"),   # bayes-bca p + ci
        cell_bbt_p(r, "compute_bayesian_bootstrap_two_sided_pval_studentized"), cell_bbt_c(r, "compute_bayesian_bootstrap_confidence_interval_studentized"), # bayes-stud p + ci
        # Param Boot (5): parametric-bootstrap-estimate (est/pval/ci) then LR (pval/ci)
        cell_pbest_estimate(r), cell_pbest_p(r), cell_pbest_c(r),
        cell_pb_p(r, "compute_lik_ratio_bootstrap_two_sided_pval"), cell_pb_c(r, "compute_lik_ratio_bootstrap_confidence_interval"),
        # Jackknife (3)
        cell_jack_estimate(r),
        cell_jack_p(r, "compute_jackknife_wald_two_sided_pval"), cell_jack_c(r, "compute_jackknife_wald_confidence_interval"),
        # Boot-Rand (9): pctile(pval+delta+ci), stud(pval+ci), symm-t(pval+ci), smooth(pval+ci)
        cell_brt_p(r, "compute_rand_bootstrap_two_sided_pval"), cell_brt_p(r, "compute_rand_bootstrap_two_sided_pval(delta=0.5)"), cell_brt_c(r, "compute_rand_bootstrap_confidence_interval"), # pctile
        cell_brt_p(r, "compute_rand_bootstrap_two_sided_pval_studentized"), cell_brt_c(r, "compute_rand_bootstrap_confidence_interval_studentized"),                                             # stud
        cell_brt_p(r, "compute_rand_bootstrap_two_sided_pval_symmetric-percentile-t"), cell_brt_c(r, "compute_rand_bootstrap_confidence_interval_symmetric-percentile-t"),                         # symm-t
        cell_brt_p(r, "compute_rand_bootstrap_two_sided_pval_smoothed"), cell_brt_c(r, "compute_rand_bootstrap_confidence_interval_smoothed"),                                                     # smooth
        # Exact other (2)
        cell_exact_p(r), cell_exact_c(r),
        # Randomization (5): vanilla pval + pval(δ) + ci, then custom pval + ci
        cell_rand_p(r, "compute_rand_two_sided_pval"), cell_rand_p(r, "compute_rand_two_sided_pval(delta=0.5)"), cell_rand_c(r, "compute_rand_confidence_interval"),
        cell_rand_custom_p(r, "compute_rand_two_sided_pval(custom)"), cell_rand_custom_c(r, "compute_rand_confidence_interval(custom)"),
        "</tr>\n")
    }
  }

  legend = '<div style="font-family:sans-serif;font-size:12px;margin-bottom:8px">
    <b>Legend:</b>
    <span style="background:#2e7d32;color:#fff;padding:2px 6px">✓ theoretically guaranteed</span>
    <span style="background:#a5d6a7;padding:2px 6px">100% observed estimable</span>
    <span style="background:#c8e6c9;padding:2px 4px">other</span><span style="background:#dce775;padding:2px 4px">varying</span><span style="background:#fff59d;padding:2px 4px">degrees</span><span style="background:#ffcc80;padding:2px 4px">of</span><span style="background:#ffb74d;padding:2px 4px">estimability</span>
    <span style="background:#ffcdd2;padding:2px 6px">SLOW skipped/lightly tested</span>
    <span style="background:#333333;color:#fff;padding:2px 6px">NTS not theoretically supported</span>
    <span style="background:#e0e0e0;padding:2px 6px">NI not implemented yet</span>
  </div>'

  reliability_note = '<p style="color:#555"><b>Note:</b> Completion rates are empirical stress-test rates, not mathematical guarantees.
  Low completion usually reflects finite-sample pathologies such as sparse data, separation, censoring, boundary estimates, or resampling degeneracy.
  A low rate does not by itself mean the estimator is invalid; it means the path is numerically fragile under these comprehensive-test settings.</p>'

  html = paste0('<!DOCTYPE html><html><head><meta charset="utf-8">
  <title>EDI Comprehensive Tests Coverage Audit</title>
  <style>
    body{font-family:sans-serif;font-size:13px;margin:16px}
    .scroll-x-top{overflow-x:auto;overflow-y:hidden;border:1px solid #ccc;border-bottom:none;position:sticky;top:0;background:#fff;z-index:20}
    .scroll-x-spacer{height:1px}
    .table-wrap{overflow-x:hidden;overflow-y:auto;max-height:85vh;border:1px solid #ccc}
    table{border-collapse:collapse;font-size:11px;border-spacing:0}
    td,th{border:1px solid #ccc;padding:2px 5px}
    .hdr0 th,.hdr1 th,.hdr2 th,.hdr3 th{position:sticky;background-clip:padding-box}
    .hdr0 th{padding:3px 6px;z-index:5;top:0;background:#263238;color:white;text-align:center;border-bottom:none}
    .hdr1 th{padding:3px 6px;z-index:4;top:0;background:#37474f;color:white;text-align:center;border-bottom:none}
    .hdr2 th{padding:2px 5px;z-index:3;top:0;background:#546e7a;color:white;text-align:center;border-bottom:none}
    .hdr3 th{padding:2px 5px;z-index:2;top:0;background:#607d8b;color:white;text-align:center}
    .hdr0 th[rowspan]{z-index:6;top:0;background:#263238;border-bottom:1px solid #ccc;vertical-align:bottom;text-align:left}
    td[style*="monospace"]{position:sticky;left:0;background:#fff;z-index:1;border-right:2px solid #999}
    .low-estimability{margin-top:8px;font-size:12px}
    .low-estimability th{background:#ef6c00;color:#fff;text-align:left;padding:4px 6px;position:static}
    .low-estimability td{position:static;padding:3px 6px}
    .low-estimability td[style*="monospace"]{position:static;left:auto;background:transparent;z-index:auto;border-right:1px solid #ccc}
  </style>
  <script>
    document.addEventListener("DOMContentLoaded", function() {
      var off = 0;
      ["hdr0","hdr1","hdr2","hdr3"].forEach(function(cls) {
        var tr = document.querySelector("tr." + cls);
        if (!tr) return;
        tr.querySelectorAll("th").forEach(function(th) { th.style.top = off + "px"; });
        off += tr.getBoundingClientRect().height;
      });
      var topScroll = document.querySelector(".scroll-x-top");
      var wrap = document.querySelector(".table-wrap");
      var table = wrap ? wrap.querySelector("table") : null;
      var spacer = document.querySelector(".scroll-x-spacer");
      if (topScroll && wrap && table && spacer) {
        var syncWidth = function() { spacer.style.width = table.scrollWidth + "px"; };
        syncWidth();
        window.addEventListener("resize", syncWidth);
        topScroll.addEventListener("scroll", function() { wrap.scrollLeft = topScroll.scrollLeft; });
        wrap.addEventListener("scroll", function() { topScroll.scrollLeft = wrap.scrollLeft; });
      }
    });
  </script>
  </head><body>
  <h2>EDI Comprehensive Tests Coverage Audit</h2>
  <p style="color:#555">Generated ', format(Sys.Date()), '. Green-to-yellow attempted cells use empirical explicit-nonestimability rates from main
  <code>comprehensive_tests_results*.csv</code> files in <code>package_tests/</code>; focused <code>_filtered_</code> runs are excluded.
  Dark green means the method is theoretically guaranteed under the comprehensive-test contract.
  rand/rci/brt columns: NTS for response types where randomization is not theoretically supported (incidence/ordinal: rand=NTS; count: rci/brt_ci=NTS).</p>
  ', reliability_note, legend, '<div class="scroll-x-top"><div class="scroll-x-spacer"></div></div><div class="table-wrap"><table>', hdr, body, '</table></div>', low_estimability_html(), '</body></html>')
  writeLines(html, outfile)
  invisible(outfile)
}

# Run from repo root: Rscript package_tests/path_audits_source.R
this_dir = tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) "package_tests")
html_from_audit(audit_classes, file.path(this_dir, "path_audits.html"))
