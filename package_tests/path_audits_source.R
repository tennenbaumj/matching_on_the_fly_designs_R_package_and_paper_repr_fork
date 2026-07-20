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

audit_classes = list(

  # ── GLOBAL ──────────────────────────────────────────────────────────────────
  list(name="InferenceAllSimpleMeanDiff",           section="Global",     resp="all",   kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, notes="rci for cont only (non-cont AllSimpleMeanDiff skipped)"),
  list(name="InferenceAllSimpleMeanDiffPooledVar",  section="Global",     resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, notes="inherits AllSimpleMeanDiff; tested cont only"),
  list(name="InferenceAllSimpleWilcox",             section="Global",     resp="all",   kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=NA,    jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA, notes="boot enabled; Bayesian bootstrap structurally unsupported (fractional weights undefined for rank statistic); BRT CI smoothed avg 94s slow"),

  # ── CONTINUOUS KK ────────────────────────────────────────────────────────────
  list(name="InferenceContinKKGLMM",               section="Continuous", resp="cont",  kk=TRUE,  types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp",  pboot=TRUE,  notes="InferenceParamBootstrap; all resampling too slow"),
  list(name="InferenceContinKKOLSOneLik",          section="Continuous", resp="cont",  kk=TRUE,  types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  notes="KKPassThroughCompound → ParamBootstrap; explicit pboot=TRUE"),
  list(name="InferenceContinKKRobustRegrOneLik",   section="Continuous", resp="cont",  kk=TRUE,  types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA,    notes="KKPassThroughCompoundNoParamBootstrap; wald only"),
  list(name="InferenceContinKKQuantileRegrOneLik", section="Continuous", resp="cont",  kk=TRUE,  types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA,    notes="AbstractQuantileRandCI → KKPassThroughCompoundNoParamBootstrap"),

  # ── CONTINUOUS non-KK ────────────────────────────────────────────────────────
  list(name="InferenceContinLin",                  section="Continuous", resp="cont",  kk=FALSE, types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceContinOLS",                  section="Continuous", resp="cont",  kk=FALSE, types="full",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceContinRobustRegr",           section="Continuous", resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp",  pboot=NA,    skip_brt="ci", notes="InferenceAsymp"),
  list(name="InferenceContinQuantileRegr",         section="Continuous", resp="cont",  kk=FALSE, types="wald",          skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp",  pboot=NA,    notes="InferenceAsymp"),

  # ── INCIDENCE KK ─────────────────────────────────────────────────────────────
  list(name="InferenceIncidKKCondLogitOneLik",          section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap directly; explicit TRUE"),
  list(name="InferenceIncidKKCondLogitPlusGLMMOneLik",  section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, skip_bbt_ci=TRUE, notes="AbstractKKCondLogitPlusGLMM → InferenceParamBootstrap; supports_lik_ratio_param_bootstrap=FALSE (no simulate_under_lik_null); rand delta=0.5 pval avg 205s slow"),
  list(name="InferenceIncidKKGEE",                     section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidKKNewcombeRiskDiff",         section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="KKPassThroughCompoundNoParamBootstrap"),
  list(name="InferenceIncidKKGCompRiskDiff",            section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, notes="AbstractKKMarginalIncid → ParamBootstrap; supports_lik_ratio_param_bootstrap()=FALSE (no simulate_under_lik_null)"),
  list(name="InferenceIncidKKGCompRiskRatio",           section="Incidence", resp="incid", kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=FALSE, notes="AbstractKKMarginalIncid → ParamBootstrap; supports_lik_ratio_param_bootstrap()=FALSE (no simulate_under_lik_null)"),
  list(name="InferenceIncidKKModifiedPoisson",          section="Incidence", resp="incid", kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AbstractKKModifiedPoisson → AbstractKKMarginalIncid → ParamBootstrap; simulate_under_lik_null: Poisson draw"),

  # ── INCIDENCE non-KK ──────────────────────────────────────────────────────────
  list(name="InferenceIncidExactFisher",               section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact → InferenceJackknife; runs exact pval+CI only"),
  list(name="InferenceIncidExactBinomial",             section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact; KK14/FixedBinaryMatch"),
  list(name="InferenceIncidenceExactZhang",            section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=TRUE,  skip_ci=TRUE,  skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="",  skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, exact_c=TRUE, notes="InferenceExact"),
  list(name="InferenceIncidWald",                     section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="inherits AllSimpleMeanDiff (ParamBootstrap); explicit FALSE"),
  list(name="InferenceIncidCMH",                      section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="ParamBootstrap; explicit FALSE"),
  list(name="InferenceIncidExtendedRobins",            section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="ParamBootstrap; explicit FALSE"),
  list(name="InferenceIncidLogRegr",                  section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferenceIncidProbitRegr",               section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferenceIncidMiettinenNurminenRiskDiff",section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidNewcombeRiskDiff",         section="Incidence", resp="incid", kk=FALSE, types="none",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp; get_supported_testing_types_impl=character(0); no direct tests"),
  list(name="InferenceIncidRiskDiff",                 section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="AsympLikStdModCacheNoParamBootstrap; supports_likelihood_tests=FALSE"),
  list(name="InferenceIncidGCompRiskDiff",            section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidGCompRiskRatio",           section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceIncidModifiedPoisson",          section="Incidence", resp="incid", kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA, notes="ParamBootstrap; explicit FALSE"),
  list(name="InferenceIncidLogBinomial",              section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferenceIncidBinomialIdentityRiskDiff", section="Incidence", resp="incid", kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="i", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache"),

  # ── PROPORTION ───────────────────────────────────────────────────────────────
  list(name="InferencePropKKGEE",                   section="Proportion", resp="prop",  kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferencePropKKGLMM",                  section="Proportion", resp="prop",  kk=TRUE,  types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=TRUE,  notes="AbstractKKCondLogitPlusGLMM → InferenceParamBootstrap; simulate_under_lik_null added"),
  list(name="InferencePropKKQuantileRegrOneLik",    section="Proportion", resp="prop",  kk=TRUE,  types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=NA,    notes="AbstractKKQuantileRegrOneLik → KKPassThroughCompoundNoParamBootstrap; BRT CI studentized avg 34s + smoothed avg 33s slow"),
  list(name="InferencePropQuantileRegr",            section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=NA,    notes="InferenceAsymp; logit-scale quantile regression, mirrors InferenceContinQuantileRegr"),
  list(name="InferencePropBetaRegr",                section="Proportion", resp="prop",  kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE,  skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="cp", pboot=TRUE,  notes="AsympLikStdModCache"),
  list(name="InferencePropZeroOneInflatedBetaRegr", section="Proportion", resp="prop",  kk=FALSE, types="full",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp", pboot=TRUE,  notes="AsympLikStdModCache; pboot=TRUE (separate from bootstrap)"),
  list(name="InferencePropGCompMeanDiff",           section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=TRUE,  rand_resp="csp", skip_rpv=TRUE,  skip_rci=TRUE,  rci_resp="cp", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferencePropFractionalLogit",         section="Proportion", resp="prop",  kk=FALSE, types="wald",   skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="cp", pboot=NA,    notes="AsympLikStdModCacheNoParamBootstrap; supports_likelihood_tests=FALSE"),

  # ── COUNT ────────────────────────────────────────────────────────────────────
  list(name="InferenceCountPoissonKKGEE",           section="Count",      resp="count", kk=TRUE,  types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA,    notes="InferenceAsymp; rand/rci N/A for count"),
  list(name="InferenceCountKKGLMM",                 section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=FALSE, skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; supports_lik_ratio_param_bootstrap=TRUE when use_rcpp; NO jackknife"),
  list(name="InferenceCountKKHurdlePoissonOneLik",  section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; simulate_under_lik_null: ZTP+GLMM matched + Poisson reservoir; in skip_ci_rand; rand N/A count"),
  list(name="InferenceCountKKCondPoissonOneLik",    section="Count",      resp="count", kk=TRUE,  types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="InferenceParamBootstrap; simulate_under_lik_null: cond-Poisson pairs Binomial + reservoir Poisson; in skip_ci_rand; rand N/A count"),
  list(name="InferenceCountPoisson",                section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="CountLikelihood → AsympLikStdModCache"),
  list(name="InferenceCountRobustPoisson",          section="Count",      resp="count", kk=FALSE, types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA, notes="CountCompositeLikelihood; explicit FALSE"),
  list(name="InferenceCountQuasiPoisson",           section="Count",      resp="count", kk=FALSE, types="wald",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=NA, notes="CountCompositeLikelihood; explicit FALSE"),
  list(name="InferenceCountNegBin",                 section="Count",      resp="count", kk=FALSE, types="full",skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE, notes="InferenceCountLikelihood → InferenceParamBootstrap; simulate_under_lik_null: rnbinom draw"),
  list(name="InferenceCountZeroInflatedPoisson",    section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; use_rcpp; 'Zero-Inflated Poisson' in pboot list; in skip_ci_rand"),
  list(name="InferenceCountZeroInflatedNegBin",     section="Count",      resp="count", kk=FALSE, types="full",  skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; use_rcpp; 'Zero-Inflated Negative Binomial' in pboot list"),
  list(name="InferenceCountHurdlePoisson",          section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE,  notes="ZAAbstract; 'Hurdle Poisson' in pboot list; in skip_ci_rand"),
  list(name="InferenceCountHurdleNegBin",           section="Count",      resp="count", kk=FALSE, types="full",           skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="count",  skip_rpv=FALSE, skip_rci=TRUE,  rci_resp="", pboot=TRUE, notes="InferenceCountLikelihood → InferenceParamBootstrap; simulate_under_lik_null implemented; in skip_ci_rand"),

  # ── SURVIVAL ─────────────────────────────────────────────────────────────────
  list(name="InferenceSurvivalGehanWilcox",            section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    skip_brt="ci", notes="InferenceAsymp; all BRT CI types avg 33-39s slow"),
  list(name="InferenceSurvivalLogRank",                section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalRestrictedMeanDiff",     section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalKMDiff",                 section="Survival",  resp="surv", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceSurvivalWeibullRegr",            section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  skip_brt="ci", notes="AsympLikStdModCache → ParamBootstrap; all BRT CI types avg 42-518s slow"),
  list(name="InferenceSurvivalDepCensTransformRegr",   section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AsympLikStdModCache; simulate_under_lik_null added"),
  list(name="InferenceSurvivalCoxPHRegr",              section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AsympLikStdModCache; pboot=use_rcpp (default)"),
  list(name="InferenceSurvivalStratCoxPHRegr",         section="Survival",  resp="surv", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="InferenceParamBootstrap directly; pboot=use_rcpp"),
  list(name="InferenceSurvivalKKClaytonCopulaOneLik",  section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=FALSE, skip_rand=FALSE, rand_resp="csp", skip_rpv=TRUE,  skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  skip_pboot_ci=TRUE, notes="InferenceParamBootstrap; pboot=TRUE even though skip_boot=TRUE (separate test family)"),
  list(name="InferenceSurvivalKKLWACoxPHOneLik",       section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AbstractKKLWACoxOneLik → ParamBootstrap; explicit TRUE"),
  list(name="InferenceSurvivalKKStratCoxPHOneLik",     section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="InferenceParamBootstrap; explicit TRUE"),
  list(name="InferenceSurvivalKKWeibullFrailtyOneLik", section="Survival",  resp="surv", kk=TRUE,  types="full", skip_asymp=FALSE, skip_ci=FALSE,   skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE, skip_rand=FALSE, rand_resp="csp", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="surv", pboot=TRUE,  notes="AbstractKKWeibullFrailtyOneLik → ParamBootstrap; pboot=use_rcpp; skip_boot=TRUE"),
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
  list(name="InferenceOrdinalAdjCatLogitRegr",          section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache; simulate_under_lik_null added"),
  list(name="InferenceOrdinalStereotypeLogitRegr",      section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache; simulate_under_lik_null added"),
  list(name="InferenceOrdinalContRatioRegr",            section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache; simulate_under_lik_null added; BRT pval smoothed avg 66s slow"),
  list(name="InferenceOrdinalCloglogRegr",              section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE even with skip_boot"),
  list(name="InferenceOrdinalCauchitRegr",              section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE"),
  list(name="InferenceOrdinalOrderedProbitRegr",        section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap; pboot=TRUE"),
  list(name="InferenceOrdinalPropOddsRegr",             section="Ordinal", resp="ord", kk=FALSE, types="full", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=TRUE,  notes="AsympLikStdModCache → ParamBootstrap"),
  list(name="InferenceOrdinalPartialProportionalOddsRegr",section="Ordinal",resp="ord",kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceOrdinalGCompMeanDiff",            section="Ordinal", resp="ord", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=TRUE,  skip_bbt=TRUE,  jack=TRUE,  skip_rand=TRUE,  rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp"),
  list(name="InferenceOrdinalJonckheereTerpstraTest",   section="Ordinal", resp="ord", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    exact_p=TRUE, notes="InferenceAsymp; SPECIAL: only exact pval + estimate called"),
  list(name="InferenceOrdinalRidit",                   section="Ordinal", resp="ord", kk=FALSE, types="wald", skip_asymp=FALSE, skip_ci=FALSE, skip_boot=FALSE, skip_bbt=FALSE, jack=TRUE,  skip_rand=FALSE, rand_resp="ordinal", skip_rpv=FALSE, skip_rci=FALSE, rci_resp="", pboot=NA,    notes="InferenceAsymp")
)

# ── HTML renderer ─────────────────────────────────────────────────────────────
html_from_audit = function(classes, outfile = "path_audits.html") {
  GREEN = "#c8e6c9"; RED = "#ffcdd2"; YELLOW = "#fff9c4"
  ok  = function(s="✓")    sprintf('<td style="background:%s;text-align:center">%s</td>', GREEN,  s)
  bad = function(s="SLOW") sprintf('<td style="background:%s;text-align:center">%s</td>', YELLOW, s)
  na  = function(s="N/A")  sprintf('<td style="background:%s;text-align:center">%s</td>', RED,    s)

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
  cell_ap = function(r, ttype) if (type_ok(r, ttype)) ok() else na()
  cell_ac = function(r, ttype) {
    if (!type_ok(r, ttype)) return(na())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Nonparam bootstrap (classical) ─────────────────────────────────────────
  cell_bp = function(r) if (isTRUE(r$skip_boot)) bad() else ok()
  cell_bc = function(r) {
    if (isTRUE(r$skip_boot)) return(bad())
    if (isTRUE(r$skip_boot_ci)) return(bad())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Bayesian bootstrap ──────────────────────────────────────────────────────
  cell_bbt_p = function(r) {
    if (is.na(r$skip_bbt)) return(na())
    if (isTRUE(r$skip_bbt)) bad() else ok()
  }
  cell_bbt_c = function(r) {
    if (is.na(r$skip_bbt)) return(na())
    if (isTRUE(r$skip_bbt)) return(bad())
    if (isTRUE(r$skip_bbt_ci)) return(bad())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Studentized bootstrap (separate skip from general boot) ──────────────────
  cell_bst_p = function(r) {
    if (isTRUE(r$skip_boot)) return(bad())
    if (isTRUE(r$skip_stud)) return(bad())
    ok()
  }
  cell_bst_c = function(r) {
    if (isTRUE(r$skip_boot)) return(bad())
    if (isTRUE(r$skip_stud)) return(bad())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Parametric bootstrap ────────────────────────────────────────────────────
  cell_pb_p = function(r) {
    if (is.na(r$pboot)) return(na())
    if (isTRUE(r$pboot)) ok() else bad()
  }
  cell_pb_c = function(r) {
    if (is.na(r$pboot)) return(na())
    if (!isTRUE(r$pboot)) return(bad())
    if (isTRUE(r$skip_pboot_ci)) return(bad())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Jackknife ───────────────────────────────────────────────────────────────
  cell_jack_p = function(r) if (isTRUE(r$jack)) ok() else bad()
  cell_jack_c = function(r) {
    if (!isTRUE(r$jack)) return(bad())
    if (isTRUE(r$skip_ci)) return(bad()) else ok()
  }

  # ── Randomization ───────────────────────────────────────────────────────────
  cell_rand_p = function(r) {
    if (!nchar(r$rand_resp)) return(na())
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(bad())
    ok()
  }
  cell_rand_c = function(r) {
    if (!nchar(r$rci_resp)) return(na())
    if (isTRUE(r$skip_rci) || isTRUE(r$skip_rand_ci)) return(bad()) else ok()
  }
  cell_rand_custom_p = function(r) {  # custom statistic: same conditions as vanilla
    if (!nchar(r$rand_resp)) return(na())
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(bad())
    ok()
  }
  cell_rand_custom_c = function(r) {  # custom CI: continuous, proportion, survival
    if (!r$resp %in% c("cont", "prop", "surv")) return(na())
    if (isTRUE(r$skip_rci)) return(bad()) else ok()
  }

  # ── Bootstrap-randomization (BRT) — all types share same per-class logic ───
  # skip_brt="ci" means pval runs but CI is too slow (e.g. ContinRobustRegr)
  cell_brt_p = function(r) {
    if (!nchar(r$rand_resp)) return(na())
    if (isTRUE(r$skip_rand) || isTRUE(r$skip_rpv)) return(bad())
    if (!identical(r$skip_brt, "ci") && isTRUE(r$skip_boot)) return(bad())
    ok()
  }
  cell_brt_c = function(r) {
    if (!nchar(r$rci_resp)) return(na())
    if (identical(r$skip_brt, "ci")) return(bad())
    if (isTRUE(r$skip_boot) || isTRUE(r$skip_rand)) return(bad())
    if (isTRUE(r$skip_rand_ci)) return(bad())
    if (isTRUE(r$skip_ci) || isTRUE(r$skip_rci)) return(bad()) else ok()
  }

  # ── Exact inference (Fisher/Binomial/Zhang: pval+CI; Jonckheer: pval only) ───
  cell_exact_p = function(r) if (isTRUE(r$exact_p)) ok() else na()
  cell_exact_c = function(r) if (isTRUE(r$exact_c)) ok() else na()

  # ── Four frozen header rows ──────────────────────────────────────────────────
  # Data cols: 8 (model-based) + 2 (exact-other) + 5 (exact-rand) + 18 (npboot, incl. 10 for the 6 bayes flavors)
  #            + 2 (pboot) + 2 (jack) + 9 (brt) = 46
  NCOL = 47  # 1 (class) + 46 (data)
  hdr = paste0(
    # Row 0: meta-categories — Asymptotic (39) | Exact (7)
    '<tr class="hdr0">',
      '<th rowspan="4" style="text-align:left">Inference Type</th>',
      '<th colspan="39">Asymptotic</th>',
      '<th colspan="7">Exact</th>',
    '</tr>',
    # Row 1: categories
    '<tr class="hdr1">',
      # Under Asymptotic (39)
      '<th colspan="8">Model-Based</th>',
      '<th colspan="18">Nonparam Boot</th>',
      '<th colspan="2">Param Boot</th>',
      '<th colspan="2">Jackknife</th>',
      '<th colspan="9">Boot-Rand</th>',
      # Under Exact (7)
      '<th colspan="2">other</th>',
      '<th colspan="5">Randomization</th>',
    '</tr>',
    # Row 2: test types
    '<tr class="hdr2">',
      # Model-Based: 4 types × colspan=2
      '<th colspan="2">wald</th>',
      '<th colspan="2">score</th>',
      '<th colspan="2">likrat</th>',
      '<th colspan="2">grad</th>',
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
      # Param Boot
      '<th colspan="2">likrat</th>',
      # Jackknife (no named sub-type)
      '<th colspan="2"></th>',
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
      # Model-Based (8)
      '<th>pval</th><th>ci</th>',
      '<th>pval</th><th>ci</th>',
      '<th>pval</th><th>ci</th>',
      '<th>pval</th><th>ci</th>',
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
      # Param Boot (2)
      '<th>pval</th><th>ci</th>',
      # Jackknife (2)
      '<th>pval</th><th>ci</th>',
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
        # Model-Based (8)
        cell_ap(r,"wald"),  cell_ac(r,"wald"),
        cell_ap(r,"score"), cell_ac(r,"score"),
        cell_ap(r,"lr"),    cell_ac(r,"lr"),
        cell_ap(r,"grad"),  cell_ac(r,"grad"),
        # Nonparam Boot (18)
        cell_bp(r), cell_bc(r),        # pctile p + ci
        cell_bp(r),                     # symm p (no ci type)
        cell_bc(r),                     # basic ci (no p type)
        cell_bp(r), cell_bc(r),        # bca p + ci
        cell_bst_p(r), cell_bst_c(r),  # stud p + ci
        cell_bbt_p(r), cell_bbt_c(r),  # bayes-pctile p + ci
        cell_bbt_p(r),                  # bayes-symm p (no ci type)
        cell_bbt_c(r),                  # bayes-basic ci (no p type)
        cell_bbt_p(r), cell_bbt_c(r),  # bayes-wald p + ci
        cell_bbt_p(r), cell_bbt_c(r),  # bayes-bca p + ci
        cell_bbt_p(r), cell_bbt_c(r),  # bayes-stud p + ci
        # Param Boot (2)
        cell_pb_p(r), cell_pb_c(r),
        # Jackknife (2)
        cell_jack_p(r), cell_jack_c(r),
        # Boot-Rand (9): pctile(pval+delta+ci), stud(pval+ci), symm-t(pval+ci), smooth(pval+ci)
        cell_brt_p(r), cell_brt_p(r), cell_brt_c(r),  # pctile
        cell_brt_p(r), cell_brt_c(r),                  # stud
        cell_brt_p(r), cell_brt_c(r),                  # symm-t
        cell_brt_p(r), cell_brt_c(r),                  # smooth
        # Exact other (2)
        cell_exact_p(r), cell_exact_c(r),
        # Randomization (5): vanilla pval + pval(δ) + ci, then custom pval + ci
        cell_rand_p(r), cell_rand_p(r), cell_rand_c(r),
        cell_rand_custom_p(r), cell_rand_custom_c(r),
        "</tr>\n")
    }
  }

  legend = '<div style="font-family:sans-serif;font-size:12px;margin-bottom:8px">
    <b>Legend:</b>
    <span style="background:#c8e6c9;padding:2px 6px">✓ runs</span>
    <span style="background:#fff9c4;padding:2px 6px">SLOW too slow to test (skipped)</span>
    <span style="background:#ffcdd2;padding:2px 6px">N/A not supported by model</span>
  </div>'

  html = paste0('<!DOCTYPE html><html><head><meta charset="utf-8">
  <title>EDI Comprehensive Tests Coverage Audit</title>
  <style>
    body{font-family:sans-serif;font-size:13px;margin:16px}
    .table-wrap{overflow-x:auto;overflow-y:auto;max-height:85vh;border:1px solid #ccc}
    table{border-collapse:collapse;font-size:11px;border-spacing:0}
    td,th{border:1px solid #ccc;padding:2px 5px}
    .hdr0 th,.hdr1 th,.hdr2 th,.hdr3 th{position:sticky;background-clip:padding-box}
    .hdr0 th{padding:3px 6px;z-index:5;top:0;background:#263238;color:white;text-align:center;border-bottom:none}
    .hdr1 th{padding:3px 6px;z-index:4;top:21px;background:#37474f;color:white;text-align:center;border-bottom:none}
    .hdr2 th{padding:2px 5px;z-index:3;top:42px;background:#546e7a;color:white;text-align:center;border-bottom:none}
    .hdr3 th{padding:2px 5px;z-index:2;top:63px;background:#607d8b;color:white;text-align:center}
    .hdr0 th[rowspan]{z-index:6;top:0;background:#263238;border-bottom:1px solid #ccc;vertical-align:bottom;text-align:left}
    td[style*="monospace"]{position:sticky;left:0;background:#fff;z-index:1;border-right:2px solid #999}
  </style></head><body>
  <h2>EDI Comprehensive Tests Coverage Audit</h2>
  <p style="color:#555">Generated 2026-07-20. rand/rci/brt columns: N/A for response types where randomization is not applicable
  (incidence/ordinal: rand=N/A; count: rci/brt_ci=N/A).</p>
  ', legend, '<div class="table-wrap"><table>', hdr, body, '</table></div></body></html>')
  writeLines(html, outfile)
  invisible(outfile)
}

# Run from repo root: Rscript package_tests/path_audits_source.R
this_dir = tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) "package_tests")
html_from_audit(audit_classes, file.path(this_dir, "path_audits.html"))
