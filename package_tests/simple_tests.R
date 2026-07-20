library(EDI)

# Settings
n = 100
p = 5
data_type = "linear"
design_cls = DesignSeqOneByOneKK21
r = 50 # number of bootstrap/randomization iterations for speed

# 1. Generate Data
dataset = generate_covariate_dataset(n, p, data_type)
X = dataset$X
y_cont = dataset$y_cont

run_tests_for_response = function(response_type, inference_classes, model_formula = NULL) {
  cat("\n\n###########################################################################\n")
  cat("##### response_type =", response_type, 
      if (!is.null(model_formula)) paste0(" (formula: ", deparse(model_formula), ")") else "", "\n")
  cat("###########################################################################\n")
  
  y_base = transform_cont_y_based_on_response_type(y_cont, response_type)
  
  # Initialize Design
  des = design_cls$new(response_type = response_type, n = n, design_formula = if (is.null(model_formula)) ~ . else model_formula)
  
  # Fill Design (Sequential)
  betaT = 1
  sd_noise = 0.1
  
  for (i in 1:n) {
    w_i = des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
    bt = if (w_i == 1) betaT else 0
    eps = rnorm(1, 0, sd_noise)
    
    y_i = switch(response_type,
      continuous = y_base[i] + bt + eps,
      incidence  = {
        p_b = if (is.finite(y_base[i]) && y_base[i] >= 0 && y_base[i] <= 1) y_base[i] else plogis(y_base[i])
        p_b = pmin(0.95, pmax(0.05, p_b))
        as.numeric(rbinom(1, 1, plogis(qlogis(p_b) + bt + eps)))
      },
      proportion = pmin(1 - 1e-9, pmax(1e-9, y_base[i] + bt + eps)),
      count = {
        lam = pmax(.Machine$double.eps, y_base[i] * exp(bt + eps))
        as.numeric(rpois(1, lambda = lam))
      },
      survival = pmax(.Machine$double.eps, y_base[i] * exp(bt + eps)),
      ordinal  = as.integer(max(1, round(y_base[i] + bt + eps))),
      stop("Unknown response_type")
    )
    des$add_one_subject_response(i, y_i, dead = 1)
  }
  
  for (inf_info in inference_classes) {
    inf_cls = if (is.list(inf_info)) inf_info[[1]] else inf_info
    inf_args = if (is.list(inf_info) && length(inf_info) > 1) inf_info[-1] else list()
    
    # Filter inf_args to only include those in the initialize method
    init_params = names(formals(inf_cls$public_methods$initialize))
    
    # If model_formula is provided, ensure it's passed to inference too if it supports it
    if (!is.null(model_formula)) {
      if ("model_formula" %in% init_params) {
        inf_args$model_formula = model_formula
      }
    }
    
    inf_args_filtered = inf_args[names(inf_args) %in% init_params]
    
    cls_name = if (is.character(inf_cls)) inf_cls else inf_cls$classname
    cat("\n--- Inference:", cls_name, if (length(inf_args_filtered)) paste0("(", paste(names(inf_args_filtered), inf_args_filtered, sep="=", collapse=", "), ")") else "", "---\n")

    inf = do.call(inf_cls$new, c(list(des_obj = des), inf_args_filtered))
    
    # 1. Treatment Estimate
    tryCatch({
      est = inf$compute_estimate()
      cat("  Estimate:", est, "\n")
    }, error = function(e) cat("  Estimate Error:", e$message, "\n"))
    
    # 2. Asymptotic Inference
    if (inherits(inf, "InferenceAsymp")) {
      tryCatch({
        cat("  Asymp P-val:", inf$compute_asymp_two_sided_pval(), "\n")
        cat("  Asymp CI:", paste(inf$compute_asymp_confidence_interval(), collapse = ", "), "\n")
      }, error = function(e) cat("  Asymp Error:", e$message, "\n"))
    }
    
    # 3. Bootstrap Inference
    if (inherits(inf, "InferenceBoot")) {
      tryCatch({
        cat("  Boot P-val:", inf$compute_bootstrap_two_sided_pval(B = r), "\n")
        cat("  Boot CI:", paste(inf$compute_bootstrap_confidence_interval(B = r), collapse = ", "), "\n")
      }, error = function(e) cat("  Boot Error:", e$message, "\n"))
    }
    
    # 4. Randomization Inference
    if (inherits(inf, "InferenceRand")) {
      tryCatch({
        cat("  Rand P-val:", inf$compute_rand_two_sided_pval(r = r), "\n")
        if (inherits(inf, "InferenceRandCI")) {
          cat("  Rand CI:", paste(inf$compute_rand_confidence_interval(r = r), collapse = ", "), "\n")
        }
      }, error = function(e) cat("  Rand Error:", e$message, "\n"))
    }
  }
}

##### response_type = continuous (standard)
run_tests_for_response("continuous", list(
  InferenceAllSimpleMeanDiff,
  list(InferenceContinOLS, model_formula = ~ 1),
  list(InferenceContinOLS, model_formula = ~ .),
  list(InferenceContinLin, model_formula = ~ .),
  list(InferenceContinQuantileRegr, model_formula = ~ .),
  list(InferenceContinRobustRegr, model_formula = ~ .),
  list(InferenceContinKKRobustRegrIVWC, model_formula = ~ .),
  list(InferenceContinKKRobustRegrOneLik, model_formula = ~ .),
  list(InferenceContinKKGLMM, model_formula = ~ .)
))

##### response_type = continuous (model_formula = ~ .)
run_tests_for_response("continuous", list(
  list(InferenceContinOLS)
), model_formula = ~ .)

##### response_type = continuous (model_formula = ~ 1)
run_tests_for_response("continuous", list(
  list(InferenceContinOLS)
), model_formula = ~ 1)

##### response_type = incidence
run_tests_for_response("incidence", list(
  list(InferenceIncidLogRegr, model_formula = ~ 1),
  list(InferenceIncidLogRegr, model_formula = ~ .),
  list(InferenceIncidLogBinomial, model_formula = ~ .),
  list(InferenceIncidModifiedPoisson, model_formula = ~ .),
  list(InferenceIncidRiskDiff, model_formula = ~ .),
  list(InferenceIncidBinomialIdentityRiskDiff, model_formula = ~ .),
  list(InferenceIncidGCompRiskDiff, model_formula = ~ .),
  list(InferenceIncidGCompRiskRatio, model_formula = ~ .),
  list(InferenceIncidKKCondLogitIVWC, model_formula = ~ .),
  list(InferenceIncidKKCondLogitOneLik, model_formula = ~ .),
  list(InferenceIncidKKGEE, model_formula = ~ .),
  list(InferenceIncidKKCondLogitPlusGLMMOneLik, model_formula = ~ .)
))

##### response_type = proportion
run_tests_for_response("proportion", list(
  InferenceAllSimpleMeanDiff,
  list(InferencePropBetaRegr, model_formula = ~ .),
  list(InferencePropFractionalLogit, model_formula = ~ .),
  list(InferencePropQuantileRegr, model_formula = ~ .),
  list(InferencePropGCompMeanDiff, model_formula = ~ .),
  list(InferencePropKKGEE, model_formula = ~ .),
  list(InferencePropKKGLMM, model_formula = ~ .)
))

##### response_type = count
run_tests_for_response("count", list(
  list(InferenceCountPoisson, model_formula = ~ .),
  list(InferenceCountNegBin, model_formula = ~ .)
))

##### response_type = survival
run_tests_for_response("survival", list(
  InferenceSurvivalLogRank,
  list(InferenceSurvivalCoxPHRegr, model_formula = ~ .),
  list(InferenceSurvivalKKLWACoxPHIVWC, model_formula = ~ .)
))

##### response_type = ordinal
run_tests_for_response("ordinal", list(
  InferenceOrdinalJonckheereTerpstraTest
))
