suppressPackageStartupMessages(library(EDI))
suppressPackageStartupMessages(library(data.table))

Nrep_W = 10003L   # Monte Carlo w-assignment draws per cell

sim = SimulationFramework$new(
        Nrep_W                        = Nrep_W,
        num_cores                     = 46L,
        results_filename              = sprintf("simulations/cmh_exact_sims_plus_greedy_results_high_signal_Nrep_%d.csv.bz2", Nrep_W),
        continue_from_last_result_row = TRUE,
        response_type                 = "incidence",
        n                             = c(64L, 128L, 256L),
        p                             = c(1L, 2L, 5L, 10L),
        norm_sq_beta_vec              = 2.5, #default is 1; we're increasing it here to give the x's more signal
        random_X_draws                = FALSE,
        seed                          = 1984,
        betaT                         = c(0, 0.5), # 0 → size / type-I error;  >0 → power / coverage
        alpha                         = 0.05,
        cond_exp_func_model           = c("linear"), #the "nonlinear" cond exp model is not as interesting for this simulation
        design_classes_and_params     = list(
                                          DesignFixediBCRD,
                                          DesignFixedBinaryMatch,
                                          DesignFixedOptimalBlocks =                list(B = 2),
                                          DesignFixedOptimalBlocks =                list(B = 4),
                                          DesignFixedOptimalBlocks =                list(B = 8),
                                          DesignFixedOptimalBlocks =                list(B = 16),
                                          DesignFixedOptimalBlocks =                list(B = 32),
                                          DesignFixedOptimalBlocks =                list(B = 64),
                                          DesignFixedOptimalBlocks =                list(B = 128),
                                          DesignFixedBlocking =                     list(B_target = 2,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 4,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 8,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 16, exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 32, exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 64,  exact_num_blocks = TRUE),
                                          DesignFixedBlocking =                     list(B_target = 128, exact_num_blocks = TRUE),
                                          DesignFixedGreedy =                       list(objective = "mahal_dist"),
                                          DesignFixedMatchingGreedyPairSwitching =  list(objective = "mahal_dist"),
                                          DesignFixedRerandomization =              list(prop_acceptable = 0.01) 
                                        ),
        inference_classes_and_params  = list(
                                          InferenceIncidWald,
                                          InferenceIncidCMH,
                                          InferenceIncidExtendedRobins 
                                          #InferenceIncidExactBinomial
                                        ),
        inference_types_and_params    = list(
                                          asymp_ci   = list(),
                                          asymp_pval = list()
                                          # exact_ci   = list(),
                                          # exact_pval = list()
                                        ),
        keep_all_intermediate_data    = FALSE,
        stop_on_error                 = FALSE
      )

suppressWarnings(sim$run())