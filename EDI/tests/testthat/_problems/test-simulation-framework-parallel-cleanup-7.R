# Extracted from test-simulation-framework-parallel-cleanup.R:7

# prequel ----------------------------------------------------------------------
library(testthat)
library(EDI)

# test -------------------------------------------------------------------------
set_num_cores(2L)
initial_cores = get_num_cores()
