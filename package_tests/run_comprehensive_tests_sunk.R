zz = file("/tmp/comprehensive_count_nc2_clean_20260712.log", "wt")
sink(zz)
sink(zz, type = "message")
source("package_tests/comprehensive_tests.R")
