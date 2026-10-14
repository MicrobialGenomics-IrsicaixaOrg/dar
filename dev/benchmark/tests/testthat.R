#!/usr/bin/env Rscript
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(dirname(sub("^--file=", "", script_argument[[1L]])), mustWork = TRUE)
Sys.setenv(TESTTHAT_PARALLEL = "FALSE")
testthat::test_dir(file.path(root, "testthat"), reporter = "summary",
                   stop_on_failure = TRUE)
