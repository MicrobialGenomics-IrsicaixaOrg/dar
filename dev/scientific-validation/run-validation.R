#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Unable to locate the validation script.", call. = FALSE)
}
validation_root <- normalizePath(
  dirname(sub("^--file=", "", script_argument[[1L]])),
  mustWork = TRUE
)
repository_root <- normalizePath(file.path(validation_root, "..", ".."), mustWork = TRUE)
options(dar.validation.root = validation_root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repository_root, quiet = TRUE)
} else {
  library(dar)
}

source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "simulate.R"), local = FALSE)
source(file.path(validation_root, "R", "engines.R"), local = FALSE)
source(file.path(validation_root, "R", "metrics.R"), local = FALSE)
source(file.path(validation_root, "R", "cli.R"), local = FALSE)

arguments <- parse_validation_arguments(commandArgs(trailingOnly = TRUE))
outcome <- run_scientific_validation(
  profile_name = arguments$profile,
  engines = arguments$engine,
  scenarios = arguments$scenario,
  output_dir = arguments$output_dir,
  base_seed = arguments$base_seed,
  strict = arguments$strict
)

message("Scientific validation artifacts: ", normalizePath(arguments$output_dir))
message("Direction reversal self-test: ", if (outcome$direction_self_test) "PASS" else "FAIL")
if (nrow(outcome$runs)) {
  print(table(outcome$runs$engine, outcome$runs$status))
}
quit(status = outcome$status, save = "no")
