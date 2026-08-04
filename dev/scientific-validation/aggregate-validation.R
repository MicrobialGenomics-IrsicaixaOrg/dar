#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validation_root <- normalizePath(
  dirname(sub("^--file=", "", script_argument[[1L]])),
  mustWork = TRUE
)
options(dar.validation.root = validation_root)
source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "metrics.R"), local = FALSE)
source(file.path(validation_root, "R", "baseline.R"), local = FALSE)
source(file.path(validation_root, "R", "report.R"), local = FALSE)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: aggregate-validation.R INPUT_DIR OUTPUT_DIR", call. = FALSE)
}
output_dir <- normalizePath(arguments[[2L]], mustWork = FALSE)
outcome <- aggregate_validation_artifacts(
  arguments[[1L]], output_dir, render = TRUE, strict = TRUE,
  enforce_baseline = identical(Sys.getenv("VALIDATION_PROFILE"), "full")
)
message("Aggregate scientific validation report: ",
        file.path(output_dir, "scientific-validation.html"))
quit(status = outcome$status, save = "no")
