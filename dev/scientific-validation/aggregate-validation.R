#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validation_root <- normalizePath(
  dirname(sub("^--file=", "", script_argument[[1L]])),
  mustWork = TRUE
)
options(dar.validation.root = validation_root)
source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "report.R"), local = FALSE)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: aggregate-validation.R INPUT_DIR OUTPUT_DIR", call. = FALSE)
}
outcome <- aggregate_validation_artifacts(
  arguments[[1L]], arguments[[2L]], render = TRUE, strict = TRUE
)
message("Aggregate scientific validation report: ",
        file.path(normalizePath(arguments[[2L]]), "scientific-validation.html"))
quit(status = outcome$status, save = "no")
