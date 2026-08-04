#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validation_root <- normalizePath(
  dirname(sub("^--file=", "", script_argument[[1L]])), mustWork = TRUE
)
options(dar.validation.root = validation_root)
source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "metrics.R"), local = FALSE)
source(file.path(validation_root, "R", "baseline.R"), local = FALSE)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: promote-baseline.R REPORT_DIR SOURCE_RUN_ID", call. = FALSE)
}
report_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
summary <- utils::read.csv(file.path(report_dir, "summary.csv"))
gates <- utils::read.csv(file.path(report_dir, "gates.csv"))
baseline <- promote_validation_baseline(summary, gates, arguments[[2L]])
output_dir <- file.path(validation_root, "baselines")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  baseline, file.path(output_dir, "accepted-full.csv"),
  row.names = FALSE, na = ""
)
metadata <- list(
  schema_version = 1L,
  source_run_id = arguments[[2L]],
  promoted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  thresholds = "../thresholds.csv"
)
jsonlite::write_json(
  metadata, file.path(output_dir, "accepted-full.json"),
  auto_unbox = TRUE, pretty = TRUE
)
message("Accepted baseline written to ", output_dir)
