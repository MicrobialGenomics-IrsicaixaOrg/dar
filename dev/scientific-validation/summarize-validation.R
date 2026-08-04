#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: summarize-validation.R REPORT_DIR", call. = FALSE)
}
report_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
read_report <- function(name) {
  path <- file.path(report_dir, paste0(name, ".csv"))
  if (file.exists(path) && file.info(path)$size > 0L) {
    utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
}
execution <- read_report("execution_status")
gates <- read_report("gates")
baseline <- read_report("baseline_comparison")
count_status <- function(data, status) {
  if (nrow(data) && "status" %in% names(data)) sum(data$status == status) else 0L
}
run_url <- paste0(
  Sys.getenv("GITHUB_SERVER_URL"), "/", Sys.getenv("GITHUB_REPOSITORY"),
  "/actions/runs/", Sys.getenv("GITHUB_RUN_ID")
)
history_url <- paste0(
  Sys.getenv("GITHUB_SERVER_URL"), "/", Sys.getenv("GITHUB_REPOSITORY"),
  "/tree/validation-results"
)
lines <- c(
  "## Scientific validation",
  "",
  sprintf("- Profile: `%s`", Sys.getenv("VALIDATION_PROFILE", "unknown")),
  sprintf("- Seed: `%s`", Sys.getenv("VALIDATION_BASE_SEED", "unknown")),
  sprintf("- Execution: %d success, %d incompatible, %d failed/missing",
          count_status(execution, "success"),
          count_status(execution, "skipped_incompatible"),
          sum(grepl("^failed_|^missing_", execution$status))),
  sprintf("- Statistical gates: %d failed", count_status(gates, "fail")),
  sprintf("- Baseline comparison: %d failed, %d informational",
          count_status(baseline, "fail"), count_status(baseline, "informational")),
  sprintf("- [Workflow artifacts](%s)", run_url),
  sprintf("- [Durable validation history](%s)", history_url)
)
summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary_path)) {
  write(lines, summary_path, append = TRUE)
} else {
  writeLines(lines)
}
