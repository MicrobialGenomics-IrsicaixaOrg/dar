#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop(
    "Usage: archive-validation.R REPORT_DIR ARCHIVE_DIR RUN_ID",
    call. = FALSE
  )
}
report_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
archive_dir <- normalizePath(arguments[[2L]], mustWork = TRUE)
run_id <- arguments[[3L]]
run_dir <- file.path(archive_dir, "runs", run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
files <- c(
  "scientific-validation.html", "summary.csv", "gates.csv", "runs.csv",
  "versions.csv", "execution_status.csv", "baseline_comparison.csv"
)
available <- files[file.exists(file.path(report_dir, files))]
if (!all(file.copy(
  file.path(report_dir, available), file.path(run_dir, available),
  overwrite = TRUE
))) {
  stop("Unable to copy one or more validation reports.", call. = FALSE)
}
execution <- utils::read.csv(file.path(report_dir, "execution_status.csv"))
gates <- utils::read.csv(file.path(report_dir, "gates.csv"))
baseline <- utils::read.csv(file.path(report_dir, "baseline_comparison.csv"))
failed <- any(grepl("^failed_|^missing_", execution$status)) ||
  any(gates$status == "fail") || any(baseline$status == "fail")
record <- data.frame(
  run_id = run_id,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  git_sha = Sys.getenv("GITHUB_SHA", NA_character_),
  profile = Sys.getenv("VALIDATION_PROFILE", NA_character_),
  base_seed = Sys.getenv("VALIDATION_BASE_SEED", NA_character_),
  status = if (failed) "fail" else "pass",
  path = file.path("runs", run_id),
  stringsAsFactors = FALSE
)
history_path <- file.path(archive_dir, "history.csv")
history <- if (file.exists(history_path)) {
  utils::read.csv(history_path, stringsAsFactors = FALSE)
} else {
  record[0, ]
}
history <- rbind(history[history$run_id != run_id, , drop = FALSE], record)
utils::write.csv(history, history_path, row.names = FALSE, na = "")
jsonlite::write_json(
  as.list(record[1, ]), file.path(archive_dir, "latest.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)
