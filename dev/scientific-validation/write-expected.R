#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: write-expected.R OUTPUT_FILE", call. = FALSE)
}
matrix <- jsonlite::fromJSON(Sys.getenv("VALIDATION_MATRIX"))$include
scenarios <- c(
  "cross_sectional_null", "cross_sectional_signal", "sparse_zero_inflated",
  "fixed_confounder", "longitudinal_interaction",
  "longitudinal_random_effect", "missing_unbalanced"
)
rows <- lapply(seq_len(nrow(matrix)), function(index) {
  scenario <- matrix$scenario[[index]]
  data.frame(
    engine = matrix$engine[[index]],
    scenario = if (identical(scenario, "all")) scenarios else scenario,
    stringsAsFactors = FALSE
  )
})
utils::write.csv(
  do.call(rbind, rows), arguments[[1L]], row.names = FALSE, na = ""
)
