validation_baseline_path <- function() {
  file.path(validation_root_path(), "baselines", "accepted-full.csv")
}

read_validation_baseline <- function(path = validation_baseline_path()) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(data.frame())
  }
  baseline <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "engine", "scenario", "contrast_id", "metric", "estimate", "lower",
    "upper", "direction", "tolerance", "source_run_id"
  )
  missing <- setdiff(required, names(baseline))
  if (length(missing)) {
    stop(
      "Validation baseline is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  baseline
}

validation_regression_tolerance <- function(metric) {
  values <- c(
    type1_error = 0.03,
    fdr = 0.03,
    power = 0.05,
    direction_recovery = 0.05,
    finite_rate = 0.05,
    execution_rate = 0.05,
    coverage = 0.05,
    abs_relative_bias = 0.10
  )
  unname(values[metric])
}

compare_validation_baseline <- function(summary, baseline) {
  keys <- c("engine", "scenario", "contrast_id", "metric")
  output_columns <- c(
    keys, "estimate", "lower", "upper", "baseline_estimate",
    "baseline_lower", "baseline_upper", "direction", "tolerance",
    "status", "reason", "source_run_id"
  )
  if (!nrow(baseline)) {
    empty <- as.data.frame(stats::setNames(
      replicate(length(output_columns), logical(), simplify = FALSE),
      output_columns
    ))
    return(empty)
  }

  if (nrow(summary)) {
    current <- summary[c(keys, "estimate", "lower", "upper")]
  } else {
    current <- baseline[0, keys, drop = FALSE]
    current$estimate <- numeric()
    current$lower <- numeric()
    current$upper <- numeric()
  }
  expected <- baseline[c(
    keys, "estimate", "lower", "upper", "direction", "tolerance",
    "source_run_id"
  )]
  names(expected)[match(
    c("estimate", "lower", "upper"), names(expected)
  )] <- c("baseline_estimate", "baseline_lower", "baseline_upper")
  compared <- dplyr::left_join(expected, current, by = keys)
  missing <- is.na(compared$estimate)
  upper_regression <- !missing & compared$direction == "upper" &
    compared$lower > compared$baseline_upper + compared$tolerance
  lower_regression <- !missing & compared$direction == "lower" &
    compared$upper < compared$baseline_lower - compared$tolerance
  compared$status <- dplyr::case_when(
    missing ~ "fail",
    upper_regression | lower_regression ~ "fail",
    TRUE ~ "pass"
  )
  compared$reason <- dplyr::case_when(
    missing ~ "expected metric missing from current run",
    upper_regression ~ "upper-is-worse metric regressed beyond tolerance",
    lower_regression ~ "lower-is-worse metric regressed beyond tolerance",
    TRUE ~ "within accepted Monte Carlo tolerance"
  )

  new_rows <- dplyr::anti_join(current, baseline[keys], by = keys)
  if (nrow(new_rows)) {
    new_rows$baseline_estimate <- NA_real_
    new_rows$baseline_lower <- NA_real_
    new_rows$baseline_upper <- NA_real_
    new_rows$direction <- NA_character_
    new_rows$tolerance <- NA_real_
    new_rows$status <- "informational"
    new_rows$reason <- "new metric is not in the accepted baseline"
    new_rows$source_run_id <- NA_character_
    compared <- dplyr::bind_rows(compared, new_rows)
  }
  compared[output_columns]
}

promote_validation_baseline <- function(summary, gates, source_run_id) {
  keys <- c("engine", "scenario", "contrast_id", "metric")
  accepted <- gates[gates$status == "pass", keys, drop = FALSE]
  baseline <- dplyr::inner_join(summary, accepted, by = keys)
  thresholds <- read_validation_thresholds()
  baseline$direction <- vapply(seq_len(nrow(baseline)), function(index) {
    item <- baseline[index, , drop = FALSE]
    scenario_match <- vapply(
      thresholds$scenario_pattern,
      function(pattern) grepl(pattern, item$scenario, perl = TRUE),
      logical(1)
    )
    engine_match <- vapply(
      thresholds$engine_pattern,
      function(pattern) grepl(pattern, item$engine, perl = TRUE),
      logical(1)
    )
    contrast_match <- vapply(
      thresholds$contrast_pattern,
      function(pattern) grepl(pattern, item$contrast_id, perl = TRUE),
      logical(1)
    )
    matched <- thresholds[
      thresholds$metric == item$metric & scenario_match & engine_match &
        contrast_match,
      , drop = FALSE
    ]
    if (nrow(matched)) matched$direction[[1L]] else NA_character_
  }, character(1))
  baseline$tolerance <- validation_regression_tolerance(baseline$metric)
  baseline$source_run_id <- as.character(source_run_id)
  baseline[c(
    keys, "estimate", "lower", "upper", "n_replicates", "direction",
    "tolerance", "source_run_id"
  )]
}
