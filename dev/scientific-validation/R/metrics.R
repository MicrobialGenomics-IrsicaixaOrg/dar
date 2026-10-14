validation_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

validation_result_universe <- function(results, truth, runs) {
  run_keys <- c("engine", "scenario", "replicate")
  scientific_keys <- c(
    "engine", "scenario", "replicate", "taxa_id", "contrast_id"
  )
  truth_keys <- c("scenario", "replicate", "taxa_id", "contrast_id")
  case_keys <- c("scenario", "replicate")
  if (!nrow(runs)) {
    if (nrow(results)) {
      stop("Validation results have no recorded runs.", call. = FALSE)
    }
    return(data.frame())
  }
  if (anyNA(runs[run_keys]) || anyDuplicated(runs[run_keys])) {
    stop("Validation runs contain missing or duplicated keys.", call. = FALSE)
  }
  statuses <- c("success", "engine_error", "expected_skip", "dependency_missing")
  if (anyNA(runs$status) || any(!runs$status %in% statuses)) {
    stop("Validation runs contain unknown execution statuses.", call. = FALSE)
  }
  if (anyNA(truth[truth_keys]) || anyNA(truth$is_null) ||
      anyDuplicated(truth[truth_keys])) {
    stop("Validation truth contains missing or duplicated scientific keys.", call. = FALSE)
  }
  attempted <- runs[
    runs$status %in% c("success", "engine_error"), run_keys, drop = FALSE
  ]
  if (nrow(dplyr::anti_join(attempted, truth[case_keys], by = case_keys))) {
    stop("Validation runs contain cases absent from the truth table.", call. = FALSE)
  }
  if (!nrow(attempted)) {
    if (nrow(results)) {
      stop("Validation results have no attempted compatible runs.", call. = FALSE)
    }
    return(data.frame())
  }
  specifications <- dplyr::bind_rows(lapply(
    unique(attempted$engine), function(engine) {
      config <- validation_engine(engine)
      data.frame(
        engine = engine, effect_metric = config$effect_metric,
        truth_multiplier = config$truth_multiplier, comparable = config$comparable
      )
    }
  ))
  expected <- attempted |>
    dplyr::left_join(specifications, by = "engine") |>
    dplyr::inner_join(truth, by = case_keys, relationship = "many-to-many")
  if (!nrow(results)) {
    expected$effect <- NA_real_
    expected$padj <- NA_real_
    expected$std_error <- NA_real_
    expected$returned <- FALSE
    return(expected)
  }
  if (anyNA(results[scientific_keys]) || anyDuplicated(results[scientific_keys])) {
    stop("Validation results contain missing or duplicated scientific keys.", call. = FALSE)
  }
  if (nrow(dplyr::anti_join(results, expected, by = scientific_keys))) {
    stop("Validation results contain keys absent from the expected universe.", call. = FALSE)
  }
  metadata <- unique(results[c(
    "engine", "effect_metric", "truth_multiplier", "comparable"
  )])
  expected_metadata <- specifications[match(metadata$engine, specifications$engine), ]
  # CSV round trips need not preserve the last floating-point digits.
  if (!identical(metadata$effect_metric, expected_metadata$effect_metric) ||
      !identical(metadata$comparable, expected_metadata$comparable) ||
      !isTRUE(all.equal(
        as.numeric(metadata$truth_multiplier), expected_metadata$truth_multiplier,
        tolerance = 1e-12, check.attributes = FALSE
      ))) {
    stop("Validation results contain inconsistent effect metadata.", call. = FALSE)
  }
  observed <- results[c(scientific_keys, "effect", "padj", "std_error")]
  observed$returned <- TRUE
  expected |>
    dplyr::left_join(observed, by = scientific_keys) |>
    dplyr::mutate(returned = tidyr::replace_na(.data$returned, FALSE))
}

score_validation_results <- function(results, truth, runs, alpha = 0.05) {
  joined <- validation_result_universe(results, truth, runs)
  if (!nrow(joined)) {
    return(data.frame())
  }
  joined |>
    dplyr::group_by(
      .data$engine, .data$scenario, .data$replicate, .data$contrast_id,
      .data$effect_metric, .data$truth_multiplier, .data$comparable
    ) |>
    dplyr::group_modify(function(data, key) {
      discovered <- is.finite(data$padj) & data$padj <= alpha
      signal <- !data$is_null
      null <- data$is_null
      finite <- is.finite(data$effect)
      direction_ok <- finite & signal &
        sign(data$effect) == data$expected_direction
      scaled_truth <- data$true_effect * key$truth_multiplier
      comparable <- isTRUE(key$comparable) && is.finite(key$truth_multiplier)
      estimated_signal <- signal & finite
      errors <- data$effect[estimated_signal] - scaled_truth[estimated_signal]
      bias <- if (comparable && any(estimated_signal)) {
        validation_mean(errors)
      } else {
        NA_real_
      }
      scale <- validation_mean(abs(scaled_truth[estimated_signal]))
      intervals <- estimated_signal & is.finite(data$std_error) & data$std_error >= 0
      coverage <- if (comparable && any(intervals)) {
        keep <- intervals
        validation_mean(
          scaled_truth[keep] >= data$effect[keep] - 1.96 * data$std_error[keep] &
            scaled_truth[keep] <= data$effect[keep] + 1.96 * data$std_error[keep]
        )
      } else {
        NA_real_
      }
      data.frame(
        n_expected = nrow(data),
        n_returned = sum(data$returned),
        n_missing = sum(!data$returned),
        n_signals = sum(signal),
        n_tested = sum(is.finite(data$padj)),
        n_finite = sum(finite),
        n_signal_finite = sum(estimated_signal),
        n_intervals = if (comparable) sum(intervals) else 0L,
        type1_error = if (any(null)) validation_mean(discovered[null]) else NA_real_,
        fdr = if (any(discovered)) sum(discovered & null) / sum(discovered) else 0,
        power = if (any(signal)) validation_mean(discovered[signal]) else NA_real_,
        direction_recovery = if (any(signal)) {
          validation_mean(direction_ok[signal])
        } else {
          NA_real_
        },
        finite_rate = validation_mean(finite),
        bias = bias,
        relative_mae = if (comparable && is.finite(scale) && scale > 0) {
          validation_mean(abs(errors)) / scale
        } else {
          NA_real_
        },
        coverage = coverage
      )
    }) |>
    dplyr::ungroup() |>
    as.data.frame()
}

score_validation_runs <- function(runs) {
  if (!nrow(runs)) {
    return(data.frame())
  }
  attempted <- runs[runs$status %in% c("success", "engine_error"), , drop = FALSE]
  if (!nrow(attempted)) {
    return(data.frame())
  }
  data.frame(
    engine = attempted$engine,
    scenario = attempted$scenario,
    replicate = attempted$replicate,
    contrast_id = "__run__",
    effect_metric = "execution",
    truth_multiplier = NA_real_,
    comparable = FALSE,
    type1_error = NA_real_, fdr = NA_real_, power = NA_real_,
    direction_recovery = NA_real_, finite_rate = NA_real_,
    bias = NA_real_, relative_mae = NA_real_, coverage = NA_real_,
    execution_rate = as.numeric(attempted$status == "success"),
    stringsAsFactors = FALSE
  )
}

validation_interval <- function(values, confidence = 0.95, bounded = TRUE) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) {
    return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
  }
  estimate <- mean(values)
  if (n == 1L || stats::sd(values) == 0) {
    return(c(mean = estimate, lower = estimate, upper = estimate, n = n))
  }
  margin <- stats::qt((1 + confidence) / 2, df = n - 1L) *
    stats::sd(values) / sqrt(n)
  interval <- c(
    mean = estimate,
    lower = estimate - margin,
    upper = estimate + margin,
    n = n
  )
  if (bounded) {
    interval[["lower"]] <- max(0, interval[["lower"]])
    interval[["upper"]] <- min(1, interval[["upper"]])
  }
  interval
}

summarize_validation_metrics <- function(metrics) {
  if (!nrow(metrics)) {
    return(data.frame())
  }
  metric_names <- intersect(
    c(
      "type1_error", "fdr", "power", "direction_recovery", "finite_rate",
      "relative_mae", "coverage", "execution_rate"
    ),
    names(metrics)
  )
  metrics |>
    tidyr::pivot_longer(
      dplyr::all_of(metric_names), names_to = "metric", values_to = "value"
    ) |>
    dplyr::filter(is.finite(.data$value)) |>
    dplyr::group_by(
      .data$engine, .data$scenario, .data$contrast_id, .data$metric
    ) |>
    dplyr::group_modify(function(data, key) {
      interval <- validation_interval(
        data$value,
        bounded = !identical(key$metric, "relative_mae")
      )
      if (identical(key$metric, "relative_mae")) {
        interval[["lower"]] <- max(0, interval[["lower"]])
      }
      data.frame(
        estimate = unname(interval[["mean"]]),
        lower = unname(interval[["lower"]]),
        upper = unname(interval[["upper"]]),
        n_replicates = as.integer(interval[["n"]])
      )
    }) |>
    dplyr::ungroup() |>
    as.data.frame()
}

read_validation_thresholds <- function() {
  thresholds <- utils::read.csv(
    file.path(validation_root_path(), "thresholds.csv"),
    stringsAsFactors = FALSE
  )
  if (!"contrast_pattern" %in% names(thresholds)) {
    thresholds$contrast_pattern <- ".*"
  }
  thresholds
}

evaluate_validation_gates <- function(summary, thresholds = read_validation_thresholds()) {
  if (!nrow(summary)) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(summary)), function(index) {
    item <- summary[index, , drop = FALSE]
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
    candidates <- thresholds[
      thresholds$metric == item$metric & scenario_match & engine_match &
        contrast_match,
      , drop = FALSE
    ]
    if (!nrow(candidates)) {
      return(cbind(item, status = "not_applicable", limit = NA_real_,
                   min_replicates = NA_integer_, notes = NA_character_))
    }
    threshold <- candidates[1L, , drop = FALSE]
    status <- if (item$n_replicates < threshold$min_replicates) {
      "informational"
    } else if (threshold$direction == "upper") {
      if (item$lower > threshold$limit) "fail" else "pass"
    } else {
      if (item$upper < threshold$limit) "fail" else "pass"
    }
    cbind(
      item,
      status = status,
      limit = threshold$limit,
      min_replicates = threshold$min_replicates,
      notes = threshold$notes
    )
  })
  do.call(rbind, rows)
}

direction_reversal_self_test <- function() {
  summary <- data.frame(
    engine = "deseq",
    scenario = "cross_sectional_signal",
    contrast_id = "condition[treated-control]",
    metric = "direction_recovery",
    estimate = 0,
    lower = 0,
    upper = 0,
    n_replicates = 20L,
    stringsAsFactors = FALSE
  )
  gate <- evaluate_validation_gates(summary)
  identical(gate$status, "fail")
}

validation_version_manifest <- function(engines, profile, base_seed) {
  installed <- as.data.frame(utils::installed.packages(), stringsAsFactors = FALSE)
  dependencies <- tryCatch(
    tools::package_dependencies(
      "dar", db = utils::installed.packages(), recursive = TRUE
    )[["dar"]],
    error = function(cnd) character()
  )
  engine_packages <- unique(unlist(lapply(
    engines, function(engine) validation_engine(engine)$packages
  ), use.names = FALSE))
  packages <- unique(c("dar", "BiocVersion", engine_packages, dependencies))
  packages <- packages[packages %in% installed$Package]
  roles <- ifelse(
    packages %in% engine_packages, "engine",
    ifelse(packages == "dar", "package", "dependency")
  )
  repositories <- if ("Repository" %in% names(installed)) {
    installed[packages, "Repository"]
  } else {
    rep(NA_character_, length(packages))
  }
  data.frame(
    component = c("R", packages),
    version = c(as.character(getRversion()), installed[packages, "Version"]),
    role = c("runtime", roles),
    source = c("R runtime", repositories),
    library = c(NA_character_, installed[packages, "LibPath"]),
    engine = paste(engines, collapse = ","),
    profile = profile,
    base_seed = as.integer(base_seed),
    platform = R.version$platform,
    git_sha = Sys.getenv("GITHUB_SHA", NA_character_),
    run_id = Sys.getenv("GITHUB_RUN_ID", NA_character_),
    stringsAsFactors = FALSE
  )
}

write_validation_artifacts <- function(output_dir, results, runs, truth,
                                       manifests, metrics, summary, gates,
                                       versions = data.frame()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  objects <- list(
    results = results, runs = runs, truth = truth, manifests = manifests,
    metrics = metrics, summary = summary, gates = gates, versions = versions
  )
  for (name in names(objects)) {
    utils::write.csv(
      objects[[name]], file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE, na = ""
    )
  }
  saveRDS(objects, file.path(output_dir, "validation-results.rds"))
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(output_dir, "session-info.txt")
  )
  invisible(objects)
}
