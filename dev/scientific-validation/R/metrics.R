validation_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

score_validation_results <- function(results, truth, alpha = 0.05) {
  if (!nrow(results)) {
    return(data.frame())
  }
  scientific_keys <- c(
    "engine", "scenario", "replicate", "taxa_id", "contrast_id"
  )
  truth_keys <- c("scenario", "replicate", "taxa_id", "contrast_id")
  if (anyDuplicated(results[scientific_keys])) {
    stop("Validation results contain duplicated scientific keys.", call. = FALSE)
  }
  joined <- dplyr::left_join(results, truth, by = truth_keys)
  if (any(is.na(joined$is_null))) {
    stop("Validation results contain keys absent from the truth table.", call. = FALSE)
  }

  joined |>
    dplyr::group_by(
      .data$engine, .data$scenario, .data$replicate, .data$contrast_id,
      .data$effect_metric, .data$truth_multiplier, .data$comparable
    ) |>
    dplyr::group_modify(function(data, key) {
      discovered <- !is.na(data$padj) & data$padj <= alpha
      signal <- !data$is_null
      null <- data$is_null
      finite <- is.finite(data$effect)
      direction_ok <- finite & signal &
        sign(data$effect) == data$expected_direction
      scaled_truth <- data$true_effect * key$truth_multiplier
      comparable <- isTRUE(key$comparable) && is.finite(key$truth_multiplier)
      bias <- if (comparable && any(signal & finite)) {
        validation_mean(data$effect[signal & finite] - scaled_truth[signal & finite])
      } else {
        NA_real_
      }
      scale <- validation_mean(abs(scaled_truth[signal & finite]))
      coverage <- if (comparable && any(signal & finite & is.finite(data$std_error))) {
        keep <- signal & finite & is.finite(data$std_error)
        validation_mean(
          scaled_truth[keep] >= data$effect[keep] - 1.96 * data$std_error[keep] &
            scaled_truth[keep] <= data$effect[keep] + 1.96 * data$std_error[keep]
        )
      } else {
        NA_real_
      }
      data.frame(
        type1_error = if (any(null)) validation_mean(discovered[null]) else NA_real_,
        fdr = if (any(discovered)) sum(discovered & null) / sum(discovered) else 0,
        power = if (any(signal)) validation_mean(discovered[signal]) else NA_real_,
        direction_recovery = if (any(signal & finite)) {
          validation_mean(direction_ok[signal & finite])
        } else {
          NA_real_
        },
        finite_rate = validation_mean(finite),
        bias = bias,
        abs_relative_bias = if (is.finite(bias) && is.finite(scale) && scale > 0) {
          abs(bias) / scale
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
    bias = NA_real_, abs_relative_bias = NA_real_, coverage = NA_real_,
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
      "abs_relative_bias", "coverage", "execution_rate"
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
        bounded = !identical(key$metric, "abs_relative_bias")
      )
      if (identical(key$metric, "abs_relative_bias")) {
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
