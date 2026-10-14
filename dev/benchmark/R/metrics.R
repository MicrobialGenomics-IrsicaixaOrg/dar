benchmark_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

benchmark_scientific_keys <- function(data) {
  paste(data$taxa_id, data$contrast_id, sep = "\r")
}

# Primary selection metrics use the original J-taxon truth. The U-conditioned
# columns remain visible when removing all-zero taxa changes the estimand.
benchmark_score <- function(results, truth, procedure, cell_id, replicate,
                            status = "success") {
  truth_columns <- c(
    "taxa_id", "contrast_id", "beta_latent", "beta_clr_U",
    "beta_median_U", "retained"
  )
  if (!is.data.frame(truth) || !all(truth_columns %in% names(truth)) ||
      !nrow(truth) || anyNA(truth[c("taxa_id", "contrast_id", "retained")]) ||
      !is.logical(truth$retained) || anyDuplicated(benchmark_scientific_keys(truth)) ||
      length(unique(truth$contrast_id)) != 1L ||
      any(!is.finite(truth$beta_latent)) ||
      any(!is.finite(truth$beta_clr_U[truth$retained])) ||
      any(!is.finite(truth$beta_median_U[truth$retained]))) {
    stop("Benchmark truth has an invalid schema or scientific keys.", call. = FALSE)
  }
  for (column in c("beta_clr_original", "beta_median_original")) {
    if (!column %in% names(truth)) {
      source <- sub("_original$", "_U", column)
      if (any(!is.finite(truth[[source]]))) {
        stop("Original truth is required for removed taxa.", call. = FALSE)
      }
      truth[[column]] <- truth[[source]]
    }
    if (any(!is.finite(truth[[column]]))) {
      stop("Original truth must be finite for every input taxon.", call. = FALSE)
    }
  }
  statuses <- c("success", "engine_error", "timeout", "contract_invalid", "unavailable")
  if (length(status) != 1L || !status %in% statuses ||
      length(procedure) != 1L || is.na(procedure) ||
      length(cell_id) != 1L || is.na(cell_id) ||
      length(replicate) != 1L || !is.finite(replicate)) {
    stop("Benchmark score requires one recorded procedure and execution status.", call. = FALSE)
  }
  required <- c("taxa_id", "contrast_id", "effect_size", "significant")
  reason <- ""
  if (is.null(results)) results <- data.frame()
  schema_ok <- is.data.frame(results) &&
    (!nrow(results) || (all(required %in% names(results)) &&
      is.logical(results$significant) && is.numeric(results$effect_size)))
  if (!schema_ok) {
    reason <- "Invalid result schema."
  } else if (nrow(results)) {
    if (anyNA(results[c("taxa_id", "contrast_id")]) ||
        anyDuplicated(benchmark_scientific_keys(results))) {
      reason <- "Missing or duplicate scientific result keys."
    } else if (any(!benchmark_scientific_keys(results) %in%
                  benchmark_scientific_keys(truth[truth$retained, ]))) {
      reason <- "Result keys are outside the common retained universe."
    } else if (any(results$significant %in% TRUE &
                  (!is.finite(results$effect_size) | results$effect_size == 0))) {
      reason <- "A significant result has no finite non-neutral direction."
    }
  }
  if (nzchar(reason)) status <- "contract_invalid"
  observed <- rep(FALSE, nrow(truth))
  significant <- rep(NA, nrow(truth))
  effect <- rep(NA_real_, nrow(truth))
  qvalue <- rep(NA_real_, nrow(truth))
  has_q <- FALSE
  keys_ok <- schema_ok && nrow(results) > 0 &&
    !anyNA(results[c("taxa_id", "contrast_id")]) &&
    !anyDuplicated(benchmark_scientific_keys(results)) &&
    all(benchmark_scientific_keys(results) %in%
          benchmark_scientific_keys(truth[truth$retained, ]))
  if (keys_ok) {
    index <- match(benchmark_scientific_keys(results), benchmark_scientific_keys(truth))
    observed[index] <- TRUE
    significant[index] <- results$significant
    effect[index] <- results$effect_size
    q_column <- intersect(
      c("adj_p_value", "padj", "p_adj", "qvalue", "p_adjust", "adjusted_pvalue"),
      names(results)
    )
    if (length(q_column)) {
      has_q <- TRUE
      qvalue[index] <- as.numeric(results[[q_column[[1L]]]])
    }
  }
  available <- identical(status, "success")
  selected <- significant %in% TRUE
  objectives <- list(
    truth_latent = c("beta_latent", "beta_latent"),
    truth_clr = c("beta_clr_original", "beta_clr_U"),
    truth_median = c("beta_median_original", "beta_median_U")
  )
  rows <- lapply(names(objectives), function(objective) {
    original <- truth[[objectives[[objective]][[1L]]]]
    current <- truth[[objectives[[objective]][[2L]]]]
    signal <- abs(original) > 1e-10
    signal_U <- truth$retained & is.finite(current) & abs(current) > 1e-10
    selected_signal <- selected & signal
    correct <- selected_signal & is.finite(effect) & sign(effect) == sign(original)
    incorrect <- selected_signal & !correct
    correct_U <- selected & signal_U & is.finite(effect) & sign(effect) == sign(current)
    incorrect_U <- selected & signal_U & !correct_U
    discoveries <- sum(selected)
    false <- sum(selected & !signal)
    false_U <- sum(selected & !signal_U)
    power <- if (any(signal) && available) sum(selected_signal) / sum(signal) else NA_real_
    directional <- if (any(signal) && available) sum(correct) / sum(signal) else NA_real_
    data.frame(
      cell_id = cell_id, replicate = replicate, procedure = procedure,
      truth = objective, contrast_id = unique(truth$contrast_id),
      status = status, available = available, reason = reason,
      n_expected = nrow(truth), n_universe = sum(truth$retained),
      n_output_rows = if (is.data.frame(results)) nrow(results) else NA_integer_,
      n_returned = sum(observed), n_missing = sum(truth$retained & !observed),
      n_significance_missing = sum(observed & is.na(significant)),
      n_effect_finite = sum(observed & is.finite(effect)),
      n_effect_nonfinite = sum(observed & !is.finite(effect)),
      n_q_finite = if (has_q) sum(observed & is.finite(qvalue)) else NA_integer_,
      n_signals = sum(signal), n_signals_U = sum(signal_U),
      n_signal_retained = sum(signal & truth$retained),
      n_signal_returned = sum(signal & observed),
      signal_retention = if (any(signal)) sum(signal & truth$retained) / sum(signal) else NA_real_,
      is_null = !any(signal),
      estimand_shift = if (any(truth$retained)) {
        max(abs(original[truth$retained] - current[truth$retained]))
      } else NA_real_,
      R = if (available) discoveries else NA_integer_,
      FP = if (available) false else NA_integer_,
      TP = if (available) sum(selected_signal) else NA_integer_,
      n_direction_wrong = if (available) sum(incorrect) else NA_integer_,
      fdp = if (available) false / max(discoveries, 1L) else NA_real_,
      fdp_U = if (available) false_U / max(discoveries, 1L) else NA_real_,
      power_endtoend = power,
      power_U = if (available && any(signal_U)) sum(selected & signal_U) / sum(signal_U) else NA_real_,
      directional_power = directional,
      directional_power_U = if (available && any(signal_U)) sum(correct_U) / sum(signal_U) else NA_real_,
      directional_fdp = if (available) (false + sum(incorrect)) / max(discoveries, 1L) else NA_real_,
      directional_fdp_U = if (available) (false_U + sum(incorrect_U)) / max(discoveries, 1L) else NA_real_,
      error_directional = if (available && any(selected_signal)) sum(incorrect) / sum(selected_signal) else NA_real_,
      conservative_power = if (!any(signal)) NA_real_ else if (available) power else 0,
      conservative_directional_power = if (!any(signal)) NA_real_ else if (available) directional else 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

benchmark_null_interval <- function(events, confidence = 0.95) {
  events <- events[!is.na(events)]
  if (!length(events)) return(c(lower = NA_real_, upper = NA_real_, upper_one_sided = NA_real_))
  exact <- stats::binom.test(sum(events), length(events), conf.level = confidence)$conf.int
  one_sided <- stats::binom.test(
    sum(events), length(events), alternative = "less", conf.level = confidence
  )$conf.int[[2L]]
  c(lower = exact[[1L]], upper = exact[[2L]], upper_one_sided = one_sided)
}

benchmark_paired_interval <- function(values, bootstrap) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
  boot <- replicate(bootstrap, mean(sample(values, length(values), replace = TRUE)))
  c(mean = mean(values), lower = unname(stats::quantile(boot, 0.025)),
    upper = unname(stats::quantile(boot, 0.975)), n = length(values))
}

# Separate cells are never pooled. Planned-but-unavailable datasets widen the
# identification range rather than becoming zero-FDP successes.
benchmark_summarize <- function(metrics, bootstrap = 2000L,
                                planned_replicates = 20L) {
  if (!is.data.frame(metrics) || !nrow(metrics)) return(data.frame())
  required <- c("cell_id", "replicate", "procedure", "truth", "available", "is_null",
                "fdp", "directional_power", "power_endtoend", "R", "signal_retention",
                "directional_fdp", "error_directional", "estimand_shift")
  if (!all(required %in% names(metrics)) ||
      anyNA(metrics[c("cell_id", "replicate", "procedure", "truth")]) ||
      anyDuplicated(metrics[c("cell_id", "replicate", "procedure", "truth")]) ||
      length(bootstrap) != 1L || !is.finite(bootstrap) || bootstrap < 1 ||
      length(planned_replicates) != 1L || planned_replicates < 1L) {
    stop("Invalid benchmark metric table or summary settings.", call. = FALSE)
  }
  bounded <- c("fdp", "directional_power", "power_endtoend", "signal_retention",
               "directional_fdp", "error_directional")
  if (any(vapply(metrics[bounded], function(x) {
    !is.numeric(x) || any(is.finite(x) & (x < 0 | x > 1))
  }, logical(1)))) {
    stop("Benchmark selection metrics must be probabilities or NA.", call. = FALSE)
  }
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  keys <- unique(metrics[c("cell_id", "procedure", "truth")])
  keys <- keys[order(keys$cell_id, keys$procedure, keys$truth), ]
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    data <- metrics[metrics$cell_id == key$cell_id & metrics$procedure == key$procedure &
                      metrics$truth == key$truth, ]
    data <- data[order(data$replicate), ]
    if (nrow(data) > planned_replicates || !is.logical(data$available) || anyNA(data$available) ||
        anyNA(data$is_null) || length(unique(data$is_null)) != 1L ||
        any(data$available & (!is.finite(data$fdp) | !is.finite(data$R))) ||
        any(data$available & !data$is_null &
              (!is.finite(data$directional_power) | !is.finite(data$power_endtoend)))) {
      stop("Inconsistent or excessive dataset records in summary.", call. = FALSE)
    }
    valid <- data$available
    null <- isTRUE(data$is_null[[1L]])
    n_valid <- sum(valid)
    null_ci <- if (null) benchmark_null_interval(data$R[valid] > 0) else rep(NA_real_, 3L)
    cell_index <- suppressWarnings(as.integer(sub("^P", "", key$cell_id)))
    if (is.na(cell_index)) cell_index <- match(key$cell_id, sort(unique(metrics$cell_id)))
    set.seed(20271014L + cell_index)
    ref <- metrics[metrics$cell_id == key$cell_id & metrics$procedure == "linda_pilot" &
                     metrics$truth == key$truth, ]
    paired <- merge(data, ref, by = "replicate", suffixes = c("", "_reference"))
    paired <- paired[paired$available & paired$available_reference, ]
    fdp_pair <- benchmark_paired_interval(paired$fdp - paired$fdp_reference, bootstrap)
    power_pair <- benchmark_paired_interval(
      paired$directional_power - paired$directional_power_reference, bootstrap
    )
    data.frame(
      cell_id = key$cell_id, procedure = key$procedure, truth = key$truth,
      n_planned = planned_replicates, n_recorded = nrow(data), n_valid = n_valid,
      n_unavailable = planned_replicates - n_valid, is_null = null,
      fdr = benchmark_mean(data$fdp[valid]),
      fdr_identification_lower = sum(data$fdp[valid]) / planned_replicates,
      fdr_identification_upper = (sum(data$fdp[valid]) + planned_replicates - n_valid) / planned_replicates,
      null_events = if (null) sum(data$R[valid] > 0) else NA_integer_,
      null_lower = null_ci[[1L]], null_upper = null_ci[[2L]],
      null_upper_one_sided = null_ci[[3L]],
      power_endtoend = benchmark_mean(data$power_endtoend[valid]),
      directional_power = benchmark_mean(data$directional_power[valid]),
      conservative_power = if (null) NA_real_ else sum(data$power_endtoend[valid]) / planned_replicates,
      conservative_directional_power = if (null) NA_real_ else sum(data$directional_power[valid]) / planned_replicates,
      directional_fdp = benchmark_mean(data$directional_fdp[valid]),
      error_directional = benchmark_mean(data$error_directional[valid]),
      mean_R = benchmark_mean(data$R[valid]),
      signal_retention = benchmark_mean(data$signal_retention),
      max_estimand_shift = if (any(is.finite(data$estimand_shift))) max(data$estimand_shift, na.rm = TRUE) else NA_real_,
      n_paired_fdp = unname(fdp_pair[["n"]]),
      delta_fdp = unname(fdp_pair[["mean"]]),
      delta_fdp_lower = unname(fdp_pair[["lower"]]),
      delta_fdp_upper = unname(fdp_pair[["upper"]]),
      n_paired_power = unname(power_pair[["n"]]),
      delta_directional_power = unname(power_pair[["mean"]]),
      delta_directional_power_lower = unname(power_pair[["lower"]]),
      delta_directional_power_upper = unname(power_pair[["upper"]]),
      paired_reference_directional_power = benchmark_mean(paired$directional_power_reference),
      paired_procedure_directional_power = benchmark_mean(paired$directional_power),
      paired_procedure_fdp = benchmark_mean(paired$fdp),
      interval_label = "exploratory", stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
