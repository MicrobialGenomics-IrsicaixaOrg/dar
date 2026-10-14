# Native fits are isolated; their decisions and effects are never rewritten.
benchmark_engines <- function() {
  c("ancom_pilot", "linda_pilot", "aldex_pilot", "maaslin_pilot", "wilcox_pilot")
}

benchmark_panel <- function() benchmark_engines()[seq_len(4L)]

benchmark_add_engine <- function(rec, engine) {
  switch(engine,
    ancom_pilot = dar::step_ancom(
      rec, p_adj_method = "BH", prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
      struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, n_cl = 1,
      verbose = FALSE, rarefy = FALSE, id = engine,
      engine_args = list(fit = list(pseudo_sens = TRUE))
    ),
    linda_pilot = dar::step_linda(
      rec, prev_filter = 0, mean_abund_filter = 0, max_abund_filter = 0,
      winsorize = TRUE, outlier_pct = 0.03, adaptive = TRUE,
      zero_handling = "pseudo-count", pseudo_count = 0.5, corr_cut = 0.1,
      p_adj_method = "BH", alpha = 0.05, n_cpus = 1, rarefy = FALSE, id = engine
    ),
    aldex_pilot = dar::step_aldex(
      rec, max_significance = 0.05, mc.samples = 128, denom = "all",
      rarefy = FALSE, id = engine
    ),
    maaslin_pilot = dar::step_maaslin(
      rec, min_abundance = 0, min_prevalence = 0, min_variance = 0,
      normalization = "TSS", transform = "LOG", max_significance = 0.05,
      correction = "BH", standardize = TRUE, median_comparison_abundance = TRUE,
      rarefy = FALSE, id = engine,
      engine_args = list(fit = list(cores = 1), contrast = list(subtract_median = TRUE))
    ),
    wilcox_pilot = dar::step_wilcox(
      rec, norm_method = "clr", max_significance = 0.05, p_adj_method = "BH",
      rarefy = FALSE, id = engine, engine_args = list(test = list(exact = FALSE))
    ),
    stop("Unknown benchmark engine.", call. = FALSE)
  )
}

benchmark_engine_recipe <- function(rec, engine) {
  if (length(rec@steps)) stop("Input recipe must not contain DA steps.", call. = FALSE)
  benchmark_add_engine(rec, engine)
}

benchmark_empty_results <- function() {
  data.frame(taxa_id = character(), contrast_id = character(), effect_size = numeric(),
             significant = logical(), adj_p_value = numeric(), stringsAsFactors = FALSE)
}

benchmark_validate_calls <- function(results, simulation) {
  required <- c("taxa_id", "contrast_id", "effect_size", "significant")
  if (!is.data.frame(results) || !all(required %in% names(results)) ||
      !is.logical(results$significant) || !is.numeric(results$effect_size) ||
      anyNA(results[c("taxa_id", "contrast_id")]) ||
      anyDuplicated(results[c("taxa_id", "contrast_id")]) ||
      any(!results$taxa_id %in% simulation$truth$taxa_id[simulation$truth$retained]) ||
      any(!results$contrast_id %in% simulation$contrast_plan$contrast_id)) {
    stop("Native result keys do not match the common scientific universe.", call. = FALSE)
  }
  selected <- results$significant %in% TRUE
  if (any(selected & (!is.finite(results$effect_size) | results$effect_size == 0))) {
    stop("Significant native calls lack a finite non-neutral direction.", call. = FALSE)
  }
  invisible(TRUE)
}

# Compare the dispatched LinDA body with its dependency without rewriting either.
# Only the verified capitalization repair may differ; NULL and missing arguments
# are compared in place, since deleting either changes the native algorithm.
benchmark_linda_correction_count <- function(native, effective) {
  difference <- function(original, dispatched) {
    if (identical(original, dispatched)) return(0L)
    if (identical(original, quote(zero.handling <- "Imputation")) &&
        identical(dispatched, quote(zero.handling <- "imputation"))) return(1L)
    if (!is.call(original) || !is.call(dispatched) ||
        length(original) != length(dispatched) ||
        !identical(attributes(original), attributes(dispatched))) return(Inf)
    changes <- vapply(seq_along(original), function(index) {
      original_missing <- rlang::is_missing(original[[index]])
      dispatched_missing <- rlang::is_missing(dispatched[[index]])
      if (original_missing || dispatched_missing) {
        return(if (original_missing && dispatched_missing) 0 else Inf)
      }
      difference(original[[index]], dispatched[[index]])
    }, numeric(1))
    total <- sum(changes)
    if (is.finite(total)) as.integer(total) else Inf
  }
  difference(body(native), body(effective))
}

benchmark_linda_adaptive_signature <- function(fun) {
  counts <- c(bug = 0L, corrected = 0L, comparison = 0L)
  inspect <- function(node) {
    if (identical(node, quote(zero.handling <- "Imputation"))) {
      counts[["bug"]] <<- counts[["bug"]] + 1L
    } else if (identical(node, quote(zero.handling <- "imputation"))) {
      counts[["corrected"]] <<- counts[["corrected"]] + 1L
    } else if (identical(node, quote(zero.handling == "imputation"))) {
      counts[["comparison"]] <<- counts[["comparison"]] + 1L
    }
    if (is.call(node)) {
      for (index in seq_along(node)) {
        if (!rlang::is_missing(node[[index]])) inspect(node[[index]])
      }
    }
    invisible(NULL)
  }
  inspect(body(fun))
  list(bug_present = counts[["bug"]] > 0L,
       corrected = counts[["bug"]] == 0L && counts[["corrected"]] == 1L &&
         counts[["comparison"]] == 1L)
}

benchmark_engine_implementation <- function(
    engine, stage, effective = dar:::engine_stage_function(engine, stage)) {
  spec <- dar:::engine_stage_spec(engine, stage)
  native <- getExportedValue(spec$package, spec$fun)
  signature <- function(fun) list(formals = formals(fun), body = body(fun))
  changes <- if (identical(engine, "linda") && identical(stage, "fit")) {
    benchmark_linda_correction_count(native, effective)
  } else 0L
  adaptive <- if (identical(engine, "linda") && identical(stage, "fit")) {
    benchmark_linda_adaptive_signature(effective)
  } else list(bug_present = NA, corrected = NA)
  unchanged_formals <- identical(formals(native), formals(effective))
  unchanged_namespace <- identical(environment(native), environment(effective))
  version <- as.character(utils::packageVersion(spec$package))
  implementation <- if (identical(signature(native), signature(effective))) {
    "native"
  } else if (identical(engine, "linda") && identical(stage, "fit") &&
             identical(version, "1.4") && changes == 1L &&
             unchanged_formals && unchanged_namespace) {
    "dar_local_adaptive_zero_fix"
  } else "unverified"
  list(
    package = spec$package, function_name = spec$fun, package_version = version,
    effective_sha256 = digest::digest(signature(effective), algo = "sha256"),
    native_sha256 = digest::digest(signature(native), algo = "sha256"),
    formals_unchanged = unchanged_formals, namespace_unchanged = unchanged_namespace,
    implementation = implementation, linda_correction_count = changes,
    linda_adaptive_bug_present = adaptive$bug_present,
    linda_adaptive_corrected = adaptive$corrected
  )
}

# Trace observes the registered dispatcher in the worker's namespace only.
# It does not replace the native function or modify any argument/result.
benchmark_trace_enter <- function(frame) {
  capture <- get(".benchmark_capture", envir = .GlobalEnv)
  index <- length(capture$stages) + 1L
  frame$.benchmark_stage_index <- index
  extra <- frame$engine_args[[frame$stage]]
  if (is.null(extra)) extra <- list()
  defaults <- frame$defaults[setdiff(names(frame$defaults), names(extra))]
  spec <- dar:::engine_stage_spec(frame$engine, frame$stage)
  effective <- dar:::engine_stage_function(frame$engine, frame$stage)
  capture$stages[[index]] <- list(
    engine = frame$engine, stage = frame$stage, package = spec$package,
    function_name = spec$fun, arguments = c(frame$fixed, defaults, extra),
    native_formals = formals(effective),
    effective_function = benchmark_engine_implementation(frame$engine, frame$stage, effective),
    completed = FALSE
  )
}

benchmark_trace_exit <- function(frame, value) {
  capture <- get(".benchmark_capture", envir = .GlobalEnv)
  index <- frame$.benchmark_stage_index
  capture$stages[[index]]$result <- value
  capture$stages[[index]]$completed <- TRUE
}

benchmark_maaslin_contrast_exit <- function(value) {
  capture <- get(".benchmark_capture", envir = .GlobalEnv)
  capture$maaslin_contrasts[[length(capture$maaslin_contrasts) + 1L]] <- value
}

benchmark_stage <- function(stages, stage) {
  matching <- Filter(function(x) identical(x$stage, stage), stages)
  if (length(matching) != 1L || !isTRUE(matching[[1L]]$completed)) {
    stop("Expected one completed native stage: ", stage, call. = FALSE)
  }
  matching[[1L]]
}

# Verify the executed zero transformation from the returned processed counts.
# Reuse the dependency's winsorization, but never refit its depth or DA models.
benchmark_linda_zero_handling <- function(fit) {
  arguments <- fit$arguments
  original <- arguments$feature.dat
  observed <- fit$result$feature.dat.use
  invalid <- function(detail) list(passed = FALSE, mode = "unverified", detail = detail)
  if (!is.matrix(original) || !is.matrix(observed) || !length(observed) ||
      is.null(rownames(observed)) || is.null(colnames(observed)) ||
      anyDuplicated(rownames(observed)) || anyDuplicated(colnames(observed)) ||
      any(!rownames(observed) %in% rownames(original)) ||
      any(!colnames(observed) %in% colnames(original)) ||
      !identical(arguments$feature.dat.type, "count")) {
    return(invalid("LinDA returned no identifiable processed count matrix."))
  }
  before <- original[rownames(observed), colnames(observed), drop = FALSE]
  if (isTRUE(arguments$is.winsor)) {
    before <- getFromNamespace("winsor.fun", "MicrobiomeStat")(
      before, 1 - arguments$outlier.pct, "count"
    )
  }
  matches <- function(expected) isTRUE(all.equal(
    observed, expected, tolerance = 1e-10, check.attributes = FALSE
  ))
  n_zeros <- sum(before == 0)
  if (!n_zeros) {
    return(list(passed = matches(before), mode = "no_zeros",
                detail = "No zeros after filtering/winsorization; zero handling was inactive."))
  }
  depth <- colSums(before)
  if (any(!is.finite(depth)) || any(depth <= 0) ||
      !is.numeric(arguments$pseudo.cnt) || length(arguments$pseudo.cnt) != 1L ||
      !is.finite(arguments$pseudo.cnt) || arguments$pseudo.cnt <= 0) {
    return(invalid("LinDA zero-handling candidates cannot be reconstructed safely."))
  }
  zero_depth <- matrix(rep(depth, nrow(before)), nrow = nrow(before), byrow = TRUE)
  zero_depth[before > 0] <- 0
  # Tied maxima have the same depth, so a deterministic tie avoids audit RNG use.
  denominator <- depth[max.col(zero_depth, ties.method = "first")]
  imputed <- before + zero_depth / denominator
  pseudo_counted <- before + arguments$pseudo.cnt
  modes <- c(imputation = matches(imputed), pseudo_count = matches(pseudo_counted))
  if (sum(modes) != 1L) {
    return(invalid("Returned LinDA counts do not identify one native zero transformation."))
  }
  mode <- names(modes)[modes]
  list(passed = TRUE, mode = mode,
       detail = sprintf("Returned processed counts match native %s (%d pre-treatment zeros).",
                        mode, n_zeros))
}

benchmark_native_checks <- function(engine, stages, results, simulation,
                                    maaslin_contrasts = list()) {
  checks <- list()
  add <- function(name, ok, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, passed = isTRUE(ok), detail = detail, stringsAsFactors = FALSE
    )
  }
  if (engine == "ancom_pilot") {
    fit <- benchmark_stage(stages, "fit")
    add("ancom_robust_calls", identical(fit$arguments$pseudo_sens, TRUE),
        "Sensitivity filtering is explicitly enabled; dar uses native diff_robust.")
    raw <- fit$result$res
    basic <- "diff_conditiontreated"
    sensitivity <- "passed_ss_conditiontreated"
    robust <- "diff_robust_conditiontreated"
    if (!all(c("taxon", basic, sensitivity, robust) %in% names(raw))) {
      stop("Cannot identify native ANCOM robust decision and sensitivity flag.", call. = FALSE)
    }
    add("ancom_robust_conjunction", identical(as.logical(raw[[robust]]),
          as.logical(raw[[basic]]) & as.logical(raw[[sensitivity]])),
        "Native diff_robust equals native diff AND passed_ss; q-values are not re-thresholded.")
    idx <- match(results$taxa_id, raw$taxon)
    add("ancom_native_significance", !anyNA(idx) &&
          identical(as.logical(results$significant), as.logical(raw[[robust]][idx])),
        "Tidy significance equals the native coefficient diff_robust decision.")
  } else if (engine == "linda_pilot") {
    fit <- benchmark_stage(stages, "fit")
    expected <- benchmark_engine_implementation("linda", "fit")
    add("linda_corrected_engine", identical(fit$effective_function, expected) &&
          identical(expected$implementation, "dar_local_adaptive_zero_fix") &&
          expected$linda_correction_count == 1L &&
          !isTRUE(expected$linda_adaptive_bug_present) &&
          isTRUE(expected$linda_adaptive_corrected) &&
          isTRUE(expected$formals_unchanged) && isTRUE(expected$namespace_unchanged),
        "The dispatched MicrobiomeStat 1.4 function differs only by dar's local adaptive-zero repair.")
    add("linda_native_controls", isTRUE(fit$arguments$is.winsor) &&
          isTRUE(fit$arguments$adaptive), "Native winsorization and adaptive zeros retained.")
    zeros <- benchmark_linda_zero_handling(fit)
    add("linda_executed_zero_handling", zeros$passed, zeros$detail)
  } else if (engine == "aldex_pilot") {
    clr <- benchmark_stage(stages, "clr")
    fit <- benchmark_stage(stages, "fit")
    add("aldex_monte_carlo", identical(clr$arguments$mc.samples, 128) &&
          identical(clr$arguments$denom, "all"), "128 draws, all-feature CLR denominator.")
    add("aldex_BH", identical(fit$arguments$fdr.method, "BH"), "Native GLM uses BH.")
  } else if (engine == "maaslin_pilot") {
    fit <- benchmark_stage(stages, "fit")
    contrast <- benchmark_stage(stages, "contrast")
    add("maaslin_abundance_only", identical(contrast$arguments$evaluate_only, "abundance"),
        "Only abundance enters the final contrast and voting.")
    add("maaslin_subtract_median", isTRUE(contrast$arguments$subtract_median),
        "The coefficient is requested relative to its tested median.")
    tested <- contrast$result$fit_data_abundance$results
    native <- dplyr::bind_rows(Filter(is.data.frame, maaslin_contrasts))
    if ("model" %in% names(native)) native <- native[native$model == "abundance", , drop = FALSE]
    if (!nrow(native) || !all(c("feature", "coef", "rhs") %in% names(native)) ||
        !all(c("feature", "coef") %in% names(tested))) {
      stop("Cannot verify MaAsLin native contrast centering.", call. = FALSE)
    }
    idx <- match(tested$feature, native$feature)
    keep <- !is.na(idx) & is.finite(tested$coef) & is.finite(native$coef[idx]) &
      is.finite(native$rhs[idx])
    add("maaslin_centered_effect", any(keep) && isTRUE(all.equal(
      as.numeric(tested$coef[keep]), as.numeric(native$coef[idx[keep]] - native$rhs[idx[keep]]),
      tolerance = 1e-8, check.attributes = FALSE
    )), "Contrast coefficient equals the inner native coefficient minus its own returned RHS.")
  } else {
    tests <- Filter(function(x) identical(x$stage, "test"), stages)
    add("wilcox_approximate", length(tests) > 0L && all(vapply(tests, function(x) {
      isTRUE(x$completed) && identical(x$arguments$exact, FALSE)
    }, logical(1))),
        "Rank-test approximation is explicit.")
  }
  plan <- simulation$contrast_plan
  add("contrast_orientation", nrow(plan) == 1L && plan$numerator == "treated" &&
        plan$denominator == "control", "The planned numerator is treated and denominator control.")
  if (simulation$specification$mode == "preflight") {
    truth <- simulation$truth
    positive <- results$effect_size[results$taxa_id %in% truth$taxa_id[truth$beta_latent > 1e-10]]
    negative <- results$effect_size[results$taxa_id %in% truth$taxa_id[truth$beta_latent < -1e-10]]
    positive <- positive[is.finite(positive)]; negative <- negative[is.finite(negative)]
    add("strong_signal_orientation", length(positive) > 0L && length(negative) > 0L &&
          stats::median(positive) > 0 && stats::median(negative) < 0,
        "Median effects follow both implanted strong-signal directions; this is a smoke test, not performance evidence.")
  }
  do.call(rbind, checks)
}

benchmark_fit_worker <- function(job) {
  suppressPackageStartupMessages(library(dar))
  simulation <- readRDS(job$simulation_file)
  rec <- benchmark_engine_recipe(simulation$recipe, job$engine)
  capture <- new.env(parent = emptyenv())
  capture$stages <- list()
  capture$maaslin_contrasts <- list()
  assign(".benchmark_capture", capture, envir = .GlobalEnv)
  suppressMessages(trace(
    "exec_engine_stage", where = asNamespace("dar"), print = FALSE,
    tracer = quote(get("benchmark_trace_enter", envir = .GlobalEnv)(environment())),
    exit = quote(get("benchmark_trace_exit", envir = .GlobalEnv)(environment(), returnValue()))
  ))
  on.exit(suppressMessages(untrace("exec_engine_stage", where = asNamespace("dar"))), add = TRUE)
  if (job$engine == "maaslin_pilot") {
    suppressMessages(trace(
      "maaslin_contrast_test_runner", where = asNamespace("maaslin3"), print = FALSE,
      exit = quote(get("benchmark_maaslin_contrast_exit", envir = .GlobalEnv)(returnValue()))
    ))
    on.exit(suppressMessages(untrace("maaslin_contrast_test_runner",
                                    where = asNamespace("maaslin3"))), add = TRUE)
  }
  warning_messages <- character()
  set.seed(job$engine_seed)
  rng_start <- .Random.seed
  start <- proc.time()
  error <- NULL
  checks <- data.frame()
  prepared <- NULL
  results <- benchmark_empty_results()
  stage <- "prep"
  outcome <- tryCatch(withCallingHandlers({
    prepared <- dar::prep(rec, parallel = FALSE)
    stage <- "contract"
    results <- as.data.frame(dar::tidy_results(prepared), stringsAsFactors = FALSE)
    benchmark_validate_calls(results, simulation)
    checks <- benchmark_native_checks(job$engine, capture$stages, results, simulation,
                                      capture$maaslin_contrasts)
    if (any(!checks$passed)) stop("Native preflight invariants failed.", call. = FALSE)
    list(prepared = prepared, results = results, checks = checks)
  }, warning = function(cnd) {
    warning_messages <<- c(warning_messages, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  }), error = function(cnd) {
    error <<- list(message = conditionMessage(cnd), class = class(cnd),
                  calls = vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), character(1)))
    NULL
  })
  timing <- proc.time() - start
  # Vector-heap high-water mark is distinct from externally sampled process RSS.
  memory <- gc()["Vcells", 6L]
  outcome <- list(
    status = if (is.null(error)) "success" else if (stage == "contract") "contract_invalid" else "engine_error",
    reason = if (is.null(error)) "" else error$message,
    error = error, warnings = unique(warning_messages),
    prepared = prepared,
    results = results,
    checks = checks,
    stages = capture$stages, maaslin_contrasts = capture$maaslin_contrasts,
    rng_start = rng_start, rng_end = .Random.seed,
    elapsed_seconds = unname(timing[["elapsed"]]),
    user_seconds = unname(timing[["user.self"]]),
    system_seconds = unname(timing[["sys.self"]]),
    peak_vector_heap_mb = as.numeric(memory), peak_rss_mb = NA_real_
  )
  saveRDS(outcome, job$outcome_file)
  list(status = outcome$status, reason = outcome$reason)
}

benchmark_run_engine <- function(job, code_root, timeout = 600, ram_limit_mb = Inf) {
  start <- proc.time()[["elapsed"]]
  peak_rss <- NA_real_
  timed_out <- FALSE
  resource_stop <- FALSE
  failure <- tryCatch({
    process <- callr::r_bg(function(job, code_root) {
      source(file.path(code_root, "R", "engines.R"), local = .GlobalEnv)
      benchmark_fit_worker(job)
    }, args = list(job = job, code_root = code_root), libpath = .libPaths(),
    user_profile = FALSE, system_profile = FALSE,
    env = c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1", BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "FALSE"),
    stdout = job$stdout_file, stderr = job$stderr_file)
    on.exit(if (process$is_alive()) process$kill_tree(), add = TRUE)
    while (process$is_alive()) {
      rss <- tryCatch(ps::ps_memory_info(ps::ps_handle(process$get_pid()))[["rss"]] / 1024^2,
                      error = function(cnd) NA_real_)
      if (is.finite(rss)) peak_rss <- max(c(peak_rss, rss), na.rm = TRUE)
      parent_rss <- tryCatch(ps::ps_memory_info()[["rss"]] / 1024^2,
                             error = function(cnd) NA_real_)
      if (is.finite(rss) && is.finite(parent_rss) && rss + parent_rss > ram_limit_mb) {
        resource_stop <- TRUE
        process$kill_tree()
        stop("Benchmark parent and worker exceeded the fixed 75% physical-RAM budget.", call. = FALSE)
      }
      if (proc.time()[["elapsed"]] - start > timeout) {
        timed_out <- TRUE
        process$kill_tree()
        stop("Native fit exceeded its fixed timeout.", call. = FALSE)
      }
      Sys.sleep(0.2)
    }
    process$get_result()
    NULL
  }, error = function(cnd) cnd)
  if (is.null(failure)) {
    outcome <- readRDS(job$outcome_file)
    outcome$peak_rss_mb <- peak_rss
    outcome$worker_elapsed_seconds <- proc.time()[["elapsed"]] - start
    outcome$resource_stop <- FALSE
    saveRDS(outcome, job$outcome_file)
    return(outcome)
  }
  outcome <- list(
    status = if (timed_out) "timeout" else "engine_error",
    reason = conditionMessage(failure), error = list(class = class(failure)),
    warnings = character(), prepared = NULL, results = benchmark_empty_results(),
    checks = data.frame(), stages = list(),
    elapsed_seconds = proc.time()[["elapsed"]] - start,
    user_seconds = NA_real_, system_seconds = NA_real_,
    peak_vector_heap_mb = NA_real_, peak_rss_mb = peak_rss,
    worker_elapsed_seconds = proc.time()[["elapsed"]] - start,
    resource_stop = resource_stop
  )
  saveRDS(outcome, job$outcome_file)
  outcome
}

benchmark_model_identity <- function(rec) {
  model <- dar::get_model(rec)
  model$formula <- paste(deparse(model$formula), collapse = " ")
  model
}

benchmark_combine_prepared <- function(outcomes, simulation, include_wilcox = TRUE) {
  ids <- if (include_wilcox) benchmark_engines() else benchmark_panel()
  if (!all(benchmark_panel() %in% names(outcomes)) ||
      any(vapply(outcomes[benchmark_panel()], function(x) x$status != "success", logical(1)))) {
    stop("The fixed four-engine panel is unavailable.", call. = FALSE)
  }
  ids <- ids[ids %in% names(outcomes) & vapply(outcomes[ids], function(x) {
    !is.null(x) && identical(x$status, "success")
  }, logical(1))]
  reference <- outcomes[[benchmark_panel()[[1L]]]]$prepared
  ref_phy <- dar:::analysis_phyloseq(reference)
  ref_counts <- dar:::recipe_count_matrix(reference, phy = ref_phy)
  rec <- methods::as(reference, "Recipe")
  rec@steps <- list()
  results <- list()
  for (id in ids) {
    prepared <- outcomes[[id]]$prepared
    phy <- dar:::analysis_phyloseq(prepared)
    if (!identical(benchmark_model_identity(prepared), benchmark_model_identity(reference)) ||
        !identical(prepared@execution$contrasts, reference@execution$contrasts) ||
        !identical(prepared@execution$dropped_samples, reference@execution$dropped_samples) ||
        !identical(dar:::recipe_count_matrix(prepared, phy = phy), ref_counts) ||
        !identical(as(phyloseq::sample_data(phy), "data.frame"),
                   as(phyloseq::sample_data(ref_phy), "data.frame")) ||
        !identical(prepared@execution$executed_steps, id) ||
        nrow(prepared@execution$skipped_steps) > 0L) {
      stop("Isolated fits do not share identical model/cohort/contrast metadata.", call. = FALSE)
    }
    rec <- benchmark_add_engine(rec, id)
    results[[id]] <- prepared@results[[id]]
  }
  execution <- reference@execution
  execution$executed_steps <- ids
  execution$skipped_steps <- reference@execution$skipped_steps
  combined <- dar:::prep_recipe(rec, results, list(), execution = execution)
  methods::validObject(combined)
  tidy <- as.data.frame(dar::tidy_results(combined))
  for (id in ids) benchmark_validate_calls(tidy[tidy$step_id == id, , drop = FALSE], simulation)
  combined
}

benchmark_explicit_consensus <- function(prepared, cutoff) {
  data <- as.data.frame(dar::tidy_results(prepared, steps = benchmark_panel(),
                                        significant_only = TRUE))
  if (!nrow(data)) return(data.frame(taxa_id = character(), contrast_id = character(),
                                   effect = character()))
  data$effect <- ifelse(data$effect_size > 0, "up", "down")
  keys <- data[c("taxa_id", "contrast_id", "effect")]
  score <- aggregate(rep(1L, nrow(keys)), keys, sum)
  score[score$x >= cutoff, c("taxa_id", "contrast_id", "effect"), drop = FALSE]
}

benchmark_consensus <- function(prepared, cutoff, simulation) {
  exclude <- setdiff(dar::steps_ids(prepared, "da"), benchmark_panel())
  baked <- dar::bake(prepared, count_cutoff = cutoff, weights = NULL,
                     exclude = exclude, id = paste0("consensus_", cutoff))
  selected <- as.data.frame(dar::cool(baked))
  expected <- benchmark_explicit_consensus(prepared, cutoff)
  columns <- c("taxa_id", "contrast_id", "effect")
  signature <- function(data) {
    if (!nrow(data)) return(character())
    sort(do.call(paste, c(data[columns], sep = "\r")))
  }
  if (!identical(signature(selected), signature(expected))) {
    stop("Actual bake/cool differs from the explicit directional check.", call. = FALSE)
  }
  results <- data.frame(
    taxa_id = as.character(selected$taxa_id), contrast_id = as.character(selected$contrast_id),
    effect_size = as.numeric(ifelse(selected$effect == "up", 1, -1)),
    significant = rep(TRUE, nrow(selected)), stringsAsFactors = FALSE
  )
  # +/-1 here encodes direction for scoring; it is not a combined effect estimate.
  benchmark_validate_calls(results, simulation)
  list(baked = baked, selected = selected, results = results, checked = TRUE)
}
