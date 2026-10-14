benchmark_parse_arguments <- function(args) {
  values <- list(mode = "preflight", output_dir = "", confirm_pilot = "false", preflight_dir = "")
  seen <- character()
  for (argument in args) {
    if (!grepl("^--[a-z_]+=.+$", argument)) {
      stop("Use --name=value syntax; empty/positional arguments are not allowed.", call. = FALSE)
    }
    name <- sub("^--([^=]+)=.*$", "\\1", argument)
    value <- sub("^--[^=]+=", "", argument)
    if (!name %in% names(values) || name %in% seen) {
      stop("Unknown or duplicate argument: ", name, call. = FALSE)
    }
    seen <- c(seen, name); values[[name]] <- value
  }
  if (!values$mode %in% c("preflight", "pilot") ||
      !values$confirm_pilot %in% c("true", "false")) {
    stop("Invalid mode or confirmation; nothing will run.", call. = FALSE)
  }
  values$confirm_pilot <- identical(values$confirm_pilot, "true")
  if (!nzchar(values$output_dir)) stop("An explicit new output_dir is required.", call. = FALSE)
  if (values$mode == "pilot" && (!values$confirm_pilot || !nzchar(values$preflight_dir))) {
    stop("Pilot is locked: explicit confirm_pilot=true and preflight_dir required.", call. = FALSE)
  }
  if (values$mode == "preflight" && (values$confirm_pilot || nzchar(values$preflight_dir))) {
    stop("Pilot authorization arguments are not applicable to preflight.", call. = FALSE)
  }
  values
}

benchmark_job_plan <- function(preflight = TRUE) {
  ledger <- benchmark_seed_ledger(preflight = preflight)
  rows <- list()
  for (i in seq_len(nrow(ledger))) {
    ids <- if (ledger$cell_id[[i]] %in% c("P5", "P6")) benchmark_panel() else benchmark_engines()
    for (attempt in seq_len(if (preflight) 2L else 1L)) {
      for (engine in ids) {
        index <- match(engine, benchmark_engines())
        rows[[length(rows) + 1L]] <- data.frame(
          mode = ledger$mode[[i]], cell_id = ledger$cell_id[[i]],
          replicate = ledger$replicate[[i]], attempt = attempt, engine = engine,
          data_seed = ledger$data_seed[[i]], engine_seed = ledger$data_seed[[i]] + 100000L * index,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

benchmark_run_row <- function(job, simulation, outcome) {
  data.frame(
    mode = job$mode, cell_id = job$cell_id, generator = simulation$specification$generator,
    replicate = job$replicate, attempt = job$attempt, procedure = job$engine,
    status = outcome$status, reason = outcome$reason,
    warnings = paste(outcome$warnings, collapse = " | "),
    elapsed_seconds = outcome$elapsed_seconds, user_seconds = outcome$user_seconds,
    system_seconds = outcome$system_seconds, peak_rss_mb = outcome$peak_rss_mb,
    worker_elapsed_seconds = if (is.null(outcome$worker_elapsed_seconds)) NA_real_ else outcome$worker_elapsed_seconds,
    n_expected = nrow(simulation$truth), n_returned = nrow(outcome$results),
    stringsAsFactors = FALSE
  )
}

benchmark_dataset_consensus <- function(outcomes, simulation, attempt) {
  records <- metrics <- list()
  completed <- all(benchmark_panel() %in% names(outcomes)) &&
    all(vapply(outcomes[benchmark_panel()], function(x) x$status == "success", logical(1)))
  prepared <- if (completed) tryCatch(
    benchmark_combine_prepared(outcomes, simulation, include_wilcox = !simulation$specification$adjust_batch),
    error = function(cnd) cnd
  ) else NULL
  for (cutoff in c(3L, 4L)) {
    procedure <- if (cutoff == 3L) "majority_3of4" else "unanimity_4of4"
    answer <- if (inherits(prepared, "error")) prepared else if (is.null(prepared)) NULL else
      tryCatch(benchmark_consensus(prepared, cutoff, simulation), error = function(cnd) cnd)
    status <- if (!completed) "unavailable" else if (inherits(answer, "error")) "contract_invalid" else "success"
    reason <- if (!completed) "One or more fixed panel engines failed; no reduced-panel consensus." else
      if (inherits(answer, "error")) conditionMessage(answer) else ""
    records[[procedure]] <- data.frame(
      mode = simulation$specification$mode, cell_id = simulation$specification$cell_id,
      generator = simulation$specification$generator, replicate = simulation$manifest$replicate,
      attempt = attempt, procedure = procedure, status = status, reason = reason,
      warnings = "", elapsed_seconds = NA_real_, user_seconds = NA_real_, system_seconds = NA_real_,
      peak_rss_mb = NA_real_, worker_elapsed_seconds = NA_real_,
      n_expected = nrow(simulation$truth),
      n_returned = if (status == "success") nrow(answer$results) else 0L,
      stringsAsFactors = FALSE
    )
    metrics[[procedure]] <- benchmark_score(
      if (status == "success") answer$results else benchmark_empty_results(),
      simulation$truth, procedure, simulation$specification$cell_id,
      simulation$manifest$replicate, status = status
    )
    if (status == "success") {
      saveRDS(answer, file.path(simulation$artifact_dir, paste0(procedure, "__", attempt, ".rds")))
    }
  }
  list(runs = dplyr::bind_rows(records), metrics = dplyr::bind_rows(metrics))
}

benchmark_write_report <- function(output_dir, runs, checks, mode, decision = NULL,
                                   method_profile = NULL) {
  status <- if (mode == "preflight") if (all(checks$passed)) "PASS" else "FAIL" else "PILOT"
  lines <- c("# DAR benchmark: technical execution", "",
    paste0("Mode: **", mode, "**. Status: **", status, "**."), "",
    if (!is.null(method_profile)) c(
      paste0("Protocol revision: **", method_profile$protocol_revision, "**. ",
             "ANCOM-BC2: `pseudo_sens = TRUE`, native `diff_robust`. ",
             "LinDA: `", method_profile$linda$implementation$implementation, "`."), ""
    ),
    if (mode == "preflight") "No full-pilot datasets were fitted. These checks do not establish scientific performance or FDR control." else
      "Practical exploratory pilot criteria are not statistical guarantees.", "",
    "## Native attempts", "", "| Dataset | Attempt | Engine | Status | Fit seconds | Sampled peak RSS MiB |",
    "|---|---:|---|---|---:|---:|")
  native <- runs[runs$procedure %in% benchmark_engines(), ]
  for (i in seq_len(nrow(native))) {
    row <- native[i, ]
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %.2f | %.1f |", row$cell_id,
      row$attempt, row$procedure, row$status, row$elapsed_seconds, row$peak_rss_mb))
  }
  if (nrow(checks)) {
    lines <- c(lines, "", "## Checks", "")
    for (i in seq_len(nrow(checks))) lines <- c(lines, paste0(
      "- ", if (checks$passed[[i]]) "PASS" else "FAIL", ": ", checks$check[[i]], ". ", checks$detail[[i]]))
  }
  failed <- unique(runs[runs$status != "success", c("procedure", "reason")])
  if (nrow(failed)) {
    lines <- c(lines, "", "## Failures", "")
    for (i in seq_len(nrow(failed))) lines <- c(lines, paste0(
      "- ", failed$procedure[[i]], ": ", gsub("\n", " ", failed$reason[[i]], fixed = TRUE)))
  }
  lines <- c(lines, "", "## Artifacts", "",
    "`runs.csv`, `metrics.csv`, `checks.csv`, the frozen contract/seed ledger, and per-fit RDS/logs retain successful and failed attempts.",
    "", "The full pilot remains a separate, explicit action. No automatic expansion is performed.")
  if (!is.null(decision)) lines <- c(lines, "", paste("Scientific decision:", decision$decision))
  writeLines(lines, file.path(output_dir, "report.md"))
}

benchmark_run <- function(arguments, code_root, repository_root) {
  if (dir.exists(arguments$output_dir) &&
      length(list.files(arguments$output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Refusing to overwrite a nonempty output directory.", call. = FALSE)
  }
  dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(arguments$output_dir, mustWork = TRUE)
  library <- benchmark_prepare_runtime(repository_root, output_dir)
  reservoir <- benchmark_prepare_reservoir(file.path(output_dir, "reservoir"))
  frozen <- benchmark_contract(repository_root, code_root, reservoir)
  benchmark_check_pilot_authority(arguments, frozen$fingerprint)
  ram_limit <- benchmark_ram_limit()
  saveRDS(frozen, file.path(output_dir, "contract.rds"))
  benchmark_write_json(frozen, file.path(output_dir, "contract.json"))
  utils::write.csv(frozen$contract$packages, file.path(output_dir, "versions.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
  provenance <- list(
    git_commit = unname(system2("git", c("-C", shQuote(repository_root), "rev-parse", "HEAD"), stdout = TRUE)),
    installed_dar = find.package("dar"), isolated_library = library,
    dar_install_build = utils::packageDescription("dar")$Built,
    method_profile = frozen$contract$method_profile,
    ram_limit_mb = ram_limit,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), mode = arguments$mode
  )
  benchmark_write_json(provenance, file.path(output_dir, "provenance.json"))
  preflight <- arguments$mode == "preflight"
  jobs <- benchmark_job_plan(preflight = preflight)
  if ((preflight && nrow(jobs) != 20L) || (!preflight && nrow(jobs) != 860L)) {
    stop("Frozen fit budget disagrees with the protocol.", call. = FALSE)
  }
  jobs$status <- "planned"
  ledger_file <- file.path(output_dir, "fit-ledger.csv")
  utils::write.csv(jobs, ledger_file, row.names = FALSE)
  simulations <- outcomes <- all_runs <- all_metrics <- list()
  resource_paused <- FALSE
  for (i in seq_len(nrow(jobs))) {
    job <- as.list(jobs[i, ])
    dataset_key <- paste(job$cell_id, job$replicate, sep = "__")
    if (is.null(simulations[[dataset_key]])) {
      simulation <- benchmark_simulate(job$cell_id, job$replicate,
        reservoir = reservoir, preflight = preflight)
      simulation$artifact_dir <- file.path(output_dir, "datasets", dataset_key)
      dir.create(simulation$artifact_dir, recursive = TRUE)
      saveRDS(simulation, file.path(simulation$artifact_dir, "simulation.rds"))
      utils::write.csv(simulation$truth, file.path(simulation$artifact_dir, "truth.csv"), row.names = FALSE)
      benchmark_write_json(simulation$manifest, file.path(simulation$artifact_dir, "manifest.json"))
      simulations[[dataset_key]] <- simulation
    }
    simulation <- simulations[[dataset_key]]
    fit_key <- paste(job$cell_id, job$replicate, job$engine, job$attempt, sep = "__")
    directory <- file.path(output_dir, "fits", fit_key)
    dir.create(directory, recursive = TRUE)
    job$simulation_file <- file.path(simulation$artifact_dir, "simulation.rds")
    job$outcome_file <- file.path(directory, "outcome.rds")
    job$stdout_file <- file.path(directory, "stdout.log")
    job$stderr_file <- file.path(directory, "stderr.log")
    saveRDS(job, file.path(directory, "job.rds"))
    jobs$status[[i]] <- "started"
    utils::write.csv(jobs, ledger_file, row.names = FALSE)
    message(sprintf("[%d/%d] %s %s attempt %d", i, nrow(jobs), job$cell_id, job$engine, job$attempt))
    outcome <- benchmark_run_engine(job, code_root, ram_limit_mb = ram_limit)
    jobs$status[[i]] <- outcome$status
    utils::write.csv(jobs, ledger_file, row.names = FALSE)
    key <- paste(job$cell_id, job$replicate, job$engine, job$attempt, sep = "__")
    # Store paths, not all captured native model objects, across the full pilot.
    outcomes[[key]] <- job$outcome_file
    all_runs[[i]] <- benchmark_run_row(job, simulation, outcome)
    metric <- benchmark_score(outcome$results, simulation$truth, job$engine,
      job$cell_id, job$replicate, status = outcome$status)
    metric$attempt <- job$attempt; metric$mode <- job$mode
    all_metrics[[i]] <- metric
    utils::write.csv(dplyr::bind_rows(all_runs), file.path(output_dir, "runs.csv"), row.names = FALSE)
    if (preflight) message("  ", outcome$status, if (nzchar(outcome$reason)) paste0(": ", outcome$reason) else "")
    observed <- dplyr::bind_rows(all_runs)
    last_engine <- if (simulation$specification$adjust_batch) tail(benchmark_panel(), 1L) else tail(benchmark_engines(), 1L)
    projection <- NA_real_
    if (!preflight && job$replicate == 5L && job$engine == last_engine) {
      projection <- benchmark_cost_projection(observed, jobs)
      benchmark_write_json(list(cell_id = job$cell_id, projected_cpu_hours = projection,
        estimator = "Observed cell-engine means; unmeasured cells provisionally use their engine mean."),
        file.path(output_dir, paste0("cost-projection-", job$cell_id, ".json")))
    }
    spent <- sum(ifelse(is.finite(observed$user_seconds + observed$system_seconds),
                        observed$user_seconds + observed$system_seconds,
                        observed$worker_elapsed_seconds), na.rm = TRUE) / 3600
    if (!preflight && (isTRUE(outcome$resource_stop) || spent > 24 ||
                       (is.finite(projection) && projection > 24))) {
      resource_paused <- TRUE
      message("Resource ceiling reached/projected; remaining fits will not start.")
      break
    }
  }
  consensus_rows <- list()
  for (key in names(simulations)) {
    simulation <- simulations[[key]]
    for (attempt in seq_len(if (preflight) 2L else 1L)) {
      ids <- if (simulation$specification$adjust_batch) benchmark_panel() else benchmark_engines()
      current <- stats::setNames(lapply(ids, function(id) {
        fit_key <- paste(simulation$specification$cell_id, simulation$manifest$replicate,
                         id, attempt, sep = "__")
        if (is.null(outcomes[[fit_key]])) {
          return(list(status = "unavailable", reason = "Fit not attempted."))
        }
        readRDS(outcomes[[fit_key]])
      }), ids)
      consensus <- benchmark_dataset_consensus(current, simulation, attempt)
      consensus_rows[[length(consensus_rows) + 1L]] <- consensus$runs
      consensus$metrics$attempt <- attempt; consensus$metrics$mode <- arguments$mode
      all_metrics[[length(all_metrics) + 1L]] <- consensus$metrics
    }
  }
  runs <- dplyr::bind_rows(c(all_runs, consensus_rows))
  metrics <- dplyr::bind_rows(all_metrics)
  utils::write.csv(runs, file.path(output_dir, "runs.csv"), row.names = FALSE)
  utils::write.csv(metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
  if (preflight) {
    read_outcomes <- lapply(outcomes, readRDS)
    checks <- benchmark_preflight_checks(read_outcomes, dplyr::bind_rows(consensus_rows))
    for (key in names(read_outcomes)) {
      extra <- read_outcomes[[key]]$checks
      if (nrow(extra)) {
        extra$check <- paste(key, extra$check, sep = "__")
        checks <- rbind(checks, extra)
      }
    }
    receipt <- list(mode = "preflight", passed = all(checks$passed),
      fingerprint = frozen$fingerprint, fits_attempted = nrow(jobs),
      method_profile = frozen$contract$method_profile,
      completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    saveRDS(receipt, file.path(output_dir, "preflight-receipt.rds"))
    benchmark_write_json(receipt, file.path(output_dir, "preflight-receipt.json"))
    decision <- list(decision = "not_applicable", reason = "Preflight is technical, not scientific evidence.")
  } else {
    checks <- data.frame(check = character(), passed = logical(), detail = character())
    summary <- benchmark_summarize(metrics)
    utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
    # Scientific interpretation of discordance remains a documented human input,
    # not an automatic success inferred from a change in truth denominators.
    diagnostics <- list(mode = "pilot", unexplained_discrepancies = sum(runs$status == "contract_invalid"),
      hypotheses_interpretable = all(summary$max_estimand_shift[summary$cell_id %in% c("P2", "P6", "P8")] < 1e-10),
      failures_documented = all(runs$status == "success" | nzchar(runs$reason)),
      projected_cpu_hours = sum(runs$user_seconds + runs$system_seconds, na.rm = TRUE) / 3600,
      discordant_taxa = data.frame(), uncertainty_too_wide = TRUE)
    decision <- benchmark_decide(summary, runs, diagnostics)
  }
  utils::write.csv(checks, file.path(output_dir, "checks.csv"), row.names = FALSE)
  benchmark_write_json(decision, file.path(output_dir, "decision.json"))
  saveRDS(list(runs = runs, metrics = metrics, checks = checks, decision = decision,
               method_profile = frozen$contract$method_profile),
          file.path(output_dir, "benchmark-results.rds"))
  benchmark_write_report(output_dir, runs, checks, arguments$mode, decision,
                         method_profile = frozen$contract$method_profile)
  list(status = if (resource_paused) 2L else if (preflight && !all(checks$passed)) 1L else 0L,
       output_dir = output_dir, checks = checks, runs = runs, decision = decision)
}
