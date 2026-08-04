read_validation_artifact <- function(input_dir, filename) {
  paths <- list.files(
    input_dir,
    pattern = paste0("^", gsub("[.]", "[.]", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  tables <- lapply(paths, function(path) {
    tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(cnd) data.frame()
    )
  })
  dplyr::bind_rows(tables)
}

aggregate_validation_artifacts <- function(input_dir, output_dir,
                                           render = TRUE, strict = TRUE,
                                           baseline = read_validation_baseline(),
                                           enforce_baseline = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  artifact_names <- c(
    "results", "runs", "truth", "manifests", "metrics", "summary", "gates",
    "versions", "expected"
  )
  artifacts <- stats::setNames(lapply(artifact_names, function(name) {
    read_validation_artifact(input_dir, paste0(name, ".csv"))
  }), artifact_names)
  artifacts$truth <- dplyr::distinct(artifacts$truth)
  artifacts$manifests <- dplyr::distinct(artifacts$manifests)
  artifacts$versions <- dplyr::distinct(artifacts$versions)
  artifacts$execution_status <- summarize_validation_execution(
    artifacts$runs, artifacts$expected
  )
  baseline_for_run <- baseline
  if (nrow(baseline_for_run) && nrow(artifacts$expected)) {
    baseline_for_run <- dplyr::semi_join(
      baseline_for_run,
      unique(artifacts$expected[c("engine", "scenario")]),
      by = c("engine", "scenario")
    )
  }
  artifacts$baseline_comparison <- compare_validation_baseline(
    artifacts$summary, baseline_for_run
  )
  if (!enforce_baseline && nrow(artifacts$baseline_comparison)) {
    artifacts$baseline_comparison$status <- "informational"
    artifacts$baseline_comparison$reason <-
      "baseline comparison is informational outside the full profile"
  }
  artifact_names <- c(
    artifact_names, "execution_status", "baseline_comparison"
  )
  for (name in artifact_names) {
    utils::write.csv(
      artifacts[[name]], file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE, na = ""
    )
  }
  saveRDS(artifacts, file.path(output_dir, "validation-results.rds"))

  if (render) {
    rmarkdown::render(
      file.path(validation_root_path(), "report", "scientific-validation.Rmd"),
      output_file = "scientific-validation.html",
      output_dir = normalizePath(output_dir),
      params = list(data_dir = normalizePath(output_dir)),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
  }

  gate_failure <- nrow(artifacts$gates) > 0L &&
    any(artifacts$gates$status == "fail", na.rm = TRUE)
  engine_failure <- nrow(artifacts$runs) > 0L &&
    any(artifacts$runs$status == "engine_error", na.rm = TRUE)
  execution_failure <- nrow(artifacts$execution_status) > 0L &&
    any(grepl("^failed_|^missing_", artifacts$execution_status$status))
  baseline_failure <- nrow(artifacts$baseline_comparison) > 0L &&
    any(artifacts$baseline_comparison$status == "fail", na.rm = TRUE)
  list(
    status = if (strict && (
      gate_failure || engine_failure || execution_failure || baseline_failure
    )) 1L else 0L,
    gate_failure = gate_failure,
    engine_failure = engine_failure,
    execution_failure = execution_failure,
    baseline_failure = baseline_failure,
    artifacts = artifacts
  )
}

summarize_validation_execution <- function(runs, expected = data.frame()) {
  if (!nrow(expected)) {
    if (!nrow(runs)) return(data.frame())
    expected <- unique(runs[c("engine", "scenario")])
  }
  expected <- unique(expected[c("engine", "scenario")])
  rows <- lapply(seq_len(nrow(expected)), function(index) {
    item <- expected[index, , drop = FALSE]
    observed <- runs[
      runs$engine == item$engine & runs$scenario == item$scenario,
      , drop = FALSE
    ]
    statuses <- unique(observed$status)
    status <- if (!nrow(observed)) {
      "missing_artifact"
    } else if ("engine_error" %in% statuses) {
      "failed_engine"
    } else if ("dependency_missing" %in% statuses) {
      "failed_dependency"
    } else if (all(statuses == "expected_skip")) {
      "skipped_incompatible"
    } else if ("success" %in% statuses) {
      "success"
    } else {
      "failed_unknown"
    }
    data.frame(
      engine = item$engine,
      scenario = item$scenario,
      status = status,
      runs = nrow(observed),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}
