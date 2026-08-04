parse_validation_arguments <- function(args) {
  defaults <- list(
    profile = "quick",
    engine = "all",
    scenario = "all",
    output_dir = file.path("validation-results", "quick", "all"),
    base_seed = "140",
    strict = "true"
  )
  for (argument in args) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", argument, call. = FALSE)
    }
    pair <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    name <- pair[[1L]]
    value <- paste(pair[-1L], collapse = "=")
    if (!name %in% names(defaults)) {
      stop("Unknown validation argument: --", name, call. = FALSE)
    }
    defaults[[name]] <- value
  }
  defaults$base_seed <- as.integer(defaults$base_seed)
  defaults$strict <- tolower(defaults$strict) %in% c("true", "1", "yes")
  defaults
}

validation_selection <- function(value, available, allow_none = FALSE) {
  if (identical(value, "all")) {
    return(available)
  }
  if (allow_none && identical(value, "none")) {
    return(character())
  }
  selected <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(selected, available)
  if (length(unknown)) {
    stop("Unknown selection: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  selected
}

run_scientific_validation <- function(profile_name = "quick", engines = "all",
                                      scenarios = "all", output_dir,
                                      base_seed = 140L, strict = TRUE) {
  profile <- validation_profile(profile_name)
  selected_engines <- validation_selection(
    engines, names(validation_engine_registry()), allow_none = TRUE
  )
  selected_scenarios <- validation_selection(
    scenarios, names(validation_scenario_registry())
  )

  results <- list()
  runs <- list()
  truths <- list()
  manifests <- list()
  result_index <- run_index <- simulation_index <- 0L

  for (scenario_id in selected_scenarios) {
    for (replicate in seq_len(profile$replicates[[scenario_id]])) {
      simulation <- simulate_validation_scenario(
        scenario_id, replicate = replicate, base_seed = base_seed
      )
      simulation_index <- simulation_index + 1L
      truths[[simulation_index]] <- simulation$truth
      manifests[[simulation_index]] <- simulation$manifest
      for (engine in selected_engines) {
        outcome <- run_validation_case(simulation, engine, profile)
        run_index <- run_index + 1L
        runs[[run_index]] <- outcome$run
        if (nrow(outcome$results)) {
          result_index <- result_index + 1L
          results[[result_index]] <- outcome$results
        }
      }
    }
  }

  results <- dplyr::bind_rows(results)
  runs <- dplyr::bind_rows(runs)
  truth <- dplyr::bind_rows(truths)
  manifests <- dplyr::bind_rows(manifests)
  scientific_metrics <- score_validation_results(results, truth, profile$alpha)
  execution_metrics <- score_validation_runs(runs)
  metrics <- dplyr::bind_rows(scientific_metrics, execution_metrics)
  summary <- summarize_validation_metrics(metrics)
  gates <- evaluate_validation_gates(summary)
  versions <- validation_version_manifest(
    selected_engines, profile$name, base_seed
  )
  write_validation_artifacts(
    output_dir, results, runs, truth, manifests, metrics, summary, gates,
    versions
  )

  engine_failures <- nrow(runs) > 0L && any(runs$status == "engine_error")
  gate_failures <- nrow(gates) > 0L && any(gates$status == "fail")
  self_test_passed <- direction_reversal_self_test()
  status <- if (strict && (engine_failures || gate_failures || !self_test_passed)) 1L else 0L
  list(
    status = status,
    engine_failures = engine_failures,
    gate_failures = gate_failures,
    direction_self_test = self_test_passed,
    results = results,
    runs = runs,
    truth = truth,
    manifests = manifests,
    metrics = metrics,
    summary = summary,
    gates = gates,
    versions = versions
  )
}
