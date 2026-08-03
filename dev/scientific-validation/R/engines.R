validation_engine_registry <- function() {
  list(
    deseq = list(
      packages = c("DESeq2", "SummarizedExperiment"),
      confounders = TRUE, time_interaction = TRUE, random = FALSE,
      truth_multiplier = 1 / log(2), comparable = TRUE
    ),
    aldex = list(
      packages = "ALDEx2",
      confounders = TRUE, time_interaction = TRUE, random = FALSE,
      truth_multiplier = NA_real_, comparable = FALSE
    ),
    ancom = list(
      packages = c("ANCOMBC", "mia", "microbiome"),
      confounders = TRUE, time_interaction = TRUE, random = TRUE,
      truth_multiplier = 1, comparable = TRUE
    ),
    corncob = list(
      packages = "corncob",
      confounders = TRUE, time_interaction = TRUE, random = FALSE,
      truth_multiplier = NA_real_, comparable = FALSE
    ),
    linda = list(
      packages = "MicrobiomeStat",
      confounders = TRUE, time_interaction = TRUE, random = TRUE,
      truth_multiplier = 1 / log(2), comparable = TRUE
    ),
    maaslin = list(
      packages = "maaslin3",
      confounders = TRUE, time_interaction = TRUE, random = TRUE,
      truth_multiplier = NA_real_, comparable = FALSE
    ),
    wilcox = list(
      packages = c("microbiome", "rstatix"),
      confounders = FALSE, time_interaction = FALSE, random = FALSE,
      truth_multiplier = NA_real_, comparable = FALSE
    ),
    lefse = list(
      packages = c("lefser", "SummarizedExperiment"),
      confounders = FALSE, time_interaction = FALSE, random = FALSE,
      truth_multiplier = NA_real_, comparable = FALSE
    )
  )
}

validation_engine <- function(engine) {
  registry <- validation_engine_registry()
  if (!engine %in% names(registry)) {
    stop("Unknown validation engine: ", engine, call. = FALSE)
  }
  c(list(id = engine), registry[[engine]])
}

validation_engine_compatible <- function(spec, engine) {
  config <- validation_engine(engine)
  reasons <- character()
  if (spec$confounders && !config$confounders) {
    reasons <- c(reasons, "does not support adjusted confounders")
  }
  if (spec$time_interaction && !config$time_interaction) {
    reasons <- c(reasons, "does not support condition-time interactions")
  }
  if (spec$random && !config$random) {
    reasons <- c(reasons, "does not support random effects")
  }
  list(
    compatible = length(reasons) == 0L,
    reason = if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_
  )
}

validation_missing_packages <- function(engine) {
  config <- validation_engine(engine)
  config$packages[!vapply(config$packages, requireNamespace, logical(1), quietly = TRUE)]
}

add_validation_engine <- function(rec, engine, profile) {
  id <- paste0(engine, "__scientific_validation")
  switch(
    engine,
    deseq = dar::step_deseq(
      rec, fitType = "mean", type = "normal", id = id
    ),
    aldex = dar::step_aldex(
      rec, mc.samples = profile$aldex_mc_samples, id = id
    ),
    ancom = dar::step_ancom(
      rec, p_adj_method = "BH", prv_cut = 0, n_cl = 1, id = id
    ),
    corncob = dar::step_corncob(
      rec, filter_discriminant = FALSE, id = id
    ),
    linda = dar::step_linda(
      rec, prev_filter = 0, winsorize = FALSE, adaptive = FALSE,
      id = id
    ),
    maaslin = dar::step_maaslin(
      rec, min_prevalence = 0, standardize = FALSE, id = id
    ),
    wilcox = dar::step_wilcox(rec, id = id),
    lefse = dar::step_lefse(
      rec, lda.threshold = 0, rarefy = FALSE, id = id
    ),
    stop("Unsupported validation engine: ", engine, call. = FALSE)
  )
}

validation_numeric_column <- function(result, candidates) {
  candidate <- candidates[candidates %in% names(result)][1L]
  if (length(candidate) == 0L || is.na(candidate)) {
    return(rep(NA_real_, nrow(result)))
  }
  as.numeric(result[[candidate]])
}

normalize_validation_result <- function(result, engine, simulation,
                                        raw_result = NULL) {
  required <- c(
    "taxa_id", "contrast_id", "effect_size", "effect_metric", "adj_p_value"
  )
  if (!all(required %in% names(result))) {
    stop(
      "Engine result is missing required modeled columns: ",
      paste(setdiff(required, names(result)), collapse = ", "),
      call. = FALSE
    )
  }
  keys <- result[c("taxa_id", "contrast_id")]
  if (anyDuplicated(keys)) {
    stop("Engine result contains duplicated taxon-contrast keys.", call. = FALSE)
  }
  config <- validation_engine(engine)
  standard_error <- rep(NA_real_, nrow(result))
  if (!is.null(raw_result)) {
    raw_result <- dar:::flatten_model_output(raw_result)
    raw_standard_error <- validation_numeric_column(
      raw_result, c("lfcSE", "stderr", "std_error", "se")
    )
    raw_keys <- data.frame(
      taxa_id = as.character(raw_result$taxa_id),
      contrast_id = as.character(raw_result$contrast_id),
      std_error = raw_standard_error,
      stringsAsFactors = FALSE
    )
    matched <- match(
      paste(result$taxa_id, result$contrast_id, sep = "\r"),
      paste(raw_keys$taxa_id, raw_keys$contrast_id, sep = "\r")
    )
    standard_error <- raw_keys$std_error[matched]
  }
  data.frame(
    engine = engine,
    scenario = simulation$manifest$scenario,
    replicate = simulation$manifest$replicate,
    taxa_id = as.character(result$taxa_id),
    contrast_id = as.character(result$contrast_id),
    effect = as.numeric(result$effect_size),
    padj = as.numeric(result$adj_p_value),
    std_error = standard_error,
    effect_metric = as.character(result$effect_metric),
    truth_multiplier = config$truth_multiplier,
    comparable = config$comparable,
    stringsAsFactors = FALSE
  )
}

validation_run_record <- function(engine, simulation, status, reason = NA_character_,
                                  warnings = character(), duration = 0) {
  data.frame(
    engine = engine,
    scenario = simulation$manifest$scenario,
    replicate = simulation$manifest$replicate,
    status = status,
    reason = reason,
    warnings = paste(unique(warnings), collapse = " | "),
    duration_seconds = as.numeric(duration),
    stringsAsFactors = FALSE
  )
}

run_validation_case <- function(simulation, engine, profile) {
  compatibility <- validation_engine_compatible(simulation$specification, engine)
  if (!compatibility$compatible) {
    return(list(
      results = data.frame(),
      run = validation_run_record(
        engine, simulation, "expected_skip", compatibility$reason
      )
    ))
  }
  missing <- validation_missing_packages(engine)
  if (length(missing)) {
    return(list(
      results = data.frame(),
      run = validation_run_record(
        engine, simulation, "dependency_missing",
        paste("Missing optional packages:", paste(missing, collapse = ", "))
      )
    ))
  }

  warnings <- character()
  started <- proc.time()[["elapsed"]]
  outcome <- tryCatch(
    withCallingHandlers(
      {
        rec <- add_validation_engine(simulation$recipe, engine, profile)
        set.seed(simulation$manifest$seed + match(
          engine, names(validation_engine_registry())
        ) * 1000L)
        prepared <- dar::prep(rec, parallel = FALSE)
        step_id <- paste0(engine, "__scientific_validation")
        normalize_validation_result(
          dar::tidy_results(prepared, steps = step_id),
          engine,
          simulation,
          raw_result = prepared@results[[step_id]]
        )
      },
      warning = function(cnd) {
        warnings <<- c(warnings, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(cnd) cnd
  )
  duration <- proc.time()[["elapsed"]] - started
  if (inherits(outcome, "error")) {
    return(list(
      results = data.frame(),
      run = validation_run_record(
        engine, simulation, "engine_error",
        paste(class(outcome)[[1L]], conditionMessage(outcome), sep = ": "),
        warnings, duration
      )
    ))
  }
  list(
    results = outcome,
    run = validation_run_record(
      engine, simulation, "success", warnings = warnings, duration = duration
    )
  )
}
