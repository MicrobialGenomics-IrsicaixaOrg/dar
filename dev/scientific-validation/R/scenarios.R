validation_scenario_registry <- function() {
  list(
    cross_sectional_null = list(
      formula = ~ condition,
      time = NULL,
      contrasts = "simple",
      na_action = "error",
      design = "cross_sectional",
      signal = "null",
      confounders = FALSE,
      time_interaction = FALSE,
      random = FALSE,
      sparse = FALSE,
      missing = FALSE
    ),
    cross_sectional_signal = list(
      formula = ~ condition,
      time = NULL,
      contrasts = "simple",
      na_action = "error",
      design = "cross_sectional",
      signal = "condition",
      confounders = FALSE,
      time_interaction = FALSE,
      random = FALSE,
      sparse = FALSE,
      missing = FALSE
    ),
    sparse_zero_inflated = list(
      formula = ~ condition,
      time = NULL,
      contrasts = "simple",
      na_action = "error",
      design = "cross_sectional",
      signal = "condition",
      confounders = FALSE,
      time_interaction = FALSE,
      random = FALSE,
      sparse = TRUE,
      missing = FALSE
    ),
    fixed_confounder = list(
      formula = ~ condition + batch,
      time = NULL,
      contrasts = "simple",
      na_action = "error",
      design = "confounded",
      signal = "condition",
      confounders = TRUE,
      time_interaction = FALSE,
      random = FALSE,
      sparse = FALSE,
      missing = FALSE
    ),
    longitudinal_interaction = list(
      formula = ~ condition * time + batch,
      time = "time",
      contrasts = "all",
      na_action = "error",
      design = "longitudinal",
      signal = "longitudinal",
      confounders = TRUE,
      time_interaction = TRUE,
      random = FALSE,
      sparse = FALSE,
      missing = FALSE
    ),
    longitudinal_random_effect = list(
      formula = ~ condition * time + batch + (time | subject),
      time = "time",
      contrasts = "all",
      na_action = "error",
      design = "longitudinal",
      signal = "longitudinal",
      confounders = TRUE,
      time_interaction = TRUE,
      random = TRUE,
      sparse = FALSE,
      missing = FALSE
    ),
    missing_unbalanced = list(
      formula = ~ condition * time + batch,
      time = "time",
      contrasts = "all",
      na_action = "drop",
      design = "longitudinal_unbalanced",
      signal = "longitudinal",
      confounders = TRUE,
      time_interaction = TRUE,
      random = FALSE,
      sparse = FALSE,
      missing = TRUE
    )
  )
}

validation_scenario <- function(id) {
  scenarios <- validation_scenario_registry()
  if (!id %in% names(scenarios)) {
    stop("Unknown validation scenario: ", id, call. = FALSE)
  }
  c(list(id = id), scenarios[[id]])
}

validation_root_path <- function() {
  root <- getOption("dar.validation.root", Sys.getenv("DAR_VALIDATION_ROOT", ""))
  if (!nzchar(root)) {
    root <- file.path("dev", "scientific-validation")
  }
  normalizePath(root, mustWork = TRUE)
}

validation_profile <- function(profile = c("quick", "full")) {
  profile <- match.arg(profile)
  config <- utils::read.csv(
    file.path(validation_root_path(), "profiles.csv"),
    stringsAsFactors = FALSE
  )
  config <- config[config$profile == profile, , drop = FALSE]
  scenarios <- names(validation_scenario_registry())
  if (!setequal(config$scenario, scenarios)) {
    stop("The profile does not configure every validation scenario.", call. = FALSE)
  }
  replicates <- stats::setNames(as.integer(config$replicates), config$scenario)
  list(
    name = profile,
    replicates = replicates,
    aldex_mc_samples = unique(as.integer(config$aldex_mc_samples)),
    alpha = unique(config$alpha)
  )
}

validation_seed <- function(scenario_id, replicate, base_seed = 140L) {
  scenario_index <- match(scenario_id, names(validation_scenario_registry()))
  if (is.na(scenario_index) || length(replicate) != 1L || replicate < 1L) {
    stop("Invalid scenario or replicate for seed generation.", call. = FALSE)
  }
  as.integer(base_seed + scenario_index * 100000L + replicate)
}
