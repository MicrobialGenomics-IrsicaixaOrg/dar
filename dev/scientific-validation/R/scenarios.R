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

validation_profile <- function(profile = c("quick", "full")) {
  profile <- match.arg(profile)
  scenarios <- names(validation_scenario_registry())
  replicates <- if (identical(profile, "quick")) {
    stats::setNames(rep(2L, length(scenarios)), scenarios)
  } else {
    c(
      cross_sectional_null = 100L,
      cross_sectional_signal = 50L,
      sparse_zero_inflated = 40L,
      fixed_confounder = 40L,
      longitudinal_interaction = 30L,
      longitudinal_random_effect = 30L,
      missing_unbalanced = 30L
    )
  }
  list(
    name = profile,
    replicates = replicates,
    aldex_mc_samples = if (identical(profile, "quick")) 8L else 32L,
    alpha = 0.05
  )
}

validation_seed <- function(scenario_id, replicate, base_seed = 140L) {
  scenario_index <- match(scenario_id, names(validation_scenario_registry()))
  if (is.na(scenario_index) || length(replicate) != 1L || replicate < 1L) {
    stop("Invalid scenario or replicate for seed generation.", call. = FALSE)
  }
  as.integer(base_seed + scenario_index * 100000L + replicate)
}
