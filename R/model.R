# CENTRAL STATISTICAL MODEL ---------------------------------------------------

#' Get the centralized statistical model from a recipe
#'
#' @param rec A [Recipe()] or [PrepRecipe()] object.
#'
#' @return `NULL` when no model is defined, otherwise a list containing the
#'   formula and its contrast, reference, missing-value and time controls.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   add_model(~ RiskGroup2)
#' get_model(rec)
get_model <- function(rec) {
  check_any_recipe(rec)
  rec@model
}

#' Add a centralized statistical model to a recipe
#'
#' The model is the single source of truth for statistical design arguments.
#' Variables in `var_info` are the targets of inference; other fixed terms are
#' adjustment variables. When `time` is supplied, numeric values are treated as
#' ordered categorical time points while compiling the model.
#'
#' @param rec An unprepared [Recipe()] object.
#' @param formula A one-sided formula or `NULL` to remove the model.
#' @param time Optional single metadata column identifying time.
#' @param reference Optional named character vector defining factor reference
#'   levels.
#' @param contrasts Either `"simple"` for pairwise/simple effects or `"all"`
#'   to additionally request difference-in-differences contrasts.
#' @param na_action Either `"error"` or `"drop"`. Dropping is performed once
#'   for all engines during [prep()].
#' @param max_time_levels Maximum number of time points that may be created by
#'   automatic factorization.
#'
#' @return The updated `Recipe`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   add_model(~ RiskGroup2 + Cluster)
#' get_model(rec)
add_model <- function(rec,
                      formula,
                      time = NULL,
                      reference = NULL,
                      contrasts = c("simple", "all"),
                      na_action = c("error", "drop"),
                      max_time_levels = 12) {
  check_recipe(rec)

  if (is.null(formula)) {
    rec@model <- NULL
    validate_recipe_object(rec)
    return(rec)
  }

  contrasts <- match.arg(contrasts)
  na_action <- match.arg(na_action)
  formula <- coerce_model_formula(formula)

  candidate <- list(
    formula = formula,
    time = time,
    reference = reference,
    contrasts = contrasts,
    na_action = na_action,
    max_time_levels = as.integer(max_time_levels)
  )

  validate_model_spec(rec, candidate)
  rec@model <- candidate
  validate_recipe_object(rec)
  report_model_step_warnings(rec)
  rec
}

#' @noRd
coerce_model_formula <- function(formula) {
  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    formula <- trimws(formula)
    if (!startsWith(formula, "~")) {
      formula <- paste("~", formula)
    }
    formula <- stats::as.formula(formula, env = parent.frame())
  }

  if (!inherits(formula, "formula") || length(formula) != 2L) {
    cli::cli_abort(
      c(
        "x" = "{.arg formula} must be a one-sided formula.",
        "i" = "For example: {.code ~ condition * time + batch + (1 | subject)}."
      ),
      class = "dar_error_invalid_model"
    )
  }
  formula
}

#' @noRd
model_metadata <- function(rec) {
  rec@phyloseq |>
    phyloseq::sample_data() |>
    as("data.frame") |>
    tibble::rownames_to_column("sample_id")
}

#' Return model validity problems without throwing
#' @noRd
model_validity_problems <- function(object) {
  model <- object@model
  if (is.null(model)) {
    return(character())
  }

  required <- c(
    "formula", "time", "reference", "contrasts", "na_action",
    "max_time_levels"
  )
  problems <- character()

  if (!identical(names(model), required)) {
    problems <- c(
      problems,
      paste0("`model` must contain exactly: ", paste(required, collapse = ", "))
    )
  }
  if (!inherits(model$formula, "formula") || length(model$formula) != 2L) {
    problems <- c(problems, "`model$formula` must be a one-sided formula")
  }
  if (!is.null(model$time) &&
      (!is.character(model$time) || length(model$time) != 1L || is.na(model$time))) {
    problems <- c(problems, "`model$time` must be NULL or one column name")
  }
  if (!is.null(model$reference) &&
      (!is.character(model$reference) || is.null(names(model$reference)) ||
       any(names(model$reference) == ""))) {
    problems <- c(problems, "`model$reference` must be a named character vector")
  }
  if (!identical(model$contrasts, "simple") && !identical(model$contrasts, "all")) {
    problems <- c(problems, "`model$contrasts` must be `simple` or `all`")
  }
  if (!identical(model$na_action, "error") && !identical(model$na_action, "drop")) {
    problems <- c(problems, "`model$na_action` must be `error` or `drop`")
  }
  if (length(model$max_time_levels) != 1L || is.na(model$max_time_levels) ||
      model$max_time_levels < 2L) {
    problems <- c(problems, "`model$max_time_levels` must be at least 2")
  }

  if (length(problems) == 0L) {
    semantic_problem <- tryCatch(
      {
        validate_model_spec(object, model)
        NULL
      },
      error = conditionMessage
    )
    if (!is.null(semantic_problem)) {
      problems <- c(problems, semantic_problem)
    }
  }
  problems
}

#' @noRd
validate_model_spec <- function(rec, model) {
  formula <- model$formula
  metadata <- model_metadata(rec)
  variables <- all.vars(formula)
  missing_vars <- setdiff(variables, names(metadata))

  if (length(missing_vars) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Model variable{?s} {.val {missing_vars}} {?is/are} missing from sample metadata.",
        "i" = "Available columns: {.val {setdiff(names(metadata), 'sample_id')}}."
      ),
      class = c("dar_error_missing_vars", "dar_error_invalid_model")
    )
  }

  if (!is.null(model$time)) {
    if (!is.character(model$time) || length(model$time) != 1L || is.na(model$time)) {
      cli::cli_abort(
        "{.arg time} must be NULL or a single metadata column name.",
        class = "dar_error_invalid_model"
      )
    }
    if (!model$time %in% variables) {
      cli::cli_abort(
        c(
          "x" = "The time variable {.var {model$time}} is not present in {.arg formula}.",
          "i" = "Include it as a fixed term, normally {.code target * time}."
        ),
        class = "dar_error_invalid_model"
      )
    }
  }

  fixed_formula <- reformulas::nobars(formula)
  fixed_terms <- stats::terms(fixed_formula)
  labels <- attr(fixed_terms, "term.labels")
  if (length(labels) == 0L || attr(fixed_terms, "intercept") == 0L) {
    cli::cli_abort(
      "The model must contain an intercept and at least one fixed term.",
      class = "dar_error_invalid_model"
    )
  }

  targets <- rec@var_info$vars
  if (length(targets) > 0L) {
    missing_targets <- setdiff(targets, all.vars(fixed_formula))
    if (length(missing_targets) > 0L) {
      cli::cli_abort(
        c(
          "x" = "Target variable{?s} {.val {missing_targets}} from {.arg var_info} {?is/are} absent from the fixed model.",
          "i" = "Every target must appear in {.arg formula}."
        ),
        class = "dar_error_invalid_model"
      )
    }
  }

  for (label in labels) {
    parts <- strsplit(label, ":", fixed = TRUE)[[1]]
    if (length(parts) == 1L) {
      if (!label %in% names(metadata)) {
        cli::cli_abort(
          c(
            "x" = "Unsupported model term {.code {label}}.",
            "i" = "Only metadata columns, target-time interactions and random effects are supported."
          ),
          class = "dar_error_invalid_model"
        )
      }
    } else if (length(parts) != 2L || is.null(model$time) ||
               !model$time %in% parts ||
               length(intersect(parts, targets)) != 1L) {
      cli::cli_abort(
        c(
          "x" = "Unsupported interaction {.code {label}}.",
          "i" = "Interactions must be two-way interactions between a target and {.arg time}."
        ),
        class = "dar_error_invalid_model"
      )
    }
  }

  random_terms <- reformulas::findbars(formula)
  for (term in random_terms) {
    grouping <- all.vars(term[[3]])
    slopes <- setdiff(all.vars(term[[2]]), "1")
    if (length(grouping) != 1L || !grouping %in% names(metadata) ||
        length(setdiff(slopes, model$time)) > 0L) {
      cli::cli_abort(
        c(
          "x" = "Unsupported random-effect term {.code {paste(deparse(term), collapse = ' ')}}.",
          "i" = "Use a random intercept or a time slope grouped by one metadata column."
        ),
        class = "dar_error_invalid_model"
      )
    }
  }

  if (!is.null(model$reference)) {
    if (!is.character(model$reference) || is.null(names(model$reference)) ||
        any(names(model$reference) == "") || anyDuplicated(names(model$reference))) {
      cli::cli_abort(
        "{.arg reference} must be a named character vector with unique names.",
        class = "dar_error_invalid_model"
      )
    }
    unknown_reference_vars <- setdiff(
      names(model$reference), all.vars(fixed_formula)
    )
    if (length(unknown_reference_vars) > 0L) {
      cli::cli_abort(
        "Reference variable{?s} {.val {unknown_reference_vars}} {?is/are} not in the model.",
        class = "dar_error_invalid_model"
      )
    }
    numeric_references <- names(model$reference)[vapply(
      names(model$reference),
      function(var) is.numeric(metadata[[var]]) && !identical(var, model$time),
      logical(1)
    )]
    if (length(numeric_references) > 0L) {
      cli::cli_abort(
        "Reference levels can only be set for categorical variables or {.arg time}; invalid: {.val {numeric_references}}.",
        class = "dar_error_invalid_model"
      )
    }
    bad_levels <- names(model$reference)[vapply(names(model$reference), function(var) {
      !model$reference[[var]] %in% as.character(stats::na.omit(metadata[[var]]))
    }, logical(1))]
    if (length(bad_levels) > 0L) {
      cli::cli_abort(
        "Invalid reference level for variable{?s} {.val {bad_levels}}.",
        class = "dar_error_invalid_model"
      )
    }
  }

  if (length(model$max_time_levels) != 1L || is.na(model$max_time_levels) ||
      model$max_time_levels < 2L) {
    cli::cli_abort(
      "{.arg max_time_levels} must be a single integer greater than one.",
      class = "dar_error_invalid_model"
    )
  }

  invisible(TRUE)
}

#' @noRd
ordered_model_levels <- function(x, numeric_order = FALSE) {
  x <- x[!is.na(x)]
  if (is.factor(x)) {
    return(levels(droplevels(x)))
  }
  if (numeric_order || is.numeric(x) || is.integer(x)) {
    return(as.character(sort(unique(x))))
  }
  sort(unique(as.character(x)))
}

#' Resolve a model against the current samples
#' @noRd
resolve_model <- function(rec, check_design = TRUE) {
  model <- get_model(rec)
  if (is.null(model)) {
    return(NULL)
  }
  validate_model_spec(rec, model)

  metadata <- model_metadata(rec)
  model_vars <- unique(all.vars(model$formula))
  complete <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
  dropped <- metadata$sample_id[!complete]
  if (length(dropped) > 0L && identical(model$na_action, "error")) {
    cli::cli_abort(
      c(
        "x" = "The model contains missing values in {length(dropped)} sample{?s}.",
        "i" = "Use {.code na_action = 'drop'} to explicitly use complete cases."
      ),
      class = "dar_error_model_missing_values"
    )
  }
  metadata <- metadata[complete, , drop = FALSE]

  targets <- get_var(rec)$vars
  for (target in targets) {
    if (is.numeric(metadata[[target]]) || is.integer(metadata[[target]])) {
      cli::cli_abort(
        "Target {.var {target}} must be categorical for automatic comparisons.",
        class = "dar_error_invalid_model"
      )
    }
    levs <- ordered_model_levels(metadata[[target]])
    if (length(levs) < 2L) {
      cli::cli_abort(
        "Target {.var {target}} must retain at least two levels after applying the model cohort.",
        class = "dar_error_non_estimable_model"
      )
    }
    metadata[[target]] <- factor(as.character(metadata[[target]]), levels = levs)
  }

  if (!is.null(model$time)) {
    time_levels <- ordered_model_levels(metadata[[model$time]], numeric_order = TRUE)
    if (length(time_levels) > model$max_time_levels) {
      cli::cli_abort(
        c(
          "x" = "Time variable {.var {model$time}} has {length(time_levels)} levels; the configured maximum is {model$max_time_levels}.",
          "i" = "Increase {.arg max_time_levels} explicitly if all pairwise contrasts are intended."
        ),
        class = "dar_error_too_many_time_levels"
      )
    }
    if (length(time_levels) < 2L) {
      cli::cli_abort(
        "Time variable {.var {model$time}} must retain at least two time points.",
        class = "dar_error_non_estimable_model"
      )
    }
    metadata[[model$time]] <- factor(
      as.character(metadata[[model$time]]), levels = time_levels, ordered = FALSE
    )
  }

  categorical <- names(metadata)[vapply(metadata, function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1))]
  for (var in intersect(model_vars, categorical)) {
    if (!is.factor(metadata[[var]])) {
      metadata[[var]] <- factor(metadata[[var]], levels = ordered_model_levels(metadata[[var]]))
    }
  }

  if (!is.null(model$reference)) {
    for (var in names(model$reference)) {
      metadata[[var]] <- stats::relevel(
        factor(metadata[[var]]), ref = model$reference[[var]]
      )
    }
  }

  fixed_formula <- reformulas::nobars(model$formula)
  design <- stats::model.matrix(fixed_formula, data = metadata)
  if (check_design && qr(design)$rank < ncol(design)) {
    cli::cli_abort(
      c(
        "x" = "The fixed-effect design matrix is not full rank.",
        "i" = "Remove aliased variables or empty combinations before preparing the recipe."
      ),
      class = "dar_error_non_estimable_model"
    )
  }

  interaction_targets <- model_interaction_targets(model$formula, targets, model$time)
  if (check_design && length(interaction_targets) > 0L) {
    for (target in interaction_targets) {
      cells <- table(metadata[[target]], metadata[[model$time]], useNA = "no")
      if (any(cells == 0L)) {
        cli::cli_abort(
          c(
            "x" = "The interaction {.code {target} * {model$time}} contains empty cells.",
            "i" = "Every requested condition-time combination must contain samples."
          ),
          class = "dar_error_non_estimable_model"
        )
      }
    }
  }

  random_terms <- reformulas::findbars(model$formula)
  if (check_design && length(random_terms) > 0L) {
    for (term in random_terms) {
      grouping <- all.vars(term[[3]])[[1]]
      group_sizes <- table(metadata[[grouping]])
      if (length(group_sizes) < 2L || any(group_sizes < 2L)) {
        cli::cli_abort(
          c(
            "x" = "Random effect grouped by {.var {grouping}} is not estimable.",
            "i" = "At least two groups with repeated observations are required."
          ),
          class = "dar_error_non_estimable_model"
        )
      }
      slopes <- setdiff(all.vars(term[[2]]), "1")
      for (slope in slopes) {
        distinct_per_group <- tapply(
          metadata[[slope]], metadata[[grouping]],
          function(x) length(unique(stats::na.omit(x)))
        )
        if (any(distinct_per_group < 2L)) {
          cli::cli_abort(
            c(
              "x" = "Random slope {.var {slope}} is not estimable within every {.var {grouping}} group.",
              "i" = "Each group must contain at least two distinct slope values."
            ),
            class = "dar_error_non_estimable_model"
          )
        }
      }
    }
  }

  fixed_vars <- all.vars(fixed_formula)
  confounders <- setdiff(fixed_vars, c(targets, model$time))
  resolved <- list(
    formula = model$formula,
    fixed_formula = fixed_formula,
    random_terms = random_terms,
    target_vars = targets,
    time = model$time,
    confounders = confounders,
    data = tibble::as_tibble(metadata),
    dropped_samples = dropped,
    matrix = design,
    interaction_targets = interaction_targets,
    contrast_mode = model$contrasts
  )
  resolved$contrast_plan <- build_contrast_plan(resolved)
  resolved
}

#' @noRd
model_interaction_targets <- function(formula, targets, time) {
  if (is.null(time)) {
    return(character())
  }
  labels <- attr(stats::terms(reformulas::nobars(formula)), "term.labels")
  labels <- labels[grepl(":", labels, fixed = TRUE)]
  unique(unlist(lapply(labels, function(label) {
    parts <- strsplit(label, ":", fixed = TRUE)[[1]]
    intersect(setdiff(parts, time), targets)
  }), use.names = FALSE))
}

#' @noRd
model_reference_row <- function(resolved) {
  data <- resolved$data
  row <- data[1, , drop = FALSE]
  for (var in all.vars(resolved$fixed_formula)) {
    value <- data[[var]]
    if (is.factor(value)) {
      row[[var]] <- factor(levels(value)[1], levels = levels(value))
    } else if (is.numeric(value) || is.integer(value)) {
      row[[var]] <- mean(value, na.rm = TRUE)
    } else {
      row[[var]] <- value[[1]]
    }
  }
  row
}

#' @noRd
model_design_row <- function(resolved, values = list()) {
  row <- model_reference_row(resolved)
  for (var in names(values)) {
    if (is.factor(resolved$data[[var]])) {
      row[[var]] <- factor(values[[var]], levels = levels(resolved$data[[var]]))
    } else {
      row[[var]] <- values[[var]]
    }
  }
  out <- stats::model.matrix(resolved$fixed_formula, data = row)
  out[1, colnames(resolved$matrix), drop = TRUE]
}

#' @noRd
model_pairs <- function(levels) {
  if (length(levels) < 2L) {
    return(list())
  }
  utils::combn(levels, 2L, simplify = FALSE)
}

#' @noRd
build_contrast_plan <- function(resolved) {
  rows <- list()
  add_row <- function(id, type, target, numerator, denominator,
                      at_var = NA_character_, at_level = NA_character_, weights) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      contrast_id = id,
      comparison = id,
      contrast_type = type,
      var = target,
      numerator = as.character(numerator),
      denominator = as.character(denominator),
      at_var = at_var,
      at_level = at_level,
      weights = list(as.numeric(weights))
    )
  }

  for (target in resolved$target_vars) {
    target_levels <- levels(resolved$data[[target]])
    target_pairs <- model_pairs(target_levels)
    has_time_interaction <- target %in% resolved$interaction_targets

    if (!has_time_interaction) {
      for (pair in target_pairs) {
        denominator <- pair[[1]]
        numerator <- pair[[2]]
        weights <- model_design_row(resolved, stats::setNames(list(numerator), target)) -
          model_design_row(resolved, stats::setNames(list(denominator), target))
        add_row(
          glue::glue("{target}[{numerator}-{denominator}]"),
          "main", target, numerator, denominator, weights = weights
        )
      }
      next
    }

    time_levels <- levels(resolved$data[[resolved$time]])
    time_pairs <- model_pairs(time_levels)
    for (time_level in time_levels) {
      for (pair in target_pairs) {
        denominator <- pair[[1]]
        numerator <- pair[[2]]
        num_values <- stats::setNames(list(numerator, time_level), c(target, resolved$time))
        den_values <- stats::setNames(list(denominator, time_level), c(target, resolved$time))
        weights <- model_design_row(resolved, num_values) - model_design_row(resolved, den_values)
        add_row(
          glue::glue("{target}[{numerator}-{denominator}]@{resolved$time}[{time_level}]"),
          "condition_within_time", target, numerator, denominator,
          resolved$time, time_level, weights
        )
      }
    }

    for (target_level in target_levels) {
      for (pair in time_pairs) {
        denominator <- pair[[1]]
        numerator <- pair[[2]]
        num_values <- stats::setNames(list(target_level, numerator), c(target, resolved$time))
        den_values <- stats::setNames(list(target_level, denominator), c(target, resolved$time))
        weights <- model_design_row(resolved, num_values) - model_design_row(resolved, den_values)
        add_row(
          glue::glue("{resolved$time}[{numerator}-{denominator}]@{target}[{target_level}]"),
          "time_within_condition", resolved$time, numerator, denominator,
          target, target_level, weights
        )
      }
    }

    if (identical(resolved$contrast_mode, "all")) {
      for (target_pair in target_pairs) {
        for (time_pair in time_pairs) {
          target_den <- target_pair[[1]]
          target_num <- target_pair[[2]]
          time_den <- time_pair[[1]]
          time_num <- time_pair[[2]]
          values <- function(target_level, time_level) {
            stats::setNames(list(target_level, time_level), c(target, resolved$time))
          }
          weights <-
            model_design_row(resolved, values(target_num, time_num)) -
            model_design_row(resolved, values(target_den, time_num)) -
            model_design_row(resolved, values(target_num, time_den)) +
            model_design_row(resolved, values(target_den, time_den))
          add_row(
            glue::glue("{target}[{target_num}-{target_den}]:{resolved$time}[{time_num}-{time_den}]"),
            "difference_in_differences", target,
            paste(target_num, time_num, sep = ":"),
            paste(target_den, time_den, sep = ":"),
            weights = weights
          )
        }
      }
    }
  }

  if (length(rows) == 0L) {
    return(tibble::tibble(
      contrast_id = character(), comparison = character(),
      contrast_type = character(), var = character(), numerator = character(),
      denominator = character(), at_var = character(), at_level = character(),
      weights = list()
    ))
  }
  dplyr::bind_rows(rows)
}

# ENGINE CAPABILITIES ---------------------------------------------------------

#' @noRd
model_engine_capabilities <- list(
  aldex = list(confounders = TRUE, time_interaction = TRUE, random = FALSE),
  ancom = list(confounders = TRUE, time_interaction = TRUE, random = TRUE),
  corncob = list(confounders = TRUE, time_interaction = TRUE, random = FALSE),
  deseq = list(confounders = TRUE, time_interaction = TRUE, random = FALSE),
  lefse = list(confounders = FALSE, time_interaction = FALSE, random = FALSE),
  maaslin = list(confounders = TRUE, time_interaction = TRUE, random = TRUE),
  wilcox = list(confounders = FALSE, time_interaction = FALSE, random = FALSE)
)

#' @noRd
step_engine <- function(step) {
  stringr::str_remove(class(step)[[1]], "^step_")
}

#' @noRd
is_da_step <- function(step) {
  !stringr::str_detect(step[["id"]], "subset|filter|rarefaction")
}

#' @noRd
model_step_status <- function(rec, step, resolved = NULL) {
  if (is.null(get_model(rec)) || !is_da_step(step)) {
    return(list(compatible = TRUE, reason = NA_character_))
  }
  if (is.null(resolved)) {
    resolved <- resolve_model(rec, check_design = FALSE)
  }
  engine <- step_engine(step)
  caps <- model_engine_capabilities[[engine]]
  if (is.null(caps)) {
    return(list(compatible = TRUE, reason = NA_character_))
  }
  reasons <- character()
  if (length(resolved$confounders) > 0L && !caps$confounders) {
    reasons <- c(reasons, "does not support adjusted confounders")
  }
  if (length(resolved$interaction_targets) > 0L && !caps$time_interaction) {
    reasons <- c(reasons, "does not support condition-time interactions")
  }
  if (length(resolved$random_terms) > 0L && !caps$random) {
    reasons <- c(reasons, "does not support random effects")
  }
  if (engine %in% c("wilcox", "lefse") && length(resolved$target_vars) != 1L) {
    reasons <- c(reasons, "requires exactly one simple target")
  }
  if (identical(engine, "deseq") && identical(step[["test"]], "LRT")) {
    reasons <- c(reasons, "LRT does not provide the requested pairwise contrasts")
  }
  if (identical(engine, "corncob") &&
      (!identical(step[["test"]], "Wald") || isTRUE(step[["boot"]]))) {
    reasons <- c(reasons, "central corncob contrasts require a non-bootstrap Wald test")
  }
  list(
    compatible = length(reasons) == 0L,
    reason = if (length(reasons) == 0L) NA_character_ else paste(unique(reasons), collapse = "; ")
  )
}

#' @noRd
model_steps_status <- function(rec, resolved = NULL) {
  steps <- purrr::keep(rec@steps, is_da_step)
  if (length(steps) == 0L) {
    return(tibble::tibble(
      step_id = character(), engine = character(), compatible = logical(),
      reason = character()
    ))
  }
  purrr::map_dfr(steps, function(step) {
    status <- model_step_status(rec, step, resolved)
    tibble::tibble(
      step_id = step[["id"]], engine = step_engine(step),
      compatible = status$compatible, reason = status$reason
    )
  })
}

#' @noRd
model_override_arguments <- function(step) {
  args <- switch(
    step_engine(step),
    ancom = c(
      "fix_formula", "rand_formula", "group",
      c("global", "pairwise", "dunnet", "trend")[vapply(
        c("global", "pairwise", "dunnet", "trend"),
        function(arg) isTRUE(step[[arg]]), logical(1)
      )]
    ),
    corncob = "formula_null",
    maaslin = c("random_effects", "reference"),
    lefse = "subclassCol",
    character()
  )
  args[vapply(args, function(arg) !is.null(step[[arg]]), logical(1))]
}

#' @noRd
report_model_step_warnings <- function(rec, step = NULL) {
  if (is.null(get_model(rec))) {
    return(invisible(NULL))
  }
  steps <- if (is.null(step)) purrr::keep(rec@steps, is_da_step) else list(step)
  if (length(steps) == 0L) {
    return(invisible(NULL))
  }
  resolved <- resolve_model(rec, check_design = FALSE)
  overrides <- purrr::map(steps, model_override_arguments)
  has_overrides <- lengths(overrides) > 0L
  if (any(has_overrides)) {
    details <- purrr::map2_chr(
      steps[has_overrides], overrides[has_overrides],
      ~ paste0(.x[["id"]], ": ", paste(.y, collapse = ", "))
    )
    cli::cli_warn(
      c(
        "!" = "The centralized model overrides design arguments in {sum(has_overrides)} step{?s}.",
        stats::setNames(details, rep("i", length(details))),
        "i" = "Method-specific thresholds and preprocessing controls are unchanged."
      ),
      class = "dar_warning_model_override"
    )
  }
  for (item in steps) {
    status <- model_step_status(rec, item, resolved)
    if (!status$compatible) {
      cli::cli_warn(
        c(
          "!" = "Step {.val {item[['id']]}} is incompatible with the centralized model and will be skipped by {.fun prep}.",
          "i" = status$reason
        ),
        class = "dar_warning_model_step_skipped"
      )
    }
  }
  invisible(NULL)
}

# MODEL COMPILATION -----------------------------------------------------------

#' Compile the centralized model for a DA engine
#' @noRd
compile_model <- function(rec, engine, resolved = NULL) {
  if (is.null(resolved)) {
    resolved <- resolve_model(rec)
  }
  if (is.null(resolved)) {
    return(NULL)
  }
  random_text <- if (length(resolved$random_terms) == 0L) {
    NULL
  } else {
    paste(vapply(resolved$random_terms, function(x) {
      paste0("(", paste(deparse(x), collapse = " "), ")")
    }, character(1)), collapse = " + ")
  }
  fixed_text <- paste(deparse(resolved$fixed_formula[[2]]), collapse = " ")
  full_text <- paste(deparse(resolved$formula), collapse = " ")
  list(
    engine = engine,
    formula = resolved$fixed_formula,
    formula_text = full_text,
    fix_formula = fixed_text,
    rand_formula = random_text,
    matrix = resolved$matrix,
    data = resolved$data,
    contrasts = resolved$contrast_plan,
    resolved = resolved
  )
}

#' Apply common complete-case and factor/reference handling
#' @noRd
apply_model_to_recipe <- function(rec) {
  resolved <- resolve_model(rec)
  if (is.null(resolved)) {
    return(list(rec = rec, resolved = NULL))
  }
  keep <- resolved$data$sample_id
  rec@phyloseq <- phyloseq::prune_samples(keep, rec@phyloseq)
  metadata <- resolved$data |>
    tibble::column_to_rownames("sample_id") |>
    data.frame(check.names = FALSE) |>
    phyloseq::sample_data()
  phyloseq::sample_data(rec@phyloseq) <- metadata
  list(rec = rec, resolved = resolved)
}

#' Add common model columns to a method result
#' @noRd
harmonize_model_result <- function(result, rec) {
  if (is.null(get_model(rec)) || !inherits(result, "data.frame")) {
    return(result)
  }
  result <- tibble::as_tibble(result)
  if (!"contrast_id" %in% names(result) && "comparison" %in% names(result)) {
    result$contrast_id <- result$comparison
  }
  if (!"comparison" %in% names(result) && "contrast_id" %in% names(result)) {
    result$comparison <- result$contrast_id
  }
  plan <- resolve_model(rec)$contrast_plan |>
    dplyr::select(
      "contrast_id", "contrast_type", planned_var = "var", "numerator",
      "denominator", "at_var", "at_level"
    )
  result <- dplyr::left_join(result, plan, by = "contrast_id")
  if (!"var" %in% names(result)) {
    result$var <- result$planned_var
  }
  dplyr::select(result, -dplyr::any_of("planned_var"))
}

#' Validate that an engine returned every planned contrast
#' @noRd
validate_model_result_contract <- function(result, rec, step_id) {
  if (is.null(get_model(rec))) {
    return(invisible(TRUE))
  }
  result <- harmonize_model_result(result, rec)
  expected <- resolve_model(rec)$contrast_plan$contrast_id
  observed <- unique(result$contrast_id)
  missing <- setdiff(expected, observed)
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Step {.val {step_id}} did not return every planned contrast.",
        "i" = "Missing contrast{?s}: {.val {missing}}."
      ),
      class = "dar_error_invalid_result_contract"
    )
  }
  invisible(TRUE)
}

#' @noRd
flatten_model_output <- function(result) {
  if (inherits(result, "data.frame")) {
    return(tibble::as_tibble(result))
  }
  if (is.list(result)) {
    tables <- purrr::map(result, flatten_model_output)
    tables <- purrr::keep(tables, ~ inherits(.x, "data.frame"))
    if (length(tables) > 0L) {
      return(dplyr::bind_rows(tables))
    }
  }
  tibble::tibble()
}

#' @noRd
harmonize_model_output <- function(result, rec) {
  if (inherits(result, "data.frame")) {
    return(harmonize_model_result(result, rec))
  }
  if (is.list(result)) {
    return(purrr::map(result, harmonize_model_output, rec = rec))
  }
  result
}
