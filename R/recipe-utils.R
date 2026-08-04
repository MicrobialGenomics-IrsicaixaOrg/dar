# CANONICAL ANALYSIS CONFIGURATION --------------------------------------------

#' @noRd
recipe_targets <- function(rec) {
  model <- normalize_model_spec(rec, rec@model)
  if (!is.null(model)) {
    return(model$targets)
  }
  if (is.null(rec@var_info) || !"vars" %in% names(rec@var_info)) {
    return(character())
  }
  rec@var_info$vars
}

#' @noRd
recipe_tax_level <- function(rec) {
  model <- normalize_model_spec(rec, rec@model)
  if (!is.null(model)) {
    return(model$tax_level)
  }
  if (is.null(rec@tax_info) || !"tax_lev" %in% names(rec@tax_info)) {
    return(character())
  }
  rec@tax_info$tax_lev
}

#' @noRd
warn_deprecated_selector <- function(selector) {
  cli::cli_warn(
    c(
      "!" = "{.fun {selector}} is deprecated and will become defunct in a future Bioconductor release.",
      "i" = "Define analysis targets and taxonomic resolution with {.fun add_model}."
    ),
    class = "dar_warning_deprecated_selector"
  )
}

# GETTERS ----------------------------------------------------------------------

#' Return legacy target information from a recipe
#'
#' `get_var()` is deprecated. Use `get_model(rec)$targets` for modeled recipes.
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_var
#' @return Tibble containing `var_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' suppressWarnings(get_var(rec))
get_var <- function(rec) {
  check_any_recipe(rec)
  warn_deprecated_selector("get_var")
  tibble::tibble(vars = recipe_targets(rec))
}

#' Return legacy taxonomic information from a recipe
#'
#' `get_tax()` is deprecated. Use `get_model(rec)$tax_level` for modeled
#' recipes.
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_tax
#' @return Tibble containing `tax_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' suppressWarnings(get_tax(rec))
get_tax <- function(rec) {
  check_any_recipe(rec)
  warn_deprecated_selector("get_tax")
  tibble::tibble(tax_lev = recipe_tax_level(rec))
}

#' Returns phyloseq from Recipe-class object
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_phy
#' @return Phyloseq class object
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' get_phy(rec)
get_phy <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq
}

# SETTERS ----------------------------------------------------------------------

#' Add legacy variables of interest to a recipe
#'
#' `add_var()` is deprecated. Supply `targets` to [add_model()] instead.
#' @param rec A `Recipe` object (must not be prepped).
#' @param var_info A character string of column names corresponding to variables
#'   that will be used in any context.
#'
#' @aliases add_var
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' rec <- suppressWarnings(add_var(rec, var_info = "RiskGroup2"))
#' rec
#' @tests
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#'
#' expect_warning(
#'   expect_s4_class(add_var(rec, "RiskGroup2"), "Recipe"),
#'   class = "dar_warning_deprecated_selector"
#' )
#' expect_error(
#'   suppressWarnings(add_var(rec, "missing_variable")),
#'   class = "dar_error_invalid_recipe"
#' )
add_var <- function(rec, var_info) {
  check_recipe(rec)
  warn_deprecated_selector("add_var")
  model <- normalize_model_spec(rec, rec@model)
  if (is.null(model)) {
    rec@var_info <- tibble::tibble(vars = var_info)
  } else {
    rec <- add_model(
      rec, model$formula, targets = var_info, tax_level = model$tax_level,
      time = model$time, reference = model$reference,
      contrasts = model$contrasts, na_action = model$na_action,
      max_time_levels = model$max_time_levels
    )
  }
  validate_recipe_object(rec)
  rec
}

#' Add a legacy taxonomic level of interest to a recipe
#'
#' `add_tax()` is deprecated. Supply `tax_level` to [add_model()] instead.
#' @param rec A `Recipe` object (must not be prepped).
#' @param tax_info A character string of taxonomic levels that will be used in
#'   any context.
#'
#' @aliases add_tax
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' rec <- suppressWarnings(add_tax(rec, tax_info = "Species"))
#' rec
#' @tests
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#'
#' expect_warning(
#'   expect_s4_class(add_tax(rec, "Species"), "Recipe"),
#'   class = "dar_warning_deprecated_selector"
#' )
#' expect_error(
#'   suppressWarnings(add_tax(rec, "Missing_rank")),
#'   class = "dar_error_invalid_recipe"
#' )
add_tax <- function(rec, tax_info) {
  check_recipe(rec)
  warn_deprecated_selector("add_tax")
  tax_info <- stringr::str_to_sentence(tax_info)
  model <- normalize_model_spec(rec, rec@model)
  if (is.null(model)) {
    rec@tax_info <- tibble::tibble(tax_lev = tax_info)
  } else {
    rec <- add_model(
      rec, model$formula, targets = model$targets, tax_level = tax_info,
      time = model$time, reference = model$reference,
      contrasts = model$contrasts, na_action = model$na_action,
      max_time_levels = model$max_time_levels
    )
  }
  validate_recipe_object(rec)
  rec
}

# PHYLOSEQ EXTRACTION ----------------------------------------------------------

#' Extracts tax_table from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' tax_table(rec)
tax_table <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq@tax_table %>%
    to_tibble("taxa_id") %>%
    dplyr::select(taxa_id, taxa = !!recipe_tax_level(rec))
}

#' Extracts sample_data from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' sample_data(rec)
sample_data <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq %>%
    phyloseq::sample_data() %>%
    to_tibble("sample_id") %>%
    dplyr::select(sample_id, dplyr::all_of(recipe_targets(rec)))
}

#' Extracts otu_table from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' otu_table(rec)
otu_table <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq %>%
    phyloseq::otu_table() %>%
    to_tibble("taxa_id")
}


# EXECUTION (PREP) -------------------------------------------------------------

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_prep <- function(x, ...) { c("furrr", "future") }

#' Performs all the steps defined in a Recipe
#'
#' For a Recipe with at least one preprocessing or DA operation run the steps in
#' a convenient order.
#'
#' @param rec A `Recipe` object.
#' @param parallel if FALSE, no palatalization. if TRUE, parallel execution
#'   using future and furrr packages.
#' @param workers Number of workers for palatalization.
#' @param force Force the reexecution of all steps. This remove previous
#'   results.
#'
#' @aliases prep
#' @return A `PrepRecipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <-
#'   recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Class") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = function(x) sum(x > 0) >= (0.03 * length(x))) |>
#'   step_maaslin()
#'
#' da_results <- prep(rec)
#'
#' @tests
#' data(metaHIV_phy)
#' invalid_rec <- suppressWarnings(recipe(
#'   metaHIV_phy,
#'   var_info = "RiskGroup2",
#'   tax_info = "Species"
#' ))
#' invalid_rec@var_info <- tibble::tibble(vars = "missing_variable")
#' expect_error(
#'   prep(invalid_rec, parallel = FALSE),
#'   class = "dar_error_invalid_recipe"
#' )
#'
#' empty_filter_rec <- suppressWarnings(recipe(
#'   metaHIV_phy,
#'   var_info = "RiskGroup2",
#'   tax_info = "Species"
#' )) |>
#'   step_filter_taxa(
#'     .f = function(x) FALSE,
#'     id = "filter_taxa__empty"
#'   )
#' expect_error(
#'   prep(empty_filter_rec, parallel = FALSE),
#'   regexp = "filter_taxa__empty",
#'   class = "dar_error_preprocessing_step"
#' )
prep <- function(rec,
                 parallel = TRUE,
                 workers = future::availableCores(constraints = "connections-16"),
                 force = FALSE) {

  if (inherits(rec, "PrepRecipe") && !force) {
    cli::cli_abort(
      c(
        "x" = "The input recipe has already been prepped!",
        "i" = "To force the rerun of all steps, please run {.code prep(rec, force = TRUE)}."
      ),
      class = "dar_error_recipe_already_prepped"
    )
  }

  check_any_recipe(rec)
  validate_engine_args_model(rec)

  initial_da_steps <- purrr::keep(rec@steps, is_da_step)
  warn_model_free_da(rec)
  static_status <- model_steps_status(rec)
  dependency_steps <- model_dependency_steps(rec, static_status)

  check <- utils::capture.output(required_deps(rec, dependency_steps))
  if (length(check) > 0) {
    cli::cli_abort(
      c(
        "x" = "Not all necessary dependencies are installed.",
        "i" = "Use {.code required_deps(rec)} to see how to install them."
      ),
      class = "dar_error_missing_dependencies"
    )
  }

  ## Phyloseq preprocessing steps, in their configured order
  preprocessing_steps <- purrr::keep(rec@steps, is_preprocessing_step)
  to_execute <- purrr::map(preprocessing_steps, step_to_call)

  for (.i in seq_along(to_execute)) {
    step_id <- preprocessing_steps[[.i]][["id"]]
    rec <- tryCatch(
      base::eval(to_execute[[.i]]),
      error = function(cnd) {
        cli::cli_abort(
          c(
            "x" = "Preprocessing step {.val {step_id}} failed.",
            "i" = "The preprocessing operation could not be completed."
          ),
          class = "dar_error_preprocessing_step",
          parent = cnd
        )
      }
    )
    validate_recipe_object(
      rec,
      context = glue::glue("Validation failed after preprocessing step `{step_id}`.")
    )
  }

  resolved <- NULL
  if (!is.null(get_model(rec))) {
    modeled <- apply_model_to_recipe(rec)
    rec <- modeled$rec
    resolved <- modeled$resolved
    validate_recipe_object(
      rec,
      context = "Validation failed after applying the centralized model cohort."
    )
  }

  ## DA steps
  final_status <- model_steps_status(rec, resolved)
  skipped <- dplyr::filter(final_status, !.data$compatible)
  executed_ids <- final_status$step_id[final_status$compatible]

  if (nrow(skipped) > 0L) {
    cli::cli_warn(
      c(
        "!" = "{nrow(skipped)} differential-abundance step{?s} {?was/were} skipped because {?it is/they are} incompatible with the centralized model.",
        stats::setNames(
          paste0(skipped$step_id, ": ", skipped$reason),
          rep("i", nrow(skipped))
        )
      ),
      class = "dar_warning_model_step_skipped"
    )
  }

  if (length(initial_da_steps) > 0L && length(executed_ids) == 0L) {
    cli::cli_abort(
      c(
        "x" = "No differential-abundance step can execute the centralized model.",
        "i" = "Simplify the model or add a compatible method."
      ),
      class = "dar_error_no_compatible_steps"
    )
  }

  da_steps <- rec@steps %>%
    purrr::keep(~ is_da_step(.x) && .x[["id"]] %in% executed_ids) %>%
    rarefy_msg()

  da_names <- purrr::map_chr(da_steps, ~ .x[["id"]])
  da_calls <- da_steps %>% purrr::map(step_to_call)

  if (parallel) {
    recipes_pkg_check(required_pkgs_prep(), "prep()")
    with(future::plan(future::multisession, workers = workers), local = TRUE)

    res <-
      da_calls %>%
      furrr::future_map(~ {
          rec <- rec
          base::eval(.x)
        },
        .options = furrr::furrr_options(
          seed = TRUE,
          packages = "dar",
          conditions = structure("condition", exclude = "message")
        )
      )

    names(res) <- da_names

  } else {
    res <-
      da_calls %>%
      purrr::map(~ {
        rec <- rec
        base::eval(.x)
      })

    names(res) <- da_names
  }

  for (step_id in names(res)) {
    if (!is.null(get_model(rec))) {
      res[[step_id]] <- harmonize_model_output(res[[step_id]], rec)
      validate_model_result_contract(flatten_model_output(res[[step_id]]), rec, step_id)
    }
    standardize_da_result(res[[step_id]], rec, step_id)
  }

  execution <- list(
    contrasts = if (is.null(resolved)) NULL else resolved$contrast_plan,
    dropped_samples = if (is.null(resolved)) character() else resolved$dropped_samples,
    executed_steps = names(res),
    skipped_steps = skipped
  )

  prep_recipe(rec, res, list(), execution = execution)
}


# RESULTS EXTRACTION -----------------------------------------------------------

#' Returns data.frame with OTU intersection between methods
#'
#' @param rec A `PrepRecipe` object.
#' @param steps character vector with step_ids to take in account.
#' @param tidy Boolean indicating if result must be in tidy format.
#'
#' @aliases intersection_df
#' @return data.frame class object
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' df <- intersection_df(test_prep_rec)
#' head(df)
intersection_df <- function(rec, steps = steps_ids(rec, "da"), tidy = FALSE) {

  check_prep_recipe(rec)

  if (!is.null(get_model(rec))) {
    keys <- c("taxa_id", "contrast_id", "effect")
    df <- .otu_effect_direction(rec, steps = steps) %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(c(keys, "method")))) %>%
      dplyr::mutate(value = 1L) %>%
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(keys), names_from = "method",
        values_from = "value", values_fill = 0L
      )
    for (step_id in setdiff(steps, names(df))) {
      df[[step_id]] <- 0L
    }
    df <- dplyr::select(df, dplyr::all_of(c(keys, steps)))
    if (tidy) {
      df <- tidyr::pivot_longer(
        df, dplyr::all_of(steps), names_to = "name", values_to = "value"
      )
    }
    return(as.data.frame(df))
  }

  tidy_result <- tidy_results(rec, steps = steps, significant_only = TRUE)
  df <-
    steps %>%
    purrr::set_names() %>%
    purrr::map_dfc( ~ {
      taxa <- tidy_result %>%
        dplyr::filter(.data$step_id == .x) %>%
        dplyr::pull(taxa_id)

      rownames(rec@phyloseq@otu_table) %>%
        tibble::tibble(taxa_id = .) %>%
        dplyr::mutate(!!.x := dplyr::if_else(taxa_id %in% taxa, 1, 0)) %>%
        dplyr::select(!!.x)
    }) %>%
    dplyr::mutate(taxa_id = rownames(rec@phyloseq@otu_table), .before = 1) %>%
    as.data.frame()

  if (tidy) { df <- tidyr::pivot_longer(df, -taxa_id) }
  df
}

#' Overlap of significant OTUs between tested methods.
#'
#' @param rec A `PrepRecipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param type Indicates whether to use all taxa ("all") or only those that are
#'   differentially abundant in at least one method ("da"). Default as "all".
#'
#' @aliases overlap_df
#' @return df
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' df <- overlap_df(test_prep_rec, steps_ids(test_prep_rec, "da"))
#' head(df)
overlap_df <- function(rec, steps = steps_ids(rec, "da"), type = "all") {

  check_prep_recipe(rec)

  df <-
    intersection_df(rec, steps = steps) %>%
    tibble::as_tibble() %>%
    dplyr::select(dplyr::all_of(steps))

  if (type == "da") {
    df <-
      dplyr::rowwise(df) %>%
      dplyr::mutate(sum = sum(dplyr::across(dplyr::all_of(steps))), .before = 1) %>%
      dplyr::filter(sum != 0) %>%
      dplyr::select(dplyr::all_of(steps))
  }

  names(df) %>%
    purrr::map_dfr(~ {
      names(df) %>%
        purrr::map_dfc(function(.y) {
          res <- tibble::tibble(
            var_1 = dplyr::pull(df, .x),
            var_2 = dplyr::pull(df, .y),
            sum = var_1 - var_2
          )

          tibble::tibble(
            !!dplyr::sym(.y) := nrow(dplyr::filter(res, sum == 0)) / nrow(res)
          )
        })
    }) %>%
    data.frame(row.names = names(df))
}

#' Extract results from defined bake
#'
#' @param rec A `PrepRecipe` object.
#' @param bake Name or index of the bake to extract.
#'
#' @aliases cool
#' @return tbl_df
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' rec <- bake(test_prep_rec)
#' cool(rec)
cool <- function(rec, bake = 1) {

  check_prep_recipe(rec)

  all_bakes <- rec@bakes
  all_names <- all_bakes %>% purrr::map_chr(~ as.character(.x[["id"]]))

  if (is.numeric(bake) && length(all_bakes) < bake) {
    cli::cli_abort(
      c(
       "x" = "Bake index is not defined in the {.cls PrepRecipe}",
       "i" = "Run {.code bake(PrepRecipe)} and then try with {.code cool()}."
      ),
      class = "dar_error_missing_bake_index"
    )
  }

  if (!is.numeric(bake) && !bake %in% all_names) {
    cli::cli_abort(
      c(
        "x" = "Bake name is not defined in the {.cls PrepRecipe}!",
        "i" = "Run {.code bake(PrepRecipe)} and then try with {.code cool()}."
      ),
      class = "dar_error_missing_bake_name"
    )
  }

  if (!is.numeric(bake)) {
    bake <- which(all_names == bake)
  }

  to_execute <-
    rec@bakes %>%
    .[bake] %>%
    purrr::map_chr(step_to_expr)

  eval(parse(text = to_execute))
}
