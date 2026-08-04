# Canonical recipe-step registry ---------------------------------------------

#' Supported recipe steps and their execution roles
#'
#' The concrete S3 class is the only source of truth for step semantics. User
#' identifiers are deliberately absent from this registry and never influence
#' classification or execution.
#'
#' @noRd
register_step_methods <- function(methods, role) {
  registry <- lapply(methods, function(method) {
    list(
      method = method,
      role = role,
      runner = paste0("run_", method),
      engine = if (identical(role, "da")) method else NULL
    )
  })
  stats::setNames(registry, paste0("step_", methods))
}

#' @noRd
recipe_step_registry <- c(
  register_step_methods(
    c(
      "filter_by_abundance", "filter_by_prevalence", "filter_by_rarity",
      "filter_by_variance", "filter_taxa", "subset_taxa", "rarefaction"
    ),
    "preprocessing"
  ),
  register_step_methods(
    c(
      "aldex", "ancom", "corncob", "deseq", "lefse", "linda", "maaslin",
      "wilcox"
    ),
    "da"
  ),
  register_step_methods("bake", "bake")
)

#' @noRd
step_validation_problems <- function(object) {
  problems <- character()
  if (!inherits(object, "step")) {
    return("must inherit from the `step` class")
  }

  matched_classes <- intersect(class(object), names(recipe_step_registry))
  if (length(matched_classes) == 0L) {
    problems <- c(problems, "must have one supported concrete `step_*` class")
  } else if (length(matched_classes) > 1L) {
    problems <- c(
      problems,
      "must not inherit from multiple supported concrete `step_*` classes"
    )
  }

  id <- object[["id"]]
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    problems <- c(problems, "must contain one non-empty, non-missing character `id`")
  }
  problems
}

#' @noRd
step_slot_validation_problems <- function(configured_steps, slot_name, roles) {
  problems <- character()
  for (index in seq_along(configured_steps)) {
    configured_step <- configured_steps[[index]]
    current <- step_validation_problems(configured_step)
    if (length(current) > 0L) {
      problems <- c(
        problems,
        paste0("`", slot_name, "[[", index, "]]` ", current)
      )
      next
    }
    role <- step_role(configured_step)
    if (!role %in% roles) {
      problems <- c(
        problems,
        paste0(
          "`", slot_name, "[[", index, "]]` has role `", role,
          "`; expected one of: ", paste(roles, collapse = ", ")
        )
      )
    }
  }
  problems
}

#' @noRd
valid_step_ids <- function(configured_steps) {
  unlist(
    lapply(configured_steps, function(configured_step) {
      if (!is.list(configured_step)) {
        return(NULL)
      }
      id <- configured_step[["id"]]
      if (is.character(id) && length(id) == 1L && !is.na(id) && nzchar(id)) {
        id
      } else {
        NULL
      }
    }),
    use.names = FALSE
  )
}

#' @noRd
duplicate_step_id_problems <- function(configured_steps, context = NULL) {
  ids <- valid_step_ids(configured_steps)
  duplicate_ids <- unique(ids[duplicated(ids)])
  if (length(duplicate_ids) == 0L) {
    return(character())
  }
  scope <- if (is.null(context)) "" else paste0(" across ", context)
  paste0(
    "Recipe step IDs must be unique", scope, "; duplicated: ",
    paste(duplicate_ids, collapse = ", ")
  )
}

#' @noRd
step_descriptor <- function(object) {
  problems <- step_validation_problems(object)
  if (length(problems) > 0L) {
    cli::cli_abort(
      c(
        "x" = "The supplied recipe step is invalid.",
        stats::setNames(problems, rep("!", length(problems)))
      ),
      class = "dar_error_invalid_step"
    )
  }
  concrete_class <- intersect(class(object), names(recipe_step_registry))[[1L]]
  recipe_step_registry[[concrete_class]]
}

#' @noRd
step_role <- function(object) step_descriptor(object)$role

#' @noRd
step_method <- function(object) step_descriptor(object)$method

#' @noRd
step_runner <- function(object) {
  get(step_descriptor(object)$runner, envir = asNamespace("dar"))
}

#' @noRd
is_preprocessing_step <- function(object) {
  identical(step_role(object), "preprocessing")
}

#' @noRd
is_da_step <- function(object) identical(step_role(object), "da")

#' @noRd
is_bake_step <- function(object) identical(step_role(object), "bake")

#' @noRd
validate_step_role <- function(object, expected_roles) {
  descriptor <- step_descriptor(object)
  if (!descriptor$role %in% expected_roles) {
    cli::cli_abort(
      c(
        "x" = "Step {.cls {class(object)[[1L]]}} cannot be added here.",
        "i" = "Expected role: {.val {expected_roles}}; received {.val {descriptor$role}}."
      ),
      class = "dar_error_invalid_step"
    )
  }
  invisible(descriptor)
}

#' Overall Wrappers to Make New `step_X` or `check_Y` Objects
#'
#' `step` sets the class of the `step` and `check` is for checks.
#'
#' @param subclass A character string for the resulting class. For example, if
#'   `subclass = "blah"` the step object that is returned has class `step_blah`
#'   or `check_blah` depending on the context.
#' @param ... All arguments to the operator that should be returned.
#' @param .prefix Prefix to the subclass created.
#'
#' @keywords internal
#' @return An updated step or check with the new class.
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
#'   step_maaslin()
#'     
#' expect_s4_class(rec, "Recipe")
#' expect_equal(length(rec@steps), 1)
step <- function(subclass, ..., .prefix = "step_") {
  structure(list(...), class = c(paste0(.prefix, subclass), "step"))
}

#' @rdname step
#' @keywords internal
#' @autoglobal
check <- function(subclass, ..., .prefix = "check_") {
  structure(list(...), class = c(paste0(.prefix, subclass), "check"))
}

#' Add a New Operation to the Current Recipe
#'
#' `add_step` adds a step to the last location in the Recipe. `add_check` does
#' the same for checks.
#'
#' @param rec A [Recipe()] or [PrepRecipe()].
#' @param object A step or check object.
#' @return A updated [Recipe()] or [PrepRecipe()] with the new operation in the last slot.
#' @keywords internal
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' 
#' # Internally, step_ functions use add_step to append themselves
#' rec <- step_maaslin(rec)
#' length(rec@steps) # Returns 1
#' 
#' # If we try to add the exact same step, it will inform us and skip it
#' rec <- step_maaslin(rec)
#' length(rec@steps) # Still returns 1
#' 
#' @tests
#' data(metaHIV_phy)
#' 
#' # 1. Test adding steps to a normal Recipe
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' rec <- step_maaslin(rec, id = "test_1")
#' expect_equal(length(rec@steps), 1)
#' 
#' # 2. Test duplicate skipping mechanism in Recipe
#' # Adding the exact same step should not increase the length
#' expect_snapshot({
#'  rec <- step_maaslin(rec, id = "test_2") 
#' })
#' expect_equal(length(rec@steps), 1)
#' 
#' # 3. Test adding bakes to a PrepRecipe
#' data(test_prep_rec)
#' initial_bakes <- length(test_prep_rec@bakes)
#' 
#' # Bake adds a step to the @bakes slot
#' prepped <- bake(test_prep_rec)
#' expect_equal(length(prepped@bakes), initial_bakes + 1)
#' 
#' # Duplicate bake should be skipped
#' expect_snapshot({
#'    prepped_dup <- bake(prepped)
#' })
#' expect_equal(length(prepped_dup@bakes), initial_bakes + 1)
add_step <- function(rec, object) {
  
  check_any_recipe(rec)
  is_prep <- methods::is(rec, "PrepRecipe")
  expected_roles <- if (is_prep) "bake" else c("preprocessing", "da")
  validate_step_role(object, expected_roles)
  target_list <- if (is_prep) rec@bakes else rec@steps
  obj_cmp <- object
  obj_cmp[["id"]] <- NULL
  
  dupl_rec <- any(purrr::map_lgl(target_list, ~ {
    step_cmp <- .x
    step_cmp[["id"]] <- NULL
    identical(step_cmp, obj_cmp)
  }))
  
  if (dupl_rec) {
    expr <- step_call_label(object)
    cli::cli_inform(
      c(
        "!" = "This step is already defined with the same parameters and will be skipped.",
        "i" = "{.code {expr}}"
      )
    )
    return(rec)
  }

  configured_operations <- rec@steps
  if (is_prep) {
    configured_operations <- c(configured_operations, rec@bakes)
  }
  configured_ids <- purrr::map_chr(configured_operations, ~ .x[["id"]])
  if (object[["id"]] %in% configured_ids) {
    cli::cli_abort(
      c(
        "x" = "Step ID {.val {object[['id']]}} is already configured.",
        "i" = "Every recipe step must have a unique ID."
      ),
      class = "dar_error_invalid_step"
    )
  }
  
  if (is_prep) {
    rec@bakes[[length(rec@bakes) + 1]] <- object
  } else {
    rec@steps[[length(rec@steps) + 1]] <- object
    if (is_da_step(object)) {
      report_model_step_warnings(rec, object)
    }
  }

  validate_recipe_object(rec)
  rec
}
