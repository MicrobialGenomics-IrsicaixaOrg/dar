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
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species") |> step_maaslin() 
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
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
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
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
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
  target_list <- if (is_prep) rec@bakes else rec@steps
  obj_cmp <- object
  obj_cmp[["id"]] <- NULL
  
  dupl_rec <- any(purrr::map_lgl(target_list, ~ {
    step_cmp <- .x
    step_cmp[["id"]] <- NULL
    identical(step_cmp, obj_cmp)
  }))
  
  if (dupl_rec) {
    expr <- stringr::str_replace(step_to_expr(object), 'run', 'step')
    rlang::inform(
      c(
        "!" = "This step is already defined with the same parameters and will be skipped.",
        "i" = expr
      ),
      use_cli_format = TRUE
    )
    return(rec)
  } 
  
  if (is_prep) {
    rec@bakes[[length(rec@bakes) + 1]] <- object
  } else {
    rec@steps[[length(rec@steps) + 1]] <- object
  }
  
  rec
}