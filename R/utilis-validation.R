#' @title Internal Validation Helpers
#' @description
#' Helper functions to validate arguments across the package.
#' Centralizes error messages and reduces boilerplate.
#'
#' @noRd
#' @keywords internal

#' Check if object is a valid Recipe (and NOT Prepped)
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#'
#' @return Invisible NULL if valid, throws error otherwise.
check_recipe <- function(rec, arg_name = "rec") {
  
  if (!checkmate::test_class(rec, "Recipe")) {
    rlang::abort(
      message = glue::glue("The argument `{arg_name}` must be a 'Recipe' object."),
      class = "dar_error_invalid_input",
      body = c(
        "x" = glue::glue("You supplied an object of class `{class(rec)[1]}` instead."),
        "i" = "Did you forget to initialize the recipe with `recipe()`?"
      )
    )
  }
  
  if (is(rec, "PrepRecipe")) {
    rlang::abort(
      message = glue::glue("The argument `{arg_name}` must be an unprepped Recipe."),
      class = "dar_error_invalid_input",
      body = c(
        "x" = "You passed a Recipe that has already been prepared (`prep()` was called).",
        "i" = "Please pass the Recipe object BEFORE the prep() step."
      )
    )
  }
  
  invisible(NULL)
}

#' Check if object is a Prepared Recipe (PrepRecipe)
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#'
#' @return Invisible NULL if valid, throws error otherwise.
check_prep_recipe <- function(rec, arg_name = "rec") {
  
  if (!checkmate::test_class(rec, "PrepRecipe")) {
    if (checkmate::test_class(rec, "Recipe")) {
      rlang::abort(
        message = glue::glue("The argument `{arg_name}` must be a prepared recipe ('PrepRecipe')."),
        class = "dar_error_invalid_input",
        body = c(
          "x" = "You supplied an unprepped 'Recipe'.",
          "i" = "Please run `prep()` on your recipe before using this function."
        )
      )
    }
    
    rlang::abort(
      message = glue::glue("The argument `{arg_name}` must be a 'PrepRecipe' object."),
      class = "dar_error_invalid_input",
      body = c(
        "x" = glue::glue("You supplied an object of class `{class(rec)[1]}` instead."),
        "i" = "This function expects the output of `prep(recipe)`."
      )
    )
  }
  
  invisible(NULL)
}

#' Check if object is EITHER a Recipe OR a PrepRecipe
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#'
#' @return Invisible NULL if valid, throws error otherwise.
check_any_recipe <- function(rec, arg_name = "rec") {
  
  # Since PrepRecipe inherits from Recipe, we only need to check the base class.
  if (!checkmate::test_class(rec, "Recipe")) {
    rlang::abort(
      message = glue::glue("The argument `{arg_name}` must be a 'Recipe' or 'PrepRecipe' object."),
      class = "dar_error_invalid_input",
      body = c(
        "x" = glue::glue("You supplied an object of class `{class(rec)[1]}` instead."),
        "i" = "Please supply a valid recipe object created with `recipe()`."
      )
    )
  }
  
  invisible(NULL)
}