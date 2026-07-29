# Internal validation helpers -------------------------------------------------

#' Abort when a Recipe violates its S4 invariants
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#' @param arg_name The name of the argument in the calling function.
#' @param context Optional context describing where validation failed.
#'
#' @return Invisible `NULL` if valid, otherwise throws an error.
#' @noRd
#' @keywords internal
validate_recipe_object <- function(rec, arg_name = "rec", context = NULL) {
  problems <- methods::validObject(rec, test = TRUE, complete = TRUE)

  if (isTRUE(problems)) {
    return(invisible(NULL))
  }

  bullets <- stats::setNames(as.character(problems), rep("!", length(problems)))
  messages <- c(
    "x" = "The argument {.arg {arg_name}} is an invalid {.cls Recipe} object.",
    bullets
  )
  if (!is.null(context)) {
    messages <- c(messages, "i" = context)
  }

  cli::cli_abort(
    messages,
    class = "dar_error_invalid_recipe"
  )
}

#' Check if object is a valid Recipe (and NOT Prepped)
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#'
#' @return Invisible NULL if valid, throws error otherwise.
#' @noRd
#' @keywords internal
check_recipe <- function(rec, arg_name = "rec") {
  
  if (!checkmate::test_class(rec, "Recipe")) {
    cli::cli_abort(
      c(
        "x" = "The argument {.arg {arg_name}} must be a {.cls Recipe} object.",
        "!" = "You supplied an object of class {.cls {class(rec)[1]}} instead.",
        "i" = "Did you forget to initialize the recipe with {.fun recipe}?"
      ),
      class = "dar_error_invalid_input"
    )
  }
  
  if (is(rec, "PrepRecipe")) {
    cli::cli_abort(
      c(
        "x" = "The argument {.arg {arg_name}} must be an unprepped {.cls Recipe}.",
        "!" = "You passed a {.cls Recipe} that has already been prepared ({.fun prep} was called).",
        "i" = "Please pass the {.cls Recipe} object before the {.fun prep} step."
      ),
      class = "dar_error_invalid_input"
    )
  }

  validate_recipe_object(rec, arg_name = arg_name)
  
  invisible(NULL)
}

#' Check if object is a Prepared Recipe (PrepRecipe)
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#' @noRd
#' @keywords internal
check_prep_recipe <- function(rec, arg_name = "rec") {
  
  if (!checkmate::test_class(rec, "PrepRecipe")) {
    if (checkmate::test_class(rec, "Recipe")) {
      cli::cli_abort(
        c(
          "x" = "The argument {.arg {arg_name}} must be a prepared {.cls PrepRecipe}.",
          "!" = "You supplied an unprepped {.cls Recipe}.",
          "i" = "Please run {.fun prep} on your recipe before using this function."
        ),
        class = "dar_error_invalid_input"
      )
    }
    
    cli::cli_abort(
      c(
        "x" = "The argument {.arg {arg_name}} must be a {.cls PrepRecipe} object.",
        "!" = "You supplied an object of class {.cls {class(rec)[1]}} instead.",
        "i" = "This function expects the output of {.code prep(recipe)}."
      ),
      class = "dar_error_invalid_input"
    )
  }

  validate_recipe_object(rec, arg_name = arg_name)
  
  invisible(NULL)
}

#' Check if object is EITHER a Recipe OR a PrepRecipe
#'
#' @param rec The object to check.
#' @param arg_name The name of the argument in the calling function.
#'
#' @return Invisible NULL if valid, throws error otherwise.
#' @noRd
#' @keywords internal
check_any_recipe <- function(rec, arg_name = "rec") {
  
  # Since PrepRecipe inherits from Recipe, we only need to check the base class.
  if (!checkmate::test_class(rec, "Recipe")) {
    cli::cli_abort(
      c(
        "x" = "The argument {.arg {arg_name}} must be a {.cls Recipe} or {.cls PrepRecipe} object.",
        "!" = "You supplied an object of class {.cls {class(rec)[1]}} instead.",
        "i" = "Please supply a valid recipe object created with {.fun recipe}."
      ),
      class = "dar_error_invalid_input"
    )
  }

  validate_recipe_object(rec, arg_name = arg_name)
  
  invisible(NULL)
}
