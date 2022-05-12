#' Create a recipe for preprocessing data
#'
#' A recipe is a description of the steps to be applied to a data set in
#'   order to prepare it for data analysis.
#'
#' @aliases recipe recipe.default recipe.formula
#' @export
recipe <- function(x, ...) {
  UseMethod("recipe")
}

#' @rdname recipe
#' @export
recipe.default <- function(x, ...) {
  rlang::abort("`x` should be a phyloseq")
}

#' @export
recipe.phyloseq <- function(x, vars = NULL) {

  s_data <-
    phyloseq::sample_data(x) %>%
    tibble::as_tibble(rownames = "SampleID")


  if (is.null(vars)) {
    vars <- colnames(s_data)
  }

  if (any(table(vars) > 1)) {
    rlang::abort("`vars` should have unique members")
  }
  if (any(!(vars %in% colnames(s_data)))) {
    rlang::abort("1+ elements of `vars` are not in `x`")
  }

  s_data <- s_data[, vars]

  var_info <- tibble::tibble(variable = vars)

  ## Add types
  var_info <- dplyr::full_join(get_types(s_data), var_info, by = "variable")
  var_info$source <- "original"

  ## Return final object of class `recipe`
  out <- list(
    phyloseq = x,
    var_info = var_info,
    term_info = var_info,
    steps = NULL
  )
  class(out) <- "recipe"
  out
}


#' Print a Recipe
#'
#' @aliases print.recipe
#' @param x A `recipe` object
#' @param form_width The number of characters used to print the variables or
#'   terms in a formula
#' @param ... further arguments passed to or from other methods (not currently
#'   used).
#' @return The original object (invisibly)
#'
#' @export
print.recipe <- function(x, form_width = 30, ...) {
  cat("Recipe\n\n")
  cat("Inputs:\n\n")
  cat(" ", nrow(x$var_info), "variables\n")
  if (!is.null(x$steps)) {
    cat("\nOperations:\n\n")
    for (i in seq_along(x$steps)) {
      print(x$steps[[i]], form_width = form_width)
    }
  }
  invisible(x)
}
