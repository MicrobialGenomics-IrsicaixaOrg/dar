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
#' @keywords internal
step <- function(subclass, ..., .prefix = "step_") {
  structure(list(...), class = c(paste0(.prefix, subclass), "step"))
}

#' @rdname step
#' @keywords internal
check <- function(subclass, ..., .prefix = "check_") {
  structure(list(...), class = c(paste0(.prefix, subclass), "check"))
}

#' Add a New Operation to the Current Recipe
#'
#' `add_step` adds a step to the last location in the recipe. `add_check` does
#' the same for checks.
#'
#' @param rec A [recipe()].
#' @param object A step or check object.
#' @return A updated [recipe()] with the new operation in the last slot.
#' @keywords internal
methods::setGeneric("add_step", function(rec, object)
  standardGeneric("add_step"))

#' @rdname add_step
#' @keywords internal
methods::setMethod(
  f = "add_step",
  signature = c("recipe"),
  definition = function(rec, object) {

    dupl_rec <-
      rec@steps %>%
      purrr::map_lgl(function(.x) {
        .x["id"] <- NULL
        object["id"] <- NULL
        testthat::compare(.x, object)$equal
      }) %>%
      any()

    if (dupl_rec) {
      expr <- 
        step_to_expr(object) %>% 
        stringr::str_replace('run', 'step') %>% 
        crayon::blue()
      
      text <- stringr::str_c(
        "This step is already defined with the same ",  
        "parameters and will be skipped: ", 
        collapse = ""
      )
      rlang::inform(c("!" = text, glue::glue("{expr}")))
    } else {
      rec@steps[[length(rec@steps) + 1]] <- object
    }

    rec
  }
)

#' @rdname add_step
#' @keywords internal
methods::setMethod(
  f = "add_step",
  signature = c("prep_recipe"),
  definition = function(rec, object) {

    dupl_rec <-
      rec@bakes %>%
      purrr::map_lgl(function(.x) {
        .x["id"] <- NULL
        object["id"] <- NULL
        testthat::compare(.x, object)$equal
      }) %>%
      any()

    if (dupl_rec) {
      expr <- 
        step_to_expr(object) %>% 
        stringr::str_replace('run', 'step') %>% 
        crayon::blue()
      
      text <- stringr::str_c(
        "This step is already defined with the same ",  
        "parameters and will be skipped: ", 
        collapse = ""
      )
      rlang::inform(c("!" = text, glue::glue("{expr}")))
    } else {
      rec@bakes[[length(rec@bakes) + 1]] <- object
    }

    rec
  }
)
