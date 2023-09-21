#' Update packages
#'
#' This will check to see if all required packages are installed.
#'
#' @param pkg A character string for the package being checked
#' @param step_name Name of the step.
#' @param ... Extra arguments to pass to [utils::install.packages()]
#' @return Nothing is returned but a message is printed to the console about
#'   which packages (if any) should be installed along with code to do so.
#' @keywords internal
#' @autoglobal
#' @tests testthat
#' expected <-
#'   glue::glue("1 package is needed for step_aldex() and is not installed: ", 
#'   "(randompackage). Start a clean R session then run: ", 
#'   "pak::pkg_install(c(\"bioc::randompackage\")) ")
#'   
#' res <- capture.output(
#'   dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()")
#' )
#' expect_equal(res, expected)
#' expect_invisible(
#'   dar:::recipes_pkg_check(dar:::required_pkgs_aldex(), "step_aldex()")
#' )
#' expect_invisible(
#'   dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()")
#' )
recipes_pkg_check <- function(pkg = NULL, step_name, ...) {
  good <- rep(TRUE, length(pkg))
  for (i in seq(along.with = pkg)) {
    x <- stringr::str_remove_all(pkg[i], ".*[/]|.*[:]")
    tested <- try(find.package(x), silent = TRUE)
    if (methods::is(tested, "try-error")) {
      good[i] <- FALSE
    }
  }
  
  if (any(!good)) {
    pkList <- paste(
      stringr::str_remove_all(pkg[!good], ".*[/]|.*[:]") , collapse = ", "
    )
    inst <- stringr::str_c(pkg[!good], collapse = '", "')
    inst <- glue::glue('pak::pkg_install(c("{inst}"))')
    msg <- glue::glue(
      '{sum(!good)} {ifelse(sum(!good) > 1, "packages are", "package is")} ',
      'needed for {crayon::blue(step_name)} and ',
      '{ifelse(sum(!good) > 1, "are", "is")} not installed: ',
      '({crayon::blue(pkList)}). Start a clean R session then run: ',
      '{crayon::blue(inst)}'
    )
    cat(c(msg, "\n"))
  }

  invisible()
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_error <- function(x, ...) { c("bioc::randompackage") }

#' Returns required pakcages for recipe object
#'
#' @param rec A `recipe` object
#'
#' @aliases required_deps
#' @return character
#' @export
#' @autoglobal
#' @examples 
#' data(test_rec)
#' 
#' ## The function returns instructions to install any uninstalled dependencies 
#' ## needed to run the recipe steps
#' required_deps(test_rec)
#' 
#' ## The function also works with prep_recipe-class objects
#' data(test_prep_rec)
#' required_deps(test_prep_rec)
methods::setGeneric("required_deps", function(rec)
  standardGeneric("required_deps"))

#' @rdname required_deps
#' @export
#' @autoglobal
methods::setMethod(
  f = "required_deps",
  signature = "recipe",
  definition = function(rec) {
    text <- 
      rec@steps %>%
      purrr::walk( ~ {
        id <- class(.x)[[1]]
        id_2 <- stringr::str_remove_all(id, "step_")
        glue::glue("recipes_pkg_check(required_pkgs_{id_2}(), '{id}()')") %>%
          parse(text = .) %>%
          eval()
      })
    
    invisible()
  }
)
