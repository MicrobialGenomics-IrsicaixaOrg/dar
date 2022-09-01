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
recipes_pkg_check <- function(pkg = NULL, step_name, ...) {
  good <- rep(TRUE, length(pkg))
  for (i in seq(along.with = pkg)) {
    x <- stringr::str_remove_all(pkg[i], ".*[/]|.*[:]")
    tested <- try(find.package(x), silent = TRUE)
    if (is(tested, "try-error")) {
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
required_pkgs_error <- function(x, ...) { c("bioc::randompackage") }

#' Returns required pakcages for recipe object
#'
#' @param rec A `recipe` object
#'
#' @aliases required_deps
#' @return character
#' @export
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

check_version <- function(pkg, v) {
  
  check <- utils::packageVersion(pkg) == v
  if (!check) {
    rlang::warn(c(
     "i" = glue::glue(
       "Temporarily the version of the package ", 
       "{crayon::blue(pkg)} must be {crayon::blue(paste0('v', v))}, but you ")
    ), use_cli_format = TRUE)
  }
}

