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
#' @tests
#' dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()") |> 
#'   expect_snapshot()
#'   
#' dar:::recipes_pkg_check(dar:::required_pkgs_aldex(), "step_aldex()") |> 
#'   expect_snapshot()
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
    
    to_inst <-
      pkg[!good] %>%
      stringr::str_remove_all("bioc::") %>%
      paste("\"", ., "\"", sep = "", collapse = ", ")

    inst <- glue::glue('BiocManager::install(c({to_inst}))')
    n_packages = glue::glue('{sum(!good)} {ifelse(sum(!good) > 1, "packages are", "package is")}')
    step_name <- glue::glue('{crayon::blue(step_name)}')
    to_be <- glue::glue('{ifelse(sum(!good) > 1, "are", "is")}')
    blue_pkList <- glue::glue('{crayon::blue(pkList)}')
    rlang::inform(c(
      "i" = glue::glue('{n_packages} needed for {step_name} and {to_be} not installed: ({blue_pkList})'),
      "*" = glue::glue('Start a clean R session then run: {crayon::blue(inst)}')
    ))
  }

  invisible()
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_error <- function(x, ...) { c("bioc::randompackage", "packrandom") }

#' Returns required pakcages for Recipe object
#'
#' @param rec A `Recipe` object
#'
#' @aliases required_deps
#' @return character
#' @keywords internal
#' @autoglobal
#' @examples 
#' data(test_rec)
#' 
#' ## The function returns instructions to install any uninstalled dependencies 
#' ## needed to run the Recipe steps
#' dar:::required_deps(test_rec)
#' 
#' ## The function also works with PrepRecipe-class objects
#' data(test_prep_rec)
#' dar:::required_deps(test_prep_rec)
methods::setGeneric("required_deps", function(rec)
  standardGeneric("required_deps"))

#' @rdname required_deps
#' @keywords internal
#' @autoglobal
methods::setMethod(
  f = "required_deps",
  signature = "Recipe",
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
