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
  if (is.null(pkg) || length(pkg) == 0) return(invisible())
  pkg_names <- stringr::str_remove_all(pkg, ".*[/]|.*[:]")
  is_inst <- vapply(pkg_names, rlang::is_installed, logical(1))
  
  if (any(!is_inst)) {
    missing_pkgs <- pkg[!is_inst]
    missing_names <- pkg_names[!is_inst]

    to_inst <- 
      missing_pkgs %>%
      stringr::str_remove_all("bioc::") %>%
      paste0('"', ., '"', collapse = ", ")

    inst_code <- glue::glue("BiocManager::install(c({to_inst}))")
    cli::cli_inform(
      c(
        "i" = "{length(missing_names)} package{?s} {?is/are} needed for {.val {step_name}} and {?is/are} not installed: {.pkg {missing_names}}.",
        "*" = "Start a clean R session then run: {.run {inst_code}}"
      )
    )
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
#' @param steps Optional list of steps whose dependencies should be checked.
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
required_deps <- function(rec, steps = NULL) {
  check_any_recipe(rec)
  if (is.null(steps)) {
    steps <- rec@steps
  }
  purrr::walk(steps, function(step_obj) {
    id <- step_method(step_obj)
    req_fun <-  get0(paste0("required_pkgs_", id), mode = "function")
    if (!is.null(req_fun)) {
      recipes_pkg_check(pkg = req_fun(), step_name = paste0("step_", id, "()"))
    }
  })
    
  invisible()
}

#' Require an optional package for a plotting path
#' @noRd
require_optional_package <- function(pkg, caller) {
  if (rlang::is_installed(pkg)) {
    return(invisible(TRUE))
  }
  cli::cli_abort(
    c(
      "x" = "Package {.pkg {pkg}} is required by {.fun {caller}}.",
      "i" = "Install it with {.run BiocManager::install(\"{pkg}\")} and try again."
    ),
    class = "dar_error_missing_package",
    package = pkg,
    caller = caller
  )
}
