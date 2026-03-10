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

    # Format for BiocManager::install
    to_inst <- 
      missing_pkgs %>%
      stringr::str_remove_all("bioc::") %>%
      paste0('"', ., '"', collapse = ", ")

    inst_code <- glue::glue("BiocManager::install(c({to_inst}))")
    
    n_missing <- length(missing_pkgs)
    pkg_word <- ifelse(n_missing > 1, "packages are", "package is")
    to_be <- ifelse(n_missing > 1, "are", "is")
    missing_str <- paste(missing_names, collapse = ", ")

    msg <- c(
      "i" = glue::glue("{n_missing} {pkg_word} needed for {crayon::blue(step_name)} and {to_be} not installed: ({crayon::blue(missing_str)})"),
      "*" = glue::glue("Start a clean R session then run: {crayon::blue(inst_code)}")
    )
    
    rlang::inform(msg)
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
required_deps <- function(rec) {
  check_any_recipe(rec)
  purrr::walk(rec@steps, function(step_obj) {
    id <- stringr::str_remove_all(class(step_obj)[[1]], "step_")
    req_fun <-  get0(paste0("required_pkgs_", id), mode = "function")
    if (!is.null(req_fun)) {
      recipes_pkg_check(pkg = req_fun(), step_name = paste0("step_", id, "()"))
    }
  })
    
  invisible()
}