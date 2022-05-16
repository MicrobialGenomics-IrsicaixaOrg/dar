#' Update packages
#'
#' This will check to see if all required packages are installed.
#'
#' @param pkg A character string for the package being checked
#' @param ... Extra arguments to pass to [utils::install.packages()]
#' @return Nothing is returned but a message is printed to the
#'  console about which packages (if any) should be installed along
#'  with code to do so.
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#' recipes_pkg_check("recipes")
#' }
recipes_pkg_check <- function(pkg = NULL, ...) {
  good <- rep(TRUE, length(pkg))
  for (i in seq(along.with = pkg)) {
    x <- stringr::str_remove_all(pkg[i], ".*[/]|.*[:]")
    tested <- try(find.package(x), silent = TRUE)
    if (class(tested)[1] == "try-error") {
      good[i] <- FALSE
    }
  }
  if (any(!good)) {
    pkList <- paste(stringr::str_remove_all(pkg[!good], ".*[/]|.*[:]") , collapse = ", ")
    inst <- stringr::str_c(pkg[!good], collapse = '", "')
    inst <- glue::glue('pak::pkg_install(c("{inst}"))')
    msg <- glue::glue(
      '{sum(!good)} {ifelse(sum(!good) > 1, "packages are", "package is")} needed for this, ',
      'step and {ifelse(sum(!good) > 1, " are", " is")} not installed: ',
      '({crayon::blue(pkList)}). Start a clean R session then run: {crayon::blue(inst)}'
    )
    cat(msg)
  }

  invisible()
}
