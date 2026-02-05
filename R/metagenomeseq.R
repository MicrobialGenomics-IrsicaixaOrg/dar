#' @title (DEFUNCT) Zero-inflated Log-Normal mixture model for differential abundance
#'
#' @description
#' \strong{WARNING:} This function is **DEFUNCT** and has been removed from the package.
#' It is no longer maintained.
#'
#' @param rec A Recipe object.
#' @param ... Any other arguments (ignored).
#'
#' @export
#' @keywords internal
step_metagenomeseq <- function(rec, ...) {
  .Defunct(
    new = "step_maaslin",
    msg = "step_metagenomeseq is defunct and has been removed."
  )
}