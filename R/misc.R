## Buckets variables into discrete, mutally exclusive types
get_types <- function(x) {
  var_types <-
    c(
      character = "nominal",
      factor = "nominal",
      ordered = "nominal",
      integer = "numeric",
      numeric = "numeric",
      double = "numeric",
      Surv = "censored",
      logical = "logical",
      Date = "date",
      POSIXct = "date",
      list = "list",
      textrecipes_tokenlist = "tokenlist"
    )

  classes <- lapply(x, class)
  res <- lapply(classes,
                function(x, types) {
                  in_types <- x %in% names(types)
                  if (sum(in_types) > 0) {
                    # not sure what to do with multiple matches; right now
                    ## pick the first match which favors "factor" over "ordered"
                    out <-
                      unname(types[min(which(names(types) %in% x))])
                  } else {
                    out <- "other"
                  }
                  out
                },
                types = var_types)
  res <- unlist(res)
  tibble::tibble(variable = names(res), type = unname(res))
}

#' Make a random identification field for steps
#'
#' @export
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#' @keywords internal
rand_id <- function(prefix = "step", len = 5) {
  candidates <- c(letters, LETTERS, paste(0:9))
  paste(prefix,
        paste0(sample(candidates, len, replace = TRUE), collapse = ""),
        sep = "_"
  )
}
