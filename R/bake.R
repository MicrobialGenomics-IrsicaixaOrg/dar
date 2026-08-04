#' Define consensus strategies from a Recipe
#'
#' For a prep Recipe adds a consensus strategy to use for result extraction.
#'
#' @param rec A PrepRecipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param count_cutoff Indicates the minimum number of methods in which an OTU
#'   must be present (Default: NULL). If `count_cutoff` is `NULL`,
#'   it is set to `length(steps_ids(rec, "da")) - length(exclude)`.
#' @param weights Named numeric vector with the weight for each method.
#' @param exclude Method ids to exclude.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Bake steps
#' @aliases bake
#' @return An object of class `PrepRecipe`
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#'
#' ## Test bake default values
#' testthat::expect_error(bake(test_rec))
#' testthat::expect_s4_class(bake(test_prep_rec), "PrepRecipe")
#'
#' ## Test count_cutoff parameter
#' testthat::expect_error(cool(bake(test_prep_rec, count_cutoff = 10)))
#'
#' ## Test exclude parameter
#' test <- bake(test_prep_rec, exclude = steps_ids(test_prep_rec, "da")[1])
#' cool(test) |>
#'   testthat::expect_s3_class("tbl_df") |>
#'   testthat::expect_snapshot()
#'
#' nrow(cool(test)) |>
#'   testthat::expect_equal(27) |>
#'   testthat::expect_snapshot()
#'
#' test <- bake(test_prep_rec, exclude = "force_error")
#' testthat::expect_error(cool(test))
#'
#' ## Test weights parameter
#' weights <- c(2, 1)
#' names(weights) <- steps_ids(test_prep_rec, "da")
#' test <- bake(test_prep_rec, weights = weights)
#' cool(test) |>
#'   testthat::expect_s3_class("tbl_df") |>
#'   testthat::expect_snapshot()
#'
#' nrow(cool(test)) |>
#'   testthat::expect_equal(41) |>
#'   testthat::expect_snapshot()
#'
#' @examples
#' data(test_prep_rec)
#' rec <- test_prep_rec
#'
#' ## Default bake extracts common OTUs in all DA tested methods
#' ## (In this case the Recipe contains 3 methods)
#' res <- bake(rec)
#' cool(res)
#'
#' ## bake() requires a PrepRecipe
#' data(test_rec)
#' err <- testthat::expect_error(bake(test_rec))
#' err
#'
#' ## We can use the parameter `count_cutoff` to select those OTUs
#' ## shared by at least two methods
#' res <- bake(rec, count_cutoff = 2)
#' cool(res)
#'
#' ## Furthermore, we can exclude methods from the consensus strategy
#' ## with the `exclude` parameter.
#' res <- bake(rec, exclude = steps_ids(rec, "da")[1])
#' cool(res)
#'
#' ## Finally, we can use the `weights` parameter to weigh each method.
#' weights <- c(2, 1)
#' names(weights) <- steps_ids(rec, "da")
#' res <- bake(rec, weights = weights)
#' cool(res)
bake <- function(rec,
                 count_cutoff = NULL,
                 weights = NULL,
                 exclude = NULL,
                 id = rand_id("bake")) {

  check_prep_recipe(rec)
  recipes_pkg_check(required_pkgs_bake(), "bake()")

  add_step(
    rec,
    bake_new(
      count_cutoff = count_cutoff,
      weights = weights,
      exclude = exclude,
      id = id
    )
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
bake_new <- function(count_cutoff, weights, exclude, id) {
  step(
    subclass = "bake",
    count_cutoff = count_cutoff,
    weights = weights,
    exclude = exclude,
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_bake <- function(x, ...) { c() }

#' @noRd
#' @keywords internal
#' @autoglobal
run_bake <- function(rec, count_cutoff, weights, exclude, id) {

  da_ids <- steps_ids(rec, "da")
  specification <- validate_bake_spec(da_ids, count_cutoff, weights, exclude)
  ids <- specification$ids
  exclude <- specification$exclude
  count_cutoff <- specification$count_cutoff
  weights <- specification$weights
  df_weights <- tibble::enframe(weights, name = "method", value = "ponderation")

  significant <- .all_significant(rec, steps = ids) %>%
    dplyr::left_join(df_weights, by = "method")

  if (is.null(get_model(rec))) {
    res <- significant %>%
      dplyr::distinct(.data$taxa_id, .data$taxa, .data$method, .keep_all = TRUE) %>%
      dplyr::group_by(.data$taxa_id, .data$taxa) %>%
      dplyr::summarise(
        method_count = sum(.data$ponderation),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$method_count >= count_cutoff) %>%
      dplyr::distinct(.data$taxa_id, .data$taxa)
  } else {
    res <- significant %>%
      dplyr::distinct(
        .data$taxa_id, .data$taxa, .data$contrast_id, .data$comparison,
        .data$contrast_type, .data$var, .data$effect, .data$method,
        .keep_all = TRUE
      ) %>%
      dplyr::group_by(
        .data$taxa_id, .data$taxa, .data$contrast_id, .data$comparison,
        .data$contrast_type, .data$var, .data$effect
      ) %>%
      dplyr::summarise(
        method_count = sum(.data$ponderation),
        methods = paste(sort(unique(.data$method)), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$method_count >= count_cutoff)
  }

  cli::cli_inform(c(
    "i" = "Baking with {.field count_cutoff =} {.val {count_cutoff}}",
    if (length(exclude) > 0) { "i" = "Excluding {.val {exclude}}." }
  ))

  res
}

#' @noRd
validate_bake_spec <- function(da_ids, count_cutoff, weights, exclude) {
  exclude <- exclude %||% character()
  if (!is.character(exclude) || anyNA(exclude) || any(!nzchar(exclude)) ||
      anyDuplicated(exclude) || !all(exclude %in% da_ids)) {
    cli::cli_abort(
      "{.arg exclude} must contain unique executable DA step IDs.",
      class = "dar_error_invalid_step_id"
    )
  }
  ids <- setdiff(da_ids, exclude)
  if (length(ids) == 0L) {
    cli::cli_abort(
      "At least one executable DA method must remain after exclusions.",
      class = "dar_error_invalid_consensus"
    )
  }
  if (is.null(weights)) {
    weights <- stats::setNames(rep(1, length(ids)), ids)
  } else if (!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights)) ||
             any(weights < 0) || is.null(names(weights)) ||
             any(!nzchar(names(weights))) || anyDuplicated(names(weights)) ||
             !setequal(names(weights), ids) || length(weights) != length(ids)) {
    cli::cli_abort(
      paste0(
        "{.arg weights} must be a finite, non-negative, uniquely named numeric ",
        "vector exactly matching the non-excluded executable DA methods."
      ),
      class = "dar_error_invalid_weights"
    )
  }
  weights <- weights[ids]
  maximum <- sum(weights)
  if (is.null(count_cutoff)) count_cutoff <- length(ids)
  if (!is.numeric(count_cutoff) || length(count_cutoff) != 1L ||
      is.na(count_cutoff) || !is.finite(count_cutoff) || count_cutoff <= 0 ||
      count_cutoff > maximum) {
    cli::cli_abort(
      "{.arg count_cutoff} must be one finite positive value no greater than the available consensus weight ({maximum}).",
      class = "dar_error_invalid_count_cutoff"
    )
  }
  list(
    ids = ids, exclude = exclude, weights = weights,
    count_cutoff = count_cutoff
  )
}
