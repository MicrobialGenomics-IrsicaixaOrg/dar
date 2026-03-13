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
  exclude <- exclude %||% character()
  invalid_steps <- exclude[!exclude %in% da_ids]
  if (length(invalid_steps) > 0) {
    cli::cli_abort(
      c(
        "x" = "{cli::qty(invalid_steps)}Step ID{?s} {.val {invalid_steps}} {?does/do} not exist in the {.cls PrepRecipe}.",
        "i" = "Check the available IDs with {.code steps_ids(rec, 'da')}."
      ),
      class = "dar_error_invalid_step_id"
    )
  }

  ids <- da_ids[!da_ids %in% exclude]

  if (is.null(count_cutoff)) { count_cutoff <- length(ids) }

  if (count_cutoff > length(da_ids)) {
    cli::cli_abort(
      c(
        "x" = "{.arg count_cutoff} must be less than or equal to the number of {.val da} methods.",
        "i" = "Current number of {.val da} methods: {.val {length(da_ids)}}."
      ),
      class = "dar_error_invalid_count_cutoff"
    )
  }

  if (is.null(weights)) { weights <- stats::setNames(rep(1, length(ids)), ids) }

  df_weights <- tibble::enframe(weights, name = "method", value = "ponderation")
  not_weights <- ids[!ids %in% names(weights)]
  if (length(not_weights) > 0) {
    cli::cli_abort(
      c(
        "x" = "{cli::qty(not_weights)}Some non-excluded method{?s} {?is/are} missing from the weights vector: {.val {not_weights}}.",
        "i" = "{cli::qty(not_weights)}Please provide a weight for {?this/these} method{?s}.",
        "i" = "{cli::qty(not_weights)}Alternatively, explicitly exclude {?it/them} via the {.arg exclude} argument."
      ),
      class = "dar_error_missing_weights"
    )
  }

  res <- 
    .all_significant(rec) %>% 
    dplyr::filter(.data$method %in% ids) %>% 
    dplyr::left_join(df_weights, by = "method") %>% 
    dplyr::mutate(method_count = .data$method_count * .data$ponderation) %>% 
    dplyr::filter(.data$method_count >= count_cutoff) %>% 
    dplyr::distinct(.data$taxa_id, .data$taxa)

  cli::cli_inform(c(
    "i" = "Baking with {.field count_cutoff =} {.val {count_cutoff}}",
    if (length(exclude) > 0) { "i" = "Excluding {.val {exclude}}." }
  ))

  res
}