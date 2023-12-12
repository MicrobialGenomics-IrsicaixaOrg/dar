#' Define consensus strategies from a Recipe
#'
#' For a prep Recipe adds a consensus strategies to use for result extraction. 
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param count_cutoff Indicates the minimum number of methods in which an OTU
#'   must be present (Default: NULL). If count_cutoff is NULL count_cutoff is
#'   equal to `length(steps_ids(rec, "da")) - length(exclude)`
#' @param weights Named vector with the ponderation value for each method.
#' @param exclude Method ids to exclude.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Bake steps
#' @aliases bake
#' @return An object of class `PrepRecipe`
#' @export
#' @autoglobal
#' @tests testthat
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
#' testthat::expect_s3_class(cool(test), "tbl_df")
#' testthat::expect_equal(nrow(cool(test)), 34)
#' 
#' test <- bake(test_prep_rec, exclude = "force_error")
#' testthat::expect_error(cool(test))
#' 
#' ## Test weights parameter
#' weights <- c(2, 1, 3)
#' names(weights) <- steps_ids(test_prep_rec, "da")
#' test <- bake(test_prep_rec, weights = weights)
#' testthat::expect_s3_class(cool(test), "tbl_df")
#' testthat::expect_equal(nrow(cool(test)), 34)
#' 
#' weights <- c(2, 1)
#' names(weights) <- steps_ids(test_prep_rec, "da")[1:2]
#' test <- bake(test_prep_rec, weights = weights)
#' expect_error(cool(test))
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
#' ## bake and cool methods needs a PrepRecipe. If you pass a non-PrepRecipe
#' ## the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(bake(test_rec))
#' err
#' 
#' ## We can use the parameter `cout_cutoff` to for example select those OTUs
#' ## shared with at least two methods
#' res <- bake(rec, count_cutoff = 2)
#' cool(res)
#' 
#' ## Furthermore, we can exclude methods from the consensus strategy with the 
#' ## `exclude` parameter.
#' res <- bake(rec, exclude = steps_ids(rec, "da")[1])
#' cool(res)
#' 
#' ## Finally, we can use the `weights` parameter to weigh each method.
#' weights <- c(2, 1, 1)
#' names(weights) <- steps_ids(rec, "da")
#' res <- bake(rec, weights = weights)
#' cool(res)
methods::setGeneric(
  name = "bake",
  def = function(rec,
                 count_cutoff = NULL,
                 weights = NULL,
                 exclude = NULL,
                 id = rand_id("bake")) {
    standardGeneric("bake")
  }
)

#' @rdname bake
#' @export
#' @autoglobal
methods::setMethod(
  f = "bake",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, count_cutoff, weights, exclude, id) {
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
)

#' @rdname bake
#' @export
#' @autoglobal
methods::setMethod(
  f = "bake",
  signature = "Recipe",
  definition = function(rec, count_cutoff, weights, exclude, id) {
    text <- crayon::bgMagenta('prep(rec)')
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {text} and then try with {crayon::bgMagenta('bake()')}"
      )
    ))
  }
)

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

  if (!all(exclude %in% steps_ids(rec, "da"))) {
    msg <- 
      exclude[!exclude %in% steps_ids(rec, "da")] %>% 
      stringr::str_c(collapse = ", ")
    
    rlang::abort(glue::glue(
      "{crayon::blue(msg)} must be an existent step_id inside the PrepRecipe!"
    ))
  }
  ids <- steps_ids(rec, type = "da") %>% .[!. %in% exclude]
  if (is.null(count_cutoff)) {
    count_cutoff <- length(steps_ids(rec, "da")) - length(exclude)
  }
  
  if (count_cutoff > length(steps_ids(rec, "da"))) {
    rlang::abort(
      glue::glue(
        "{crayon::blue('count_cutoff')} must be less than or equal to ",
        "the number of da methods!"
      )
    )
  }

  if (is.null(weights)) {
    ids %>%
      purrr::set_names() %>%
      purrr::map_chr( ~ as.character(1)) %>%
      as.numeric()

    weights <- rep(1, length(ids))
    names(weights) <- ids
  }
  
  df_weights <- tibble::enframe(weights, name = "method", value = "ponderation")

  not_weights <- ids[!ids %in% names(weights)]
  if (length(not_weights) > 0) {
    not_weights <- not_weights %>% stringr::str_c(collapse = ", ")
    text_1 <-
      "Please include a value for this/these step/s in the vector of weights."
    text_2 <- glue::glue(
      "Alternatively, explicitly exclude the method/s via the ", 
      "exclude parameter."
    )
    rlang::abort(c(
      glue::glue(
        "Some non-excluded methods are not present in the weights vector ",
        "({crayon::yellow(not_weights)})"
      ),
      text_1,
      text_2
    ))
  }

  # weighted matrix
  res <- .all_significant(rec) %>% 
    dplyr::filter(method %in% ids) %>% 
    dplyr::left_join(df_weights, by = "method") %>% 
    dplyr::mutate(method_count = method_count * ponderation) %>% 
    dplyr::filter(method_count >= count_cutoff) %>% 
    dplyr::distinct(taxa_id, taxa)

  not_incl <- NULL
  if (!is.null(exclude)) {
    coll <- stringr::str_c(exclude, collapse = ", ")
    not_incl <- glue::glue("Results from {crayon::blue(coll)} are excluded")
  }

  rlang::inform(c(
    i = glue::glue(
      "Bake for {crayon::blue('count_cutoff =')} {crayon::blue(count_cutoff)}"
    ),
    i = not_incl,
    ""
  ))

  res
}
