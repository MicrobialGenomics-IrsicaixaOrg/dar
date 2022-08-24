#' Define consensus strategies from a recipe
#'
#' For a prep recipe adds a consensus strategies to use for result extraction. 
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
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
#' @return An object of class `prep_recipe`
#' @export
#' @examples 
#' data(test_prep_rec)
#' rec <- test_prep_rec
#' 
#' ## Default bake extracts common OTUs in all DA tested methods 
#' ## (In this case the recipe contains 3 methods)
#' res <- bake(rec)
#' cool(res)
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
methods::setMethod(
  f = "bake",
  signature = c(rec = "prep_recipe"),
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
methods::setMethod(
  f = "bake",
  signature = "recipe",
  definition = function(rec, count_cutoff, weights, exclude, id) {
    text <- crayon::bgMagenta('prep(rec)')
    rlang::abort(c(
      "This function needs a prep recipe!",
      glue::glue(
        "Run {text} and then try with {crayon::bgMagenta('bake()')}"
      )
    ))
  }
)

#' @noRd
#' @keywords internal
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
required_pkgs_bake <- function(x, ...) { c() }

#' @noRd
#' @keywords internal
run_bake <- function(rec, count_cutoff, weights, exclude, id) {

  ids <- steps_ids(rec, type = "da") %>% .[!. %in% exclude]
  if (is.null(count_cutoff)) {
    count_cutoff <- length(steps_ids(rec, "da")) - length(exclude)
  }

  if (is.null(weights)) {
    ids %>%
      purrr::set_names() %>%
      purrr::map_chr( ~ 1) %>%
      as.numeric()

    weights <- rep(1, length(ids))
    names(weights) <- ids
  }

  not_weights <- ids[!ids %in% names(weights)]
  if (length(not_weights) > 0) {
    not_weights <- not_weights %>% stringr::str_c(collapse = ", ")
    text_1 <-
      "Please include a value for this/these step/s in the vector of weights."
    text_2 <-
      "Alternatively, explicitly exclude the method/s via the exclude parameter."
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
  res <-
    ids %>%
    purrr::map2_dfc(seq_along(.), ~ {
      intersection_df(rec, ids) %>%
        tibble::as_tibble() %>%
        dplyr::select(!!.x) %>%
        dplyr::mutate(!!(.x) := (!!dplyr::sym(.x) * weights[.y]))
    }) %>%
    dplyr::bind_cols(intersection_df(rec, ids)[1], .) %>%
    tidyr::pivot_longer(cols = -1) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(
      step_ids = 
        purrr::map_chr(.data$name, ~ .x) %>% 
        stringr::str_c(collapse = ", "),
      sum_methods = sum(.data$value)
    ) %>%
    dplyr::right_join(tax_table(rec), ., by = "taxa_id") %>%
    dplyr::arrange(-.data$sum_methods) %>%
    dplyr::filter(.data$sum_methods >= count_cutoff) %>%
    dplyr::select(.data$taxa_id, .data$taxa)

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
