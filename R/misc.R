#' Make a random identification field for steps
#'
#' @param prefix A single character string
#' @return A character string with the prefix and random letters separated by
#'   and underscore.
#'
#' @export
#' @autoglobal
#' @tests
#' set.seed(123); first_id <- rand_id()
#' set.seed(123); expect_identical(rand_id(), first_id)
#' expect_match(first_id, "^step__"); expect_true(sub("^step__", "", first_id) %in% pastry_df$id)
#' @examples
#' rand_id("step")
rand_id <- function(prefix = "step") {
  candidate <- sample(x = pastry_df$id, size = 1)
  paste(prefix, candidate, sep = "__")
}

#' Wrapper to convert phyloseq slots to tibble
#'
#' @param df output of `otu_table()`, `sample_data()` or `tax_table()` phyloseq
#'   functions.
#' @param id_name Name of the new column generated from rownames
#'
#' @return tibble
#' @keywords internal
#' @autoglobal
#'
#' @tests
#' data("metaHIV_phy")
#' phy <- metaHIV_phy
#' df <- phyloseq::otu_table(phy)
#'
#' test_1 <- to_tibble(df)
#' test_2 <- to_tibble(df, id_name = "random_id")
#'
#' expect_length(test_1, 157)
#' expect_s3_class(test_1, "tbl_df")
#' expect_s3_class(test_2, "tbl_df")
#' expect_equal(names(test_1)[1], "otu_id")
#' expect_equal(names(test_2)[1], "random_id")
#'
#' @examples
#' data(metaHIV_phy)
#' dar:::to_tibble(phyloseq::otu_table(metaHIV_phy))
to_tibble <- function(df, id_name = "otu_id") {
  df %>%
    data.frame(check.names = FALSE) %>%
    tibble::as_tibble(rownames = id_name)
}

#' Convert a step object into an unevaluated R call
#'
#' @description
#' This function takes a step from a `Recipe` and converts it into a native R
#' call (language object) using `rlang::call2`. This avoids the pitfalls of
#' pasting strings and parsing them, handling complex R objects like functions
#' or formulas naturally. The resulting call will have `rec` as its first
#' argument, followed by the parameters defined in the step.
#'
#' @param step A list or `step` class object containing the parameters for the
#'   operation, including an `id` string (e.g., `"method__randomstring"`).
#'
#' @return An unevaluated call (language object).
#'
#' @noRd
#' @keywords internal
#' @autoglobal
#' @tests
#' # 1. Test standard string and numeric parameters
#' step_standard <- step("maaslin",
#'   id = "maaslin__123",
#'   transform = "LOG",
#'   min_abundance = 0.1
#' )
#' call_standard <- step_to_call(step_standard)
#'
#' expect_true(is.call(call_standard))
#' expect_true(is.function(call_standard[[1]]))
#' expect_equal(call_standard[[1]], run_maaslin)
#' expect_equal(rlang::call_args(call_standard)$transform, "LOG")
#' expect_equal(rlang::call_args(call_standard)$min_abundance, 0.1)
#'
#' # 2. Test with a function parameter (The main reason for this refactor)
#' my_fun <- function(x) sum(x > 0) >= (0.03 * length(x))
#' step_func <- step("filter_taxa",
#'   id = "filter_taxa__abc",
#'   .f = my_fun
#' )
#' call_func <- step_to_call(step_func)
#'
#' expect_true(is.call(call_func))
#' expect_true(is.function(call_func[[1]]))
#' expect_equal(call_func[[1]], run_filter_taxa)
#' expect_true(is.function(rlang::call_args(call_func)$.f))
#'
#' # 3. Test with a formula parameter
#' step_formula <- step("deseq",
#'   id = "deseq__xyz",
#'   design = ~ RiskGroup2
#' )
#' call_formula <- step_to_call(step_formula)
#'
#' expect_true(is.call(call_formula))
#' expect_true(is.function(call_formula[[1]]))
#' expect_equal(call_formula[[1]], run_deseq)
#' expect_true(inherits(rlang::call_args(call_formula)$design, "formula"))
step_to_call <- function(step) {
  fn_obj <- step_runner(step)
  args <- step
  if (!"id" %in% names(formals(fn_obj))) { args[["id"]] <- NULL }
  rlang::call2(fn_obj, rec = quote(rec), !!!args)
}

#' Format a configured step without producing executable source text
#' @noRd
step_call_label <- function(step) {
  constructor <- paste0("step_", step_method(step))
  args <- unclass(step)
  args$id <- NULL
  if ("engine_args" %in% names(args) && length(args$engine_args) == 0L) {
    args$engine_args <- NULL
  }
  paste0("rec %>% ", rlang::expr_text(rlang::call2(constructor, !!!args)))
}

#' Finds common OTU between method results
#'
#' @param rec A Recipe object.
#' @param steps character vector with step ids to take in account
#' @param target Optional modeled target to include.
#' @param contrast_id Optional modeled contrast identifier to include.
#'
#' @return tibble
#' @export
#' @autoglobal
#'
#' @tests
#' data(test_prep_rec)
#' res_1 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")
#' )
#' res_2 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-1]
#' )
#' res_3 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-2]
#' )
#'
#' significant <- tidy_results(test_prep_rec, significant_only = TRUE) |>
#'   dplyr::mutate(effect = effect_direction(effect_size))
#' keys <- c("taxa_id", "contrast_id", "effect")
#' selected_steps <- list(
#'   steps_ids(test_prep_rec, "da"),
#'   steps_ids(test_prep_rec, "da")[-1],
#'   steps_ids(test_prep_rec, "da")[-2]
#' )
#' for (i in seq_along(selected_steps)) {
#'   expected <- significant |>
#'     dplyr::filter(step_id %in% selected_steps[[i]]) |>
#'     dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
#'     dplyr::summarise(sum_methods = dplyr::n_distinct(step_id), .groups = "drop") |>
#'     dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
#'   observed <- list(res_1, res_2, res_3)[[i]] |>
#'     dplyr::select(dplyr::all_of(keys), sum_methods) |>
#'     dplyr::arrange(dplyr::across(dplyr::all_of(keys)))
#'   expect_equal(observed, expected)
#' }
#' expect_s3_class(res_1, "tbl_df")
#'
#' @examples
#' data(test_prep_rec)
#'
#' ## Extract all intersections
#' find_intersections(test_prep_rec)
#'
#' ## Exclude the first method
#' find_intersections(
#'   test_prep_rec,
#'   steps = steps_ids(test_prep_rec, "da")[-1]
#' )
find_intersections <- function(rec, steps = steps_ids(rec, "da"),
                               target = NULL, contrast_id = NULL) {
  df <- intersection_df(
    rec,
    steps,
    target = target,
    contrast_id = contrast_id
  ) %>% tibble::as_tibble()
  df %>%
    tidyr::pivot_longer(
      dplyr::all_of(steps),
      names_to = "name",
      values_to = "value"
    ) %>%
    dplyr::filter(.data$value == 1) %>%
    dplyr::group_by(
      .data$taxa_id,
      .data$contrast_id,
      .data$comparison,
      .data$contrast_type,
      .data$var,
      .data$effect
    ) %>%
    dplyr::summarise(
      step_ids = paste(sort(unique(.data$name)), collapse = ", "),
      sum_methods = sum(.data$value),
      .groups = "drop"
    ) %>%
    dplyr::left_join(analysis_tax_table(rec), by = "taxa_id") %>%
    dplyr::arrange(dplyr::desc(.data$sum_methods))
}

#' Get step_ids from recipe
#'
#' @param rec A Recipe object.
#' @param type character vector indicating the type class. Options `c("all",
#'   "da", "prepro")`.
#' @param include_skipped For a prepared modeled recipe, include DA steps that
#'   were declared but skipped as statistically incompatible.
#'
#' @return character vector
#' @export
#' @autoglobal
#' @tests
#' data(test_prep_rec)
#' expect_snapshot(
#'   print(test_prep_rec),
#'   transform = function(lines) trimws(lines, which = "right")
#' )
#' rec <- test_prep_rec
#' expect_equal(
#'   steps_ids(rec),
#'   c("subset_taxa__Viennoiserie",
#'     "filter_taxa__Karakudamono",
#'     "maaslin__Welsh_cake",
#'     "deseq__Coussin_de_Lyon"
#'    )
#' )
#' expect_equal(
#'   steps_ids(rec, "da"),
#'   c("maaslin__Welsh_cake", "deseq__Coussin_de_Lyon")
#' )
#' expect_equal(
#'   steps_ids(rec, "prepro"),
#'   c("subset_taxa__Viennoiserie", "filter_taxa__Karakudamono")
#' )
#' expect_error(steps_ids(rec, "das"))
#' expect_type(steps_ids(rec), "character")
#'
#' @examples
#' data(test_rec)
#'
#' steps_ids(test_rec)
#'
#' ## Select differential-abundance or preprocessing steps
#' steps_ids(test_rec, type = "da")
#' steps_ids(test_rec, type = "prepro")
steps_ids <- function(rec, type = "all", include_skipped = FALSE) {
  check_any_recipe(rec)
  if (!type %in% c("all", "da", "prepro")) {
    cli::cli_abort(
      c(
        "x" = "Incorrect step type!",
        "i" = "Please use one of: {.val {c('all', 'da', 'prepro')}}."
      ),
      class = "dar_error_invalid_step_type"
    )
  }

  selected_steps <- switch(
    type,
    "all" = rec@steps,
    "da" = purrr::keep(rec@steps, is_da_step),
    "prepro" = purrr::keep(rec@steps, is_preprocessing_step)
  )
  out <- purrr::map_chr(selected_steps, ~ .x[["id"]])
  if (methods::is(rec, "PrepRecipe") && !include_skipped &&
      length(rec@execution) > 0L && type %in% c("all", "da")) {
    skipped <- rec@execution$skipped_steps$step_id %||% character()
    out <- setdiff(out, skipped)
  }
  out
}

#' @keywords internal
#' @autoglobal
tick <- function() { crayon::green(cli::symbol$tick) }
#' @keywords internal
#' @autoglobal
cross <- function() { crayon::red(cli::symbol$cross) }
#' @keywords internal
#' @autoglobal
info <- function() { crayon::blue(cli::symbol$info) }
#' @keywords internal
#' @autoglobal
dot <- function() { crayon::blue(cli::symbol$circle_filled) }
#' Checks if Recipe contains a rarefaction step
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#'
#' @export
#' @autoglobal
#' @return boolean
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' contains_rarefaction(rec)
#' contains_rarefaction(step_rarefaction(rec))
contains_rarefaction <- function(rec) {
  check_any_recipe(rec)
  any(purrr::map_lgl(rec@steps, function(configured_step) {
    identical(step_method(configured_step), "rarefaction")
  }))
}


#' Perform Rarefaction on Phyloseq Object
#'
#' This function performs rarefaction on a phyloseq object if the `rarefy`
#' parameter is set to TRUE. Rarefaction is a process that randomly subsamples
#' the data to a specified depth. This is done to account for differences in
#' sequencing depth between samples. However, this process is not without
#' controversy. Rarefaction can lead to loss of information and can also lead to
#' false positives in differential abundance testing. For more information, see
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2
#'
#' @param phy A phyloseq object.
#' @param rarefy A logical value indicating whether to perform rarefaction. If
#'   'no_seed', rarefaction is performed without a set seed. If FALSE, no
#'   rarefaction is performed.
#'
#' @return A phyloseq object after rarefaction if `rarefy` is TRUE or "no_seed",
#'   otherwise the original phyloseq object is returned.
#' @keywords internal
use_rarefy <- function(phy, rarefy) {
  if (isTRUE(rarefy)) {
    phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE)
  }
  if (rarefy == "no_seed") { phy <- phyloseq::rarefy_even_depth(phy, verbose = FALSE) }
  phy
}

#' @noRd
#' @keywords internal
#' @autoglobal
rarefy_msg <- function(steps) {
  info_rarefy <- "Rarefaction is a process that randomly subsamples the data..."
  rarefy_vals <- purrr::map(steps, ~ .x[["rarefy"]])

  has_true <- any(
    purrr::map_lgl(
      rarefy_vals,
      ~ isTRUE(.x) || identical(.x, as.name("T"))
    )
  )
  has_no_seed <- any(purrr::map_lgl(rarefy_vals, ~ identical(.x, "no_seed")))

  if (has_true) {
    cli::cli_inform(
      c(
        "!" = "Rarefaction is enabled (with set seed).",
        "i" = "This process is not without controversy. Use {.fun rarefaction_help} to get more info."
      )
    )
  }

  if (has_no_seed) {
    cli::cli_inform(
      c(
        "!" = "Rarefaction is enabled (without set seed).",
        "i" = "This process is not without controversy. Use {.fun rarefaction_help} to get more info."
      )
    )
  }

  steps
}


#' Information about the Rarefaction process
#'
#' Prints detailed information regarding the use of rarefaction in microbiome
#' analysis, its implications, and relevant literature.
#'
#' @return invisible
#' @export
#' @autoglobal
#' @examples
#' rarefaction_help()
rarefaction_help <- function() {
  cli::cli_rule(left = "{.strong Rarefaction in Microbiome Analysis}")
  cli::cli_text("")
  cli::cli_alert_info(c("{.strong What is Rarefaction?} It is a process that randomly subsamples the data to a specified depth to account for differences in sequencing depth between samples."))
  cli::cli_text("")
  cli::cli_alert_warning("{.strong The Controversy:}")
  cli::cli_bullets(
    c(
      "i" = "It can lead to an irreversible {.strong loss of information}.",
      "i" = "It can lead to {.strong false positives} in differential abundance testing."
    )
  )
  cli::cli_text("")
  cli::cli_alert_success("{.strong Reproducibility:} If performed with a {.emph set seed}, it ensures reproducible results, but removes randomness in the subsampling process.")
  cli::cli_alert_info("{.strong More info:} {.url https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2}")
  cli::cli_text("")
  invisible(TRUE)
}
