#' Filter taxa based on across-sample OTU abundance criteria
#'
#' This function is directly analogous to the genefilter function for microarray
#' filtering, but is used for filtering OTUs from phyloseq objects. It applies
#' an arbitrary set of functions — as a function list, for instance, created by
#' filterfun — as across-sample criteria, one OTU at a time. It takes as input a
#' phyloseq object, and returns a logical vector indicating whether or not each
#' OTU passed the criteria. Alternatively, if the "prune" option is set to
#' FALSE, it returns the already-trimmed version of the phyloseq object.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param .f A function or list of functions that take a vector of abundance
#'   values and return a logical. Some canned useful function types are included
#'   in the genefilter-package.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_taxa
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' rec_1 <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0 * length(x))") |> 
#'   step_metagenomeseq(rm_zeros = 0)
#'   
#' rec_2 <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0 * length(x))") |> 
#'   step_metagenomeseq(rm_zeros = 0.01)
#'   
#' rec_3 <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0 * length(x))") |> 
#'   step_metagenomeseq(rm_zeros = NULL)
#'   
#' expect_error(prep(rec_1))
#' expect_match(capture_warnings(prep(rec_2)), "NaNs produced", all = TRUE)
#' expect_match(capture_warnings(prep(rec_3)), "NaNs produced", all = TRUE)
#' expect_s4_class(suppressWarnings(prep(rec_2)), "PrepRecipe")
#' expect_s4_class(suppressWarnings(prep(rec_3)), "PrepRecipe")
#'
#' data(test_prep_rec)
#' expect_error(step_filter_taxa(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#' 
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#' 
#' ## Define filter taxa step with default parameters
#' rec <- 
#'   step_filter_taxa(rec, .f = "function(x) sum(x > 0) >= (0.03 * length(x))")
#' rec
methods::setGeneric(
  name = "step_filter_taxa",
  def = function(rec, .f, id = rand_id("filter_taxa")) {
    standardGeneric("step_filter_taxa")
  }
)

#' @rdname step_filter_taxa
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_taxa",
  signature = c(rec = "Recipe"),
  definition = function(rec, .f, id) {
    recipes_pkg_check(required_pkgs_filter_taxa(), "step_filter_taxa()")
    add_step(rec, step_filter_taxa_new(.f = .f, id = id))
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_filter_taxa_new <- function(.f, id) {
  step(subclass = "filter_taxa", .f = .f, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_filter_taxa <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_filter_taxa <- function(rec, .f) {
  
  rm_zeros <- NULL
  if (any(stringr::str_detect(steps_ids(rec), "metagenomeseq"))) {
    rm_zeros <- rec@steps %>%
      purrr::pluck(
        which(stringr::str_detect(steps_ids(rec), "metagenomeseq")),
        "rm_zeros", 
        .default = NULL
      )
  }
  
  is_metagenomeseq <- TRUE
  if (is.null(rm_zeros)) { 
    rm_zeros <- 0 
    is_metagenomeseq <- FALSE
  }
  
  rec@phyloseq <- 
    phyloseq::filter_taxa(get_phy(rec), eval(parse(text = .f)), prune = TRUE)
  
  val <- 
    zero_otu(rec) %>% 
    dplyr::filter(pct == 0) %>% 
    nrow()
  
  if (val > 0 & rm_zeros == 0 & is_metagenomeseq) {
    rlang::abort(c(
      "!" = glue::glue(
        "{crayon::bgMagenta('step_filter_taxa()')} returns a phyloseq ", 
        "object that contains taxa with values of 0 in all samples of a ", 
        "level within the variable of interest. This can cause errors during ",
        "the execution of metagenomeseq method!"
      ),
      "*" = "Please create a new Recipe using a stricter filter expression.", 
      "*" = glue::glue(
        "Alternatively, you can increase the rm_zeros value ", 
        "{crayon::bgMagenta('step_metagenomeseq(rm_zeros = 0.01)')}. This ", 
        "value indicates the minimum proportion of samples of the same level ", 
        "with more than 0 counts."
      )
    ), use_cli_format = TRUE)
  }
  
  rec
}


## Extract outs with all 0 values in at least on level of the variable ----

#' Extract outs with all 0 values in at least on level of the variable 
#'
#' @param obj A `Recipe` or `phyloseq` object.
#' @param var Variable of interest. Must be present in the metadata.
#' @param pct_cutoff Minimum of pct counts samples with counts for each taxa. 
#'
#' @aliases zero_otu
#' @return character vector
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' 
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' 
#' ## Extract outs with all 0 values
#' zero_otu(rec)
methods::setGeneric(
  name = "zero_otu",
  def = function(obj, var = NULL, pct_cutoff = 0) {
    standardGeneric("zero_otu")
  }
)

#' @rdname zero_otu
#' @export
#' @autoglobal
methods::setMethod(
  f = "zero_otu", 
  signature = "Recipe", 
  definition = function(obj, var, pct_cutoff) {
    var <- get_var(obj)[[1]]
    otu_table(obj) %>% 
      tidyr::pivot_longer(-1, names_to = "sample_id") %>% 
      dplyr::left_join(sample_data(obj), by = "sample_id") %>% 
      dplyr::mutate(no_zero = ifelse(value == 0, 0, 1)) %>% 
      dplyr::group_by(taxa_id, !!dplyr::sym(var)) %>%
      dplyr::summarise(
        no_zero = sum(no_zero), 
        total = dplyr::n(), 
        pct = no_zero / total,
        .groups = "drop"
      ) %>%
      dplyr::arrange(pct) %>% 
      dplyr::filter(pct >= pct_cutoff)
  }
)

#' @rdname zero_otu
#' @export
#' @autoglobal
methods::setMethod(
  f = "zero_otu",
  signature = "phyloseq",
  definition = function(obj, var, pct_cutoff) {
    if (is.null(pct_cutoff)) { pct_cutoff <- 0 }
    phyloseq::otu_table(obj) %>%
      to_tibble("taxa_id") %>%
      tidyr::pivot_longer(-1, names_to = "sample_id") %>%
      dplyr::left_join(
        phyloseq::sample_data(obj) %>%
          to_tibble("sample_id") %>%
          dplyr::select(1, !!var),
        by = "sample_id"
      ) %>%
      dplyr::mutate(no_zero = ifelse(value == 0, 0, 1)) %>% 
      dplyr::group_by(taxa_id, !!dplyr::sym(var)) %>%
      dplyr::summarise(
        no_zero = sum(no_zero), 
        total = dplyr::n(), 
        pct = no_zero / total,
        .groups = "drop"
      ) %>%
      dplyr::arrange(pct) %>% 
      dplyr::filter(pct >= pct_cutoff)
  }
)
