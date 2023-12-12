#' Filter taxa by abundance
#'
#' This is a convenience wrapper around the filter_taxa function. It is intended
#' to speed up filtering complex experimental objects with one function call. In
#' the case of filter_by_abundance, the filtering will be based on the relative
#' abundance of each taxon. The taxa retained in the dataset are those where the
#' sum of their abundance is greater than the product of the total abundance and
#' the provided threshold.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param threshold The relative abundance threshold for filtering taxa,
#'   expressed as a proportion of the total abundance. For example, a threshold
#'   of 0.01 means that a taxon must make up at least 1% of the total abundance
#'   to be retained. The default value is 0.01.
#' @param id A character string that is unique to this step to identify it.
#'
#' @note This function modifies `rec` in place, you might want to make a copy of
#'   `rec` before modifying it if you need to preserve the original object.
#' @details The function calculates the total abundance of all taxa in the
#'   phyloseq object. It then compares this total abundance to the abundance of
#'   each individual taxon. If a taxon's abundance is less than the threshold
#'   times the total abundance, that taxon is removed from the phyloseq object.
#' @return A Recipe object that has been filtered based on abundance.
#' @seealso \code{\link[phyloseq]{filter_taxa}}
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_by_abundance
#' @export
#' @autoglobal
#' @tests testthat
#' data(test_prep_rec)
#' expect_error(step_filter_by_abundance(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#'
#' ## Define filter_by_abundance step with default parameters
#' rec <- step_filter_by_abundance(rec, threshold = 0.01)
#' rec
methods::setGeneric(
  name = "step_filter_by_abundance",
  def = function(rec, 
                 threshold = 0.01, 
                 id = rand_id("filter_by_abundance")) {
    standardGeneric("step_filter_by_abundance")
  }
)

#' @rdname step_filter_by_abundance
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_abundance",
  signature = c(rec = "Recipe"),
  definition = function(rec, threshold = 0.01, id) {
    recipes_pkg_check(
      required_pkgs_filter_by_abundance(),
      "step_filter_by_abundance()"
    )
    add_step(
      rec,
      step_filter_by_abundance_new(threshold = threshold, id = id)
    )
  }
)

#' @rdname step_filter_by_abundance
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_abundance",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, threshold = 0.01, id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_filter_by_abundance_new <- function(threshold = 0.01, id) {
  step(subclass = "filter_by_abundance", threshold = threshold, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_filter_by_abundance <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_filter_by_abundance <- function(rec, threshold = 0.01) {
  t_abun <- sum(phyloseq::otu_table(get_phy(rec)))
  rec@phyloseq <- 
    get_phy(rec) %>%
    phyloseq::filter_taxa(function(x) sum(x) > (t_abun * threshold), TRUE) 
  
  rec
}
