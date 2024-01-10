#' Filter taxa by rarity
#'
#' This is a convenience function around the filter_taxa function. It is
#' designed to speed up filtering complex experimental objects with one function
#' call. In the case of run_filter_by_rarity, the filtering will be based on the
#' rarity of each taxon. The taxa retained in the dataset are those where the
#' sum of their rarity is less than the provided threshold.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param threshold The rarity threshold for filtering taxa, expressed as a
#'   proportion of the total number of samples. For example, a threshold of 0.01
#'   means that a taxon must be present in less than 1% of the samples to be
#'   retained. The default value is 0.01.
#' @param id A character string that is unique to this step to identify it.
#'
#' @note This function modifies `rec` in place, you might want to make a copy of
#'   `rec` before modifying it if you need to preserve the original object.
#' @details The function calculates the rarity of all taxa in the phyloseq
#'   object as the proportion of samples in which they are present. It then
#'   compares this rarity to the threshold. If a taxon's rarity is greater than
#'   the threshold, that taxon is removed from the phyloseq object.
#' @return A Recipe object that has been filtered based on rarity.
#' @seealso \code{\link[phyloseq]{filter_taxa}}
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_by_rarity
#' @export
#' @autoglobal
#' @tests
#' data(test_prep_rec)
#' expect_error(step_filter_by_rarity(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Phylum")
#' rec
#'
#' ## Define step_filter_by_rarity step with default parameters
#' rec <- step_filter_by_rarity(rec, threshold = 0.01)
#' rec
methods::setGeneric(
  name = "step_filter_by_rarity",
  def = function(rec, 
                 threshold = 0.01, 
                 id = rand_id("filter_by_rarity")) {
    standardGeneric("step_filter_by_rarity")
  }
)

#' @rdname step_filter_by_rarity
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_rarity",
  signature = c(rec = "Recipe"),
  definition = function(rec, threshold = 0.01, id) {
    recipes_pkg_check(
      required_pkgs_filter_by_rarity(),
      "step_filter_by_rarity()"
    )
    add_step(
      rec,
      step_filter_by_rarity_new(threshold = threshold, id = id)
    )
  }
)

#' @rdname step_filter_by_rarity
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_rarity",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, threshold = 0.01, id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_filter_by_rarity_new <- function(threshold = 0.01, id) {
  step(subclass = "filter_by_rarity", threshold = threshold, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_filter_by_rarity <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_filter_by_rarity <- function(rec, threshold = 0.01) {
  rec@phyloseq <- 
    get_phy(rec) %>%
    phyloseq::filter_taxa(function(x) sum(x > 0) < (threshold * length(x)), TRUE) 
  
  rec
}
