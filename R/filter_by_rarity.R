#' Filter taxa by rarity
#'
#' This is a convenience function around the filter_taxa function. It is
#' designed to speed up filtering complex experimental objects with one function
#' call. In the case of run_filter_by_rarity, the filtering will be based on the
#' rarity of each taxon. The taxa retained in the dataset are those where the
#' prevalence (proportion of samples where present) is LESS than the provided
#' threshold.
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
#' @details The function calculates the prevalence of all taxa in the phyloseq
#'   object as the proportion of samples in which they are present. It then
#'   compares this prevalence to the threshold. If a taxon's prevalence is greater
#'   than or equal to the threshold, that taxon is removed (filtered out) from 
#'   the phyloseq object, leaving only the "rare" taxa.
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
#' rec <- recipe(metaHIV_phy)
#' rec
#'
#' ## Define step_filter_by_rarity step with default parameters
#' rec <- step_filter_by_rarity(rec, threshold = 0.01)
#' rec
step_filter_by_rarity <- function(rec, 
                                  threshold = 0.01, 
                                  id = rand_id("filter_by_rarity")) {
  
  check_recipe(rec)
  recipes_pkg_check(
    required_pkgs_filter_by_rarity(), 
    "step_filter_by_rarity()"
  )
  
  add_step(
    rec,
    step(
      subclass = "filter_by_rarity", 
      threshold = threshold, 
      id = id
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_filter_by_rarity <- function(rec, threshold = 0.01, id) {
  rec@phyloseq <- 
    get_phy(rec) %>%
    phyloseq::filter_taxa(function(x) sum(x > 0) < (threshold * length(x)), TRUE) 
  
  rec
}

#' @noRd
#' @keywords internal
required_pkgs_filter_by_rarity <- function(x, ...) {  c("bioc::phyloseq") }
