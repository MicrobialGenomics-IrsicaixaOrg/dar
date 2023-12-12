#' Filter taxa by prevalence
#'
#' This is a convenience function around the filter_taxa function. It is
#' designed to speed up filtering complex experimental objects with one function
#' call. In the case of run_filter_by_prevalence, the filtering will be based on
#' the prevalence of each taxon. The taxa retained in the dataset are those
#' where the prevalence is greater than the provided threshold.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param threshold The prevalence threshold for filtering taxa, expressed as a
#'   proportion of the total number of samples. For example, a threshold of 0.01
#'   means that a taxon must be present in at least 1% of the samples to be
#'   retained. The default value is 0.01.
#' @param id A character string that is unique to this step to identify it.
#'
#' @note This function modifies `rec` in place, you might want to make a copy of
#'   `rec` before modifying it if you need to preserve the original object.
#' @details The function calculates the prevalence of all taxa in the phyloseq
#'   object as the proportion of samples in which they are present. It then
#'   compares this prevalence to the threshold. If a taxon's prevalence is less
#'   than the threshold, that taxon is removed from the phyloseq object.
#' @return A Recipe object that has been filtered based on prevalence.
#' @seealso \code{\link[phyloseq]{filter_taxa}}
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_by_prevalence
#' @export
#' @autoglobal
#' @tests testthat
#' data(test_prep_rec)
#' expect_error(step_filter_by_prevalence(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#'
#' ## Define step_filter_by_prevalence step with default parameters
#' rec <- step_filter_by_prevalence(rec, threshold = 0.01)
#' rec
methods::setGeneric(
  name = "step_filter_by_prevalence",
  def = function(rec, 
                 threshold = 0.01, 
                 id = rand_id("filter_by_prevalence")) {
    standardGeneric("step_filter_by_prevalence")
  }
)

#' @rdname step_filter_by_prevalence
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_prevalence",
  signature = c(rec = "Recipe"),
  definition = function(rec, threshold = 0.01, id) {
    recipes_pkg_check(
      required_pkgs_filter_by_prevalence(),
      "step_filter_by_prevalence()"
    )
    add_step(
      rec,
      step_filter_by_prevalence_new(threshold = threshold, id = id)
    )
  }
)

#' @rdname step_filter_by_prevalence
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_prevalence",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, threshold = 0.01, id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_filter_by_prevalence_new <- function(threshold = 0.01, id) {
  step(subclass = "filter_by_prevalence", threshold = threshold, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_filter_by_prevalence <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_filter_by_prevalence <- function(rec, threshold = 0.01) {
  rec@phyloseq <- 
    get_phy(rec) %>%
    phyloseq::filter_taxa(function(x) sum(x > 0) > (threshold * length(x)), TRUE) 
  
  rec
}
