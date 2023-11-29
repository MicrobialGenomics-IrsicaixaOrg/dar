#' Filter taxa by variance
#'
#' This is a convenience function around the filter_taxa function. It is
#' designed to speed up filtering complex experimental objects with one function
#' call. In the case of run_filter_by_variance, the filtering will be based on
#' the variance of each taxon. The taxa retained in the dataset are those where
#' the variance of their abundance is greater than the provided threshold.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param threshold The variance threshold for filtering taxa. The default value
#'   is 0.01.
#' @param id A character string that is unique to this step to identify it.
#'
#' @note This function modifies `rec` in place, you might want to make a copy of
#'   `rec` before modifying it if you need to preserve the original object.
#' @details The function calculates the variance of all taxa in the phyloseq
#'   object. It then compares this variance to the variance of each individual
#'   taxon. If a taxon's variance is less than the threshold, that taxon is
#'   removed from the phyloseq object.
#' @return A Recipe object that has been filtered based on variance.
#' @seealso \code{\link[phyloseq]{filter_taxa}}
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_by_variance
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#'
#' ## Define step_filter_by_variance step with default parameters
#' rec <- step_filter_by_variance(rec, threshold = 0.01)
#' rec
methods::setGeneric(
  name = "step_filter_by_variance",
  def = function(rec, 
                 threshold = 0.01, 
                 id = rand_id("filter_by_variance")) {
    standardGeneric("step_filter_by_variance")
  }
)

#' @rdname step_filter_by_variance
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_variance",
  signature = c(rec = "Recipe"),
  definition = function(rec, threshold = 0.01, id) {
    recipes_pkg_check(
      required_pkgs_filter_by_variance(),
      "step_filter_by_variance()"
    )
    add_step(
      rec,
      step_filter_by_variance_new(threshold = threshold, id = id)
    )
  }
)

#' @rdname step_filter_by_variance
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_filter_by_variance",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, threshold = 0.01, id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_filter_by_variance_new <- function(threshold = 0.01, id) {
  step(subclass = "filter_by_variance", threshold = threshold, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_filter_by_variance <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_filter_by_variance <- function(rec, threshold = 0.01) {
  rec@phyloseq <- 
    get_phy(rec) %>%
    phyloseq::filter_taxa(function(x) stats::var(x) > threshold, TRUE) 
  
  rec
}
