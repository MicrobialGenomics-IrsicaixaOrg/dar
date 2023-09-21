#' Subset species by taxonomic expression
#'
#' This is a convenience wrapper around the subset function. It is intended to
#' speed subsetting complex experimental objects with one function call. In the
#' case of subset_taxa, the subsetting will be based on an expression related to
#' the columns and values within the tax_table (taxonomyTable component) slot of
#' physeq. The OTUs retained in the dataset is equivalent to \code{x[subset &
#' !is.na(subset)]}, where x is the vector of OTU IDs and subset is the logical
#' that results from your subsetting expression. This is important to keep in
#' mind, as users are often unaware that this subsetting step also removes/omits
#' OTUs that have a missing value result, NA, somewhere in the expression.
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param expr The subsetting expression that should be applied to the
#'   taxonomyTable. This is passed on to subset, and more details and examples
#'   about how it functions can be found in its documentation.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_subset_taxa
#' @return An object of class `recipe`
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' 
#' ## Init recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#' 
#' ## Define subset taxa step with default parameters
#' rec <- step_subset_taxa(rec, expr = 'Kingdom %in% c("Bacteria", "Archaea")')
#' rec
methods::setGeneric(
  name = "step_subset_taxa",
  def = function(rec, expr, id = rand_id("subset_taxa")) {
    standardGeneric("step_subset_taxa")
  }
)

#' @rdname step_subset_taxa
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_subset_taxa",
  signature = c(rec = "recipe"),
  definition = function(rec, expr, id) {
    recipes_pkg_check(required_pkgs_subset_taxa(), "step_subset_taxa()")
    add_step(rec, step_subset_taxa_new(expr = expr, id = id))
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_subset_taxa_new <- function(expr, id) {
  step(subclass = "subset_taxa", expr = expr, id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_subset_taxa <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_subset_taxa <- function(rec, expr) {
  rec@phyloseq <-
    glue::glue("phyloseq::subset_taxa(get_phy(rec), {expr})") %>%
    parse(text = .) %>%
    eval()

  rec
}
