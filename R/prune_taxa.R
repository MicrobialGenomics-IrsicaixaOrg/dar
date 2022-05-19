#' Filter taxa based on across-sample OTU abundance criteria
#'
#' This function is directly analogous to the genefilter function for microarray
#' filtering, but is used for filtering OTUs from phyloseq objects. It applies an
#' arbitrary set of functions — as a function list, for instance, created by filterfun —
#' as across-sample criteria, one OTU at a time. It takes as input a phyloseq object, and
#' returns a logical vector indicating whether or not each OTU passed the criteria.
#' Alternatively, if the "prune" option is set to FALSE, it returns the already-trimmed
#' version of the phyloseq object.
#'
#' @param rec A recipe object. The step will be added to the sequence of operations for
#'   this recipe.
#' @param .f A function or list of functions that take a vector of abundance values and
#'   return a logical. Some canned useful function types are included in the
#'   genefilter-package.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_taxa
#' @return An object of class `recipe`
#' @export
methods::setGeneric(
  name = "step_filter_taxa",
  def = function(rec, .f, id = rand_id("filter_taxa")) {
    standardGeneric("step_filter_taxa")
  }
)

#' @rdname step_filter_taxa
#' @export
methods::setMethod(
  f = "step_filter_taxa",
  signature = c(rec = "recipe"),
  definition = function(rec, .f, id) {
    recipes_pkg_check(required_pkgs_filter_taxa())
    add_step(rec, step_filter_taxa_new(.f = .f, id = id))
  }
)

#' @rdname step_filter_taxa
#' @keywords internal
step_filter_taxa_new <- function(.f, id) {
  step(subclass = "filter_taxa", .f = .f, id = id)
}

#' @noRd
#' @keywords internal
required_pkgs_filter_taxa <- function(x, ...) {  c("bioc::phyloseq") }

#' @rdname step_filter_taxa
#' @export
run_filter_taxa <- function(rec, .f) {
  rec@phyloseq <- phyloseq::filter_taxa(get_phy(rec), eval(parse(text = .f)), prune = TRUE)
  rec
}
