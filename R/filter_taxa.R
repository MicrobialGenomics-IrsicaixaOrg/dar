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
#' data(test_prep_rec)
#' expect_error(step_filter_taxa(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#' 
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy)
#' rec
#' 
#' ## Define filter taxa step with default parameters
#' rec <- 
#'   step_filter_taxa(rec, .f = "function(x) sum(x > 0) >= (0.03 * length(x))")
#'   
#' rec
step_filter_taxa <- function(rec, .f, id = rand_id("filter_taxa")) {
  
  check_recipe(rec)
  recipes_pkg_check(
    required_pkgs_filter_taxa(), 
    "step_filter_taxa()"
  )
  
  add_step(
    rec,
    step(
      subclass = "filter_taxa", 
      .f = .f, 
      id = id
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_filter_taxa <- function(rec, .f, id) {
  if (is.character(.f)) { .f <- eval(parse(text = .f)) }
  rec@phyloseq <- phyloseq::filter_taxa(get_phy(rec), .f, prune = TRUE)
 
  rec
}

#' @noRd
#' @keywords internal
required_pkgs_filter_taxa <- function(x, ...) {  c("bioc::phyloseq") }
