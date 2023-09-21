#' Resample an OTU table such that all samples have the same library size.
#'
#' Please note that the authors of phyloseq do not advocate using this as a
#' normalization procedure, despite its recent popularity. Our justifications
#' for using alternative approaches to address disparities in library sizes have
#' been made available as an article in PLoS Computational Biology. See
#' phyloseq_to_deseq2 for a recommended alternative to rarefying directly
#' supported in the phyloseq package, as well as the supplemental materials for
#' the PLoS-CB article and the phyloseq extensions repository on GitHub.
#' Nevertheless, for comparison and demonstration, the rarefying procedure is
#' implemented here in good faith and with options we hope are useful. This
#' function uses the standard R sample function to resample from the abundance
#' values in the otu_table component of the first argument, physeq. Often one of
#' the major goals of this procedure is to achieve parity in total number of
#' counts between samples, as an alternative to other formal normalization
#' procedures, which is why a single value for the sample.size is expected. This
#' kind of resampling can be performed with and without replacement, with
#' replacement being the more computationally-efficient, default setting. See
#' the replace parameter documentation for more details. We recommended that you
#' explicitly select a random number generator seed before invoking this
#' function, or, alternatively, that you explicitly provide a single positive
#' integer argument as rngseed.
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family rarefaction phy steps
#' @aliases step_rarefaction
#' @return An object of class `recipe`
#' @export
#' @autoglobal
#' @examples 
#' data(metaHIV_phy)
#' 
#' ## Init recipe
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species")
#' 
#' rec
#' 
#' ## Define step with default parameters and prep
#' rec <- 
#'   step_rarefaction(rec) %>% 
#'   prep(parallel = TRUE)
#'   
#' rec
methods::setGeneric(
  name = "step_rarefaction",
  def = function(rec, id = rand_id("rarefaction")) {
    standardGeneric("step_rarefaction")
  }
)

#' @rdname step_rarefaction
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_rarefaction",
  signature = c(rec = "recipe"),
  definition = function(rec, id) {
    recipes_pkg_check(required_pkgs_rarefaction(), "step_rarefaction()")
    add_step(rec, step_rarefaction_new(id = id))
  }
)

#' @rdname step_rarefaction
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_rarefaction",
  signature = c(rec = "prep_recipe"),
  definition = function(rec, id) {
    rlang::abort("This function needs a non-prep recipe!")
  }
)


#' @noRd
#' @keywords internal
#' @autoglobal
step_rarefaction_new <- function(id) {
  step(subclass = "rarefaction", id = id)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_rarefaction <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_rarefaction <- function(rec) {
  rec@phyloseq <- 
    phyloseq::rarefy_even_depth(get_phy(rec), rngseed = 1234, verbose = FALSE)
  rec
}
