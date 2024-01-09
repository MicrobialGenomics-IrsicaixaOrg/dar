#' Subset taxa by taxonomic level
#'
#' This is a convenience function around the subset_taxa function from the
#' phyloseq package. It is designed to speed up subsetting complex experimental
#' objects with one function call. In the case of run_subset_taxa,
#' the subsetting will be based on the taxonomic level of each taxon. The taxa
#' retained in the dataset are those where the taxonomic level matches the
#' provided taxa.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param tax_level The taxonomic level for subsetting taxa.
#' @param taxa The taxa to be retained in the dataset.
#' @param id A character string that is unique to this step to identify it.
#'
#' @note This function modifies `rec` in place, you might want to make a copy of
#'   `rec` before modifying it if you need to preserve the original object.
#' @details The function subsets the taxa in the phyloseq object based on the
#'   provided taxonomic level and taxa. Only the taxa that match the provided
#'   taxa at the given taxonomic level are retained in the phyloseq object.
#' @return A Recipe object that has been subsetted based on taxonomic level.
#' @seealso \code{\link[phyloseq]{subset_taxa}}
#' @include recipe-class.R
#' @family subset phy steps
#' @aliases step_subset_taxa
#' @export
#' @autoglobal 
#' @tests
#' data(test_prep_rec)
#' expect_error(step_subset_taxa(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#'
#' ## Define step_subset_taxa step with default parameters
#' rec <- step_subset_taxa(
#'   rec, 
#'   tax_level = "Kingdom",
#'   taxa = c("Bacteria", "Archaea")
#' )
#' rec
methods::setGeneric(
  name = "step_subset_taxa",
  def = function(rec, 
                 tax_level, 
                 taxa, 
                 id = rand_id("subset_taxa")) {
    standardGeneric("step_subset_taxa")
  }
)

#' @rdname step_subset_taxa
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_subset_taxa",
  signature = c(rec = "Recipe"),
  definition = function(rec, tax_level, taxa, id) {
    recipes_pkg_check(
      required_pkgs_subset_taxa(),
      "step_subset_taxa()"
    )
    add_step(
      rec,
      step_subset_taxa_new(
        tax_level = tax_level, 
        taxa = taxa, 
        id = id
      )
    )
  }
)

#' @rdname step_subset_taxa
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_subset_taxa",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec, tax_level, taxa, id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_subset_taxa_new <- function(tax_level, taxa, id) {
  step(
    subclass = "subset_taxa", 
    tax_level = tax_level, 
    taxa = taxa, 
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_subset_taxa <- function(x, ...) {  c("bioc::phyloseq") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_subset_taxa <- function(rec, tax_level, taxa) {
  expr <- paste0(tax_level, ' %in% c(', paste(shQuote(taxa), collapse = ", "), ')')
  rec@phyloseq <-
    glue::glue("phyloseq::subset_taxa(get_phy(rec), {expr})") %>%
    parse(text = .) %>%
    eval()
  
  rec
}
