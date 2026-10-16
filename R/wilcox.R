#' Wilcox analysis
#'
#' Performs a wilcox test to determine features (be it Operational Taxanomic
#' Unit (OTU), species, etc.) that are differentially abundant between two or
#' more groups of multiple samples.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param norm_method Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#'   'hellinger', 'identity', 'clr', 'alr', or any method from the
#'   vegan::decostand function.
#' @param p_adj_method Character. Specifying the method to adjust p-values for
#'   multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
#' @param max_significance The q-value threshold for significance.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. 
#' @param id Unique, non-empty identifier for this step. DA step IDs must not
#'   use the reserved internal names `taxa_id`, `contrast_id`, `comparison`,
#'   `contrast_type`, `var`, `effect`, `sum`, `.tax_order`, `.contrast_order` or
#'   `.effect_order`. Other identifiers, including names such as `.custom`,
#'   are retained unchanged.
#' @param engine_args Named lists of advanced arguments for the native `test`
#'   stage. Arguments managed by dar or exposed above cannot be overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_wilcox
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' 
#' test <-
#'  recipe(metaHIV_phy) |>
#'  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_variance() |> 
#'  step_wilcox()
#' test <- suppressWarnings(prep(test, parallel = FALSE))
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_wilcox(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_wilcox(rec)
#'
#' # Use Holm's multiple-testing correction.
#' rec_holm <- step_wilcox(rec, p_adj_method = "holm")
step_wilcox <- function(rec,
                        norm_method = "compositional",
                        max_significance = 0.05,
                        p_adj_method = "BH",
                        rarefy = FALSE,
                        id = rand_id("wilcox"),
                        engine_args = list()) {
  
  check_recipe(rec)
  recipes_pkg_check(required_pkgs_wilcox(), "step_wilcox()")
  if (rarefy && !contains_rarefaction(rec)) {
    cli::cli_inform(
      c(
        "!" = "Running {.fun wilcox} with rarefaction is recommended for low depth samples.",
        "i" = "Metadata: {.field id} = {.val {id}}"
      )
    )
  }
  
  add_step(
    rec,
    step(
      subclass = "wilcox",
      norm_method = norm_method,
      max_significance = max_significance,
      p_adj_method = p_adj_method,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("wilcox", engine_args)
    )
  )
}

#' @noRd
#' @keywords internal
required_pkgs_wilcox <- function(x, ...) { "bioc::microbiome"}

#' @noRd
#' @keywords internal
#' @autoglobal
run_wilcox <- function(rec,
                       norm_method,
                       max_significance,
                       p_adj_method,
                       rarefy,
                       id,
                       engine_args = list()) {

  engine_args <- check_engine_args_execution(rec, "wilcox", engine_args)
  run_wilcox_model(rec, norm_method, max_significance, p_adj_method, rarefy, engine_args)
}
