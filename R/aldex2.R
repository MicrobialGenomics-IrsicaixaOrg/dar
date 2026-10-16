#' ALDEx2 analysis
#'
#' Differential abundance analysis using ALDEx2's Monte Carlo log-ratio
#' transformation and native statistical tests. Requires [add_model()],
#' fits a GLM for each planned contrast and reports its coefficient and
#' adjusted P value.
#'
#' @details With [add_model()], `effect` is the contrast-oriented coefficient
#'   from `ALDEx2::aldex.glm()`, averaged over Monte Carlo instances, and
#'   `effect_metric` is `"aldex_glm_coefficient"`. It comes from the same GLM
#'   as the reported adjusted P value. Under the default denominator, scale
#'   and Gaussian identity-link settings, the response is CLR in base 2;
#'   this coefficient is not an absolute-abundance log2 fold change or the
#'   standardized effect from `ALDEx2::aldex.effect()`. Advanced native
#'   settings can change the response or link scale.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param max_significance Adjusted P value cutoff for the GLM.
#' @param mc.samples The number of Monte Carlo instances to use to estimate the
#'   underlying distributions. The default is 128.
#' @param denom Features used as the geometric-mean denominator: `"all"`
#'   (default) or explicit numeric row indices in the aggregated count matrix.
#'   Choose indices a priori and account for native removal of all-zero taxa.
#'   Data-dependent denominators are not supported with the shared design.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native `clr`
#'   and `fit` stages. Effects come from the GLM. Arguments managed by dar or
#'   exposed above cannot be overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_aldex2
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
#'  step_filter_by_abundance() |> 
#'  step_aldex2()
#' test <- suppressWarnings(prep(test, parallel = FALSE))
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_aldex2(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_aldex2(rec)
#'
#' # Increase the number of Monte Carlo instances.
#' rec_more_instances <- step_aldex2(rec, mc.samples = 256)
step_aldex2 <- function(rec,
                       max_significance = 0.05,
                       mc.samples = 128,
                       denom = "all",
                       rarefy = FALSE,
                       id = rand_id("aldex2"),
                       engine_args = list()) {
  check_recipe(rec)
  recipes_pkg_check(c("bioc::ALDEx2"), "step_aldex2()")

  add_step(
    rec,
    step(
      subclass = "aldex2",
      max_significance = max_significance,
      mc.samples = mc.samples,
      denom = denom,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("aldex2", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_aldex2 <- function(rec,
                      max_significance,
                      mc.samples,
                      denom,
                      rarefy,
                      id,
                      engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "aldex2", engine_args)
  run_aldex2_model(rec, max_significance, mc.samples, denom, rarefy, engine_args)
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_aldex2 <- function(x, ...) { c("bioc::ALDEx2") }
