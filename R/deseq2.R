#' DESeq2 analysis
#'
#' Differential expression analysis based on the Negative Binomial (a.k.a.
#' Gamma-Poisson) distribution. This function performs a default analysis
#' through the steps: 1) estimation of size factors: `estimateSizeFactors`. 2)
#' estimation of dispersion: `estimateDispersions`. 3) Negative Binomial GLM
#' fitting and Wald statistics: nbinomWaldTest. For complete details on each
#' step, see the manual pages of the respective functions. After the DESeq
#' function returns a DESeqDataSet object, results tables (log2 fold changes and
#' p-values) can be generated using the results function. Shrunken LFC can then
#' be generated using the lfcShrink function.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param test Either "Wald" or "LRT", which will then use either Wald
#'   significance tests (defined by nbinomWaldTest), or the likelihood ratio
#'   test on the difference in deviance between a full and reduced model formula
#'   (defined by nbinomLRT).
#' @param fitType either "parametric", "local", "mean", or "glmGamPoi" for the
#'   type of fitting of dispersions to the mean intensity. See
#'   estimateDispersions for description.
#' @param betaPrior whether or not to put a zero-mean normal prior on the
#'   non-intercept coefficients See nbinomWaldTest for description of the
#'   calculation of the beta prior. In versions >=1.16, the default is set to
#'   FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
#' @param type "apeglm" is the adaptive Student's t prior shrinkage estimator
#'   from the 'apeglm' package; "ashr" is the adaptive shrinkage estimator from
#'   the 'ashr' package, using a fitted mixture of normals prior - see the
#'   Stephens (2016) reference below for citation; "normal" is the 2014 DESeq2
#'   shrinkage estimator using a Normal prior.
#' @param max_significance The q-value threshold for significance.
#' @param log2FC log2FC cutoff.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native
#'   `size_factors`, `fit`, `results`, or `shrink` stage. Arguments managed by
#'   dar or exposed above cannot be overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_deseq
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
#'  step_filter_by_rarity(0.95) |> 
#'  step_deseq()
#' test <- suppressWarnings(prep(test, parallel = FALSE))
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_deseq(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_deseq(rec)
#'
#' # Pass advanced arguments to native DESeq2 stages.
#' rec_custom <- step_deseq(
#'   rec,
#'   engine_args = list(
#'     fit = list(minReplicatesForReplace = Inf),
#'     shrink = list(lfcThreshold = 1)
#'   )
#' )
step_deseq <- function(rec,
                       test = "Wald",
                       fitType = "local",
                       betaPrior = FALSE,
                       type = "ashr",
                       max_significance = 0.05,
                       log2FC = 0,
                       rarefy = FALSE,
                       id = rand_id("deseq"),
                       engine_args = list()) {
  check_recipe(rec)
  if (type == "ashr") {
    recipes_pkg_check(required_pkgs_deseq()[-3], "step_deseq()")
  } else {
    recipes_pkg_check(required_pkgs_deseq()[-2], "step_deseq()")
  }

  add_step(
    rec,
    step(
      subclass = "deseq",
      test = test,
      fitType = fitType,
      betaPrior = betaPrior,
      type = type,
      max_significance = max_significance,
      log2FC = log2FC,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("deseq", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_deseq <- function(rec,
                      test,
                      fitType,
                      betaPrior,
                      type,
                      max_significance,
                      log2FC,
                      rarefy,
                      id,
                      engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "deseq", engine_args)

  run_deseq_model(
    rec,
    test,
    fitType,
    betaPrior,
    type,
    max_significance,
    log2FC,
    rarefy,
    engine_args
  )
}

#' @noRd
#' @keywords internal
required_pkgs_deseq <- function(x, ...) { c("bioc::DESeq2", "bioc::apeglm", "ashr") }
