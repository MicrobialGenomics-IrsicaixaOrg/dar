#' lefse analysis
#'
#' Lefser is metagenomic biomarker discovery tool that is based on LEfSe tool
#' and is published by Huttenhower et al. 2011. Lefser is the R implementation
#' of the LEfSe method. Using statistical analyses, lefser compares microbial
#' populations of healthy and diseased subjects to discover differencially
#' expressed microorganisms. Lefser than computes effect size, which estimates
#' magnitude of differential expression between the populations for each
#' differentially expressed microorganism. Subclasses of classes can also be
#' assigned and used within the analysis.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param kruskal.threshold numeric(1) The p-value for the Kruskal-Wallis Rank
#'   Sum Test (default 0.05).
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#'   when 'blockCol' is present (default 0.05).
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param subclassCol character(1) Optional column name in 'colData(expr)'
#'   indicating the blocks, usually a factor with two levels (e.g., 'c("adult",
#'   "senior")'; default NULL).
#' @param assay The i-th assay matrix in the ‘SummarizedExperiment' (’expr';
#'   default 1).
#' @param trim.names If 'TRUE' extracts the most specific taxonomic rank of
#'   organism.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed.
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native
#'   `transform` or `fit` stage. Arguments managed by dar or exposed above
#'   cannot be overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_lefse
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#'
#' expect_condition(
#'  test <-
#'    recipe(metaHIV_phy) |>
#'    add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
#'    step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'    step_filter_by_prevalence() |>
#'    step_lefse() |>
#'    step_lefse(rarefy = FALSE),
#'  "lefse.*without rarefaction"
#' )
#'
#' expect_s4_class(suppressWarnings(prep(test, parallel = FALSE)), "PrepRecipe") |>
#'   expect_snapshot()
#'
#' data(test_prep_rec)
#' expect_error(step_lefse(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_lefse(rec)
#'
#' # Require a larger LDA score.
#' rec_higher_lda <- step_lefse(rec, lda.threshold = 3)
step_lefse <- function(rec,
                       kruskal.threshold = 0.05,
                       wilcox.threshold = 0.05,
                       lda.threshold = 2,
                       subclassCol = NULL,
                       assay = 1L,
                       trim.names = FALSE,
                       rarefy = TRUE,
                       id = rand_id("lefse"),
                       engine_args = list()) {
  check_recipe(rec)
  recipes_pkg_check(required_pkgs_lefse(), "step_lefse()")

  if (!rarefy & !contains_rarefaction(rec)) {
    cli::cli_inform(
      c(
        "!" = "Running {.pkg lefse} without rarefaction is not recommended ({.arg id} = {.val {id}})."
      )
    )
  }

  add_step(
    rec,
    step(
      subclass = "lefse",
      kruskal.threshold = kruskal.threshold,
      wilcox.threshold = wilcox.threshold,
      lda.threshold = lda.threshold,
      subclassCol = subclassCol,
      assay = assay,
      trim.names = trim.names,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("lefse", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_lefse <- function(rec,
                      kruskal.threshold,
                      wilcox.threshold,
                      lda.threshold,
                      subclassCol,
                      assay,
                      trim.names,
                      rarefy,
                      id,
                      engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "lefse", engine_args)

  run_lefse_model(
    rec,
    kruskal.threshold,
    wilcox.threshold,
    lda.threshold,
    assay,
    trim.names,
    rarefy,
    engine_args
  )
}

#' @noRd
#' @keywords internal
required_pkgs_lefse <- function(x, ...) {
  c("bioc::lefser", "bioc::SummarizedExperiment")
}

#' @noRd
#' @autoglobal
#' @keywords internal
kruskal_test <- function(se, levels, assay = 1L) {
  group <- levels %>% as.factor() %>% as.numeric()
  expr <- SummarizedExperiment::assay(se, i = assay)
  apply(expr, 1L, function(x) { stats::kruskal.test(x ~ group)[["p.value"]] }) %>%
    tibble::as_tibble(rownames = "otu") %>%
    dplyr::rename(pvalue = value) %>%
    dplyr::mutate(adjp = stats::p.adjust(pvalue, method = "BH"))
}
