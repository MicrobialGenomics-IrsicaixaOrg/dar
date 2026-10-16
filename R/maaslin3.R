#' MaAsLin3 analysis
#'
#' MaAsLin 3 finds associations between microbiome meta-omics features and
#' complex metadata. It uses a unified framework to test for both abundance
#' (using linear models) and prevalence (using logistic regression).
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param min_abundance The minimum abundance for each feature.
#' @param min_prevalence The minimum percent of samples for which a feature is
#'   detected at minimum abundance.
#' @param min_variance Keep features with variance greater than.
#' @param normalization The normalization method to apply. Default: "TSS".
#'   Choices: "TSS", "CLR", "NONE".
#' @param transform The transform to apply. Default: "LOG" (Base 2). Choices: "LOG",
#'   "PLOG", "NONE".
#' @param max_significance The q-value threshold for significance.
#' @param random_effects The random effects for the model (vector of character strings).
#' @param correction The correction method for computing the q-value.
#' @param standardize Apply z-score so continuous metadata are on the same
#'   scale.
#' @param reference The factor to use as a reference for a variable with more
#'   than two levels provided as a string of 'variable,reference' semi-colon
#'   delimited.
#' @param median_comparison_abundance Test abundance coefficients against a null
#'   value corresponding to the median coefficient for a metadata variable across
#'   the features. Recommended for relative abundance (default: TRUE).
#' @param rarefy Boolean indicating if OTU counts must be rarefyed.
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native `fit` or
#'   `contrast` stage. Arguments managed by dar or exposed above cannot be
#'   overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_maaslin
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#'
#' test <-
#'   recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_by_abundance() |>
#'   step_maaslin() |>
#'   step_maaslin(rarefy = TRUE) |>
#'   step_maaslin(rarefy = "no_seed")
#'
#' expect_s4_class(suppressWarnings(prep(test, parallel = FALSE)), "PrepRecipe") |>
#'   expect_snapshot()
#'
#' data(test_prep_rec)
#' expect_error(step_maaslin(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_maaslin(rec)
#'
#' # Increase the minimum feature prevalence.
#' rec_filtered <- step_maaslin(rec, min_prevalence = 0.2)
step_maaslin <- function(rec,
                         min_abundance = 0.0,
                         min_prevalence = 0.1,
                         min_variance = 0.0,
                         normalization = "TSS",
                         transform = "LOG",
                         max_significance = 0.1,
                         random_effects = NULL,
                         correction = "BH",
                         standardize = TRUE,
                         reference = NULL,
                         median_comparison_abundance = TRUE,
                         rarefy = FALSE,
                         id = rand_id("maaslin"),
                         engine_args = list()) {
  check_recipe(rec)
  recipes_pkg_check(required_pkgs_maaslin(), "step_maaslin()")

  add_step(
    rec,
    step(
      subclass = "maaslin",
      min_abundance = min_abundance,
      min_prevalence = min_prevalence,
      min_variance = min_variance,
      normalization = normalization,
      transform = transform,
      max_significance = max_significance,
      random_effects = random_effects,
      correction = correction,
      standardize = standardize,
      reference = reference,
      median_comparison_abundance = median_comparison_abundance,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("maaslin", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_maaslin <- function(rec,
                        min_abundance,
                        min_prevalence,
                        min_variance,
                        normalization,
                        transform,
                        max_significance,
                        random_effects,
                        correction,
                        standardize,
                        reference,
                        median_comparison_abundance,
                        rarefy,
                        id,
                        engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "maaslin", engine_args)

  run_maaslin_model(
    rec = rec,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    min_variance = min_variance,
    normalization = normalization,
    transform = transform,
    max_significance = max_significance,
    correction = correction,
    standardize = standardize,
    median_comparison_abundance = median_comparison_abundance,
    rarefy = rarefy,
    engine_args = engine_args
  )
}

#' @noRd
#' @keywords internal
required_pkgs_maaslin <- function(x, ...) { c("bioc::maaslin3") }
