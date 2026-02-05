#' MaAsLin3 analysis
#'
#' MaAsLin 3 finds associations between microbiome meta-omics features and
#' complex metadata. It uses a unified framework to test for both abundance
#' (using linear models) and prevalence (using logistic regression).
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
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
#' @param id A character string that is unique to this step to identify it.
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
#'  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_abundance() |>
#'  step_maaslin() |>
#'  step_maaslin(rarefy = TRUE) |>
#'  step_maaslin(rarefy = "no_seed")
#'
#' expect_s4_class(prep(test), "PrepRecipe") |>
#'   expect_snapshot()
#'
#' data(test_prep_rec)
#' expect_error(step_maaslin(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.4 * length(x))")
#'
#' rec
#'
#' ## Define step with default parameters and prep
#' rec <-
#'   step_maaslin(rec) |>
#'   prep(parallel = FALSE)
#'
#' rec
#'
#' ## Wearing rarefaction only for this step
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_maaslin(rarefy = TRUE)
#'
#' rec
methods::setGeneric(
  name = "step_maaslin",
  def = function(rec,
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
                 id = rand_id("maaslin")) {
    standardGeneric("step_maaslin")
  }
)

#' @rdname step_maaslin
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_maaslin",
  signature = c(rec = "Recipe"),
  definition = function(rec,
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
                        id) {
    
    recipes_pkg_check(required_pkgs_maaslin(), "step_maaslin()")
    add_step(
      rec,
      step_maaslin_new(
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
        id = id
      )
    )
  }
)

#' @rdname step_maaslin
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_maaslin",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec,
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
                        id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_maaslin_new <- function(rec,
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
                             id) {
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
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_maaslin <- function(x, ...) { c("bioc::maaslin3") }

#' @noRd
#' @keywords internal
#' @autoglobal
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
                        rarefy) {

  output <- glue::glue("{tempdir()}/maaslin3_output")

  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  phy <- get_phy(rec) %>% use_rarefy(rarefy)
  phy <- phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
        purrr::map_dfr(function(comparison) {

          input_metadata <-
            sample_data(rec) %>%
            dplyr::filter(!!dplyr::sym(var) %in% comparison) %>%
            data.frame(row.names = 1)

          input_data <-
            otu_table(rec) %>%
            dplyr::select(taxa_id, dplyr::all_of(rownames(input_metadata))) %>%
            data.frame(row.names = 1) %>%
            as.matrix()

          m3_res <- maaslin3_quietly(
            input_data,
            input_metadata,
            output,
            var,
            min_abundance,
            min_prevalence,
            min_variance,
            normalization,
            transform,
            max_significance,
            random_effects,
            fixed_effects = var,
            correction,
            standardize,
            reference,
            median_comparison_abundance
          ) 
          
          m3_res %>%
            purrr::pluck("fit_data_abundance", "results") %>%
            tibble::as_tibble() %>%
            dplyr::select(
              taxa_id = feature,
              coef,
              stderr,
              pval = pval_individual,
              qval = qval_individual
            ) %>%
            dplyr::mutate(
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = !!var
            ) %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::mutate(
              effect = coef,
              signif = ifelse(qval < max_significance, TRUE, FALSE)
            )
        })
    })
}

#' Maaslin3 without console noise
#'
#' @noRd
#' @keywords internal
#' @autoglobal
maaslin3_quietly <- function(
  input_data,
  input_metadata,
  output,
  var,
  min_abundance,
  min_prevalence,
  min_variance,
  normalization,
  transform,
  max_significance,
  random_effects,
  fixed_effects,
  correction,
  standardize,
  reference,
  median_comparison_abundance
) {
  f_quietly <- function(...) {
    maaslin3::maaslin3(
      input_data = input_data,
      input_metadata = input_metadata,
      output = output,
      fixed_effects = fixed_effects,
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
      verbosity = 'ERROR'
    )
  }

  purrr::quietly(f_quietly)()$result
}