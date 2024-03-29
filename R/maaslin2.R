#' MaAsLin2 analysis
#'
#' MaAsLin2 finds associations between microbiome meta-omics features and
#' complex metadata in population-scale epidemiological studies. The software
#' includes multiple analysis methods (including support for multiple covariates
#' and repeated measures), filtering, normalization, and transform options to
#' customize analysis for your specific study.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param min_abundance The minimum abundance for each feature.
#' @param min_prevalence The minimum percent of samples for which a feature is
#'   detected at minimum abundance.
#' @param min_variance Keep features with variance greater than.
#' @param normalization The normalization method to apply. Default: "TSS".
#'   Choices: "TSS", "CLR", "CSS", "NONE", "TMM".
#' @param transform The transform to apply. Default: "LOG". Choices: "LOG",
#'   "LOGIT", "AST", "NONE".
#' @param analysis_method The analysis method to apply. Default: "LM". Choices:
#'   "LM", "CPLM", "ZICP", "NEGBIN", "ZINB".
#' @param max_significance The q-value threshold for significance.
#' @param random_effects The random effects for the model, comma-delimited for
#'   multiple effects.
#' @param correction The correction method for computing the q-value.
#' @param standardize Apply z-score so continuous metadata are on the same
#'   scale.
#' @param reference The factor to use as a reference for a variable with more
#'   than two levels provided as a string of 'variable,reference' semi-colon
#'   delimited for multiple variables.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
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
                 analysis_method = "LM",
                 max_significance = 0.25,
                 random_effects = NULL,
                 correction = "BH",
                 standardize = TRUE,
                 reference = NULL,
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
                        analysis_method,
                        max_significance,
                        random_effects,
                        correction,
                        standardize,
                        reference,
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
        analysis_method = analysis_method,
        max_significance = max_significance,
        random_effects = random_effects,
        correction = correction,
        standardize = standardize,
        reference = reference,
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
                        analysis_method,
                        max_significance,
                        random_effects,
                        correction,
                        standardize,
                        reference,
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
                             analysis_method,
                             max_significance,
                             random_effects,
                             correction,
                             standardize,
                             reference,
                             rarefy,
                             id) {
  step(
    subclass = "maaslin",
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    min_variance = min_variance,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    max_significance = max_significance,
    random_effects = random_effects,
    correction = correction,
    standardize = standardize,
    reference = reference,
    rarefy = rarefy,
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_maaslin <- function(x, ...) { c("bioc::Maaslin2") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_maaslin <- function(rec,
                        min_abundance,
                        min_prevalence,
                        min_variance,
                        normalization,
                        transform,
                        analysis_method,
                        max_significance,
                        random_effects,
                        correction,
                        standardize,
                        reference,
                        rarefy) {

  output <- glue::glue("{tempdir()}/maaslin_output")

  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  phy <- 
    get_phy(rec) %>% 
    use_rarefy(rarefy)
  
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

          maaslin2_quietly(
            input_data,
            input_metadata,
            output,
            var,
            rec,
            min_abundance,
            min_prevalence,
            min_variance,
            normalization,
            transform,
            analysis_method,
            max_significance,
            random_effects,
            fixed_effects = var,
            correction,
            standardize,
            reference
          ) %>%
            purrr::pluck("results") %>%
            tibble::as_tibble() %>%
            dplyr::select(-metadata, -value) %>%
            dplyr::mutate(
              coef = coef / log10(2),
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = !!var
            ) %>%
            dplyr::rename(taxa_id = feature) %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>% 
            dplyr::mutate(
              effect = coef,
              signif = ifelse(qval < max_significance, TRUE, FALSE)
            )
        })
    })
}

#' Maaslin2 without log infromation
#'
#' @noRd
#' @keywords internal
#' @autoglobal
maaslin2_quietly <- function(input_data,
                             input_metadata,
                             output,
                             var,
                             rec,
                             min_abundance,
                             min_prevalence,
                             min_variance,
                             normalization,
                             transform,
                             analysis_method,
                             max_significance,
                             random_effects,
                             fixed_effects,
                             correction,
                             standardize,
                             reference,
                             rarefy) {
  f_quietly <- function(input_data,
                        input_metadata,
                        output,
                        var,
                        rec,
                        min_abundance,
                        min_prevalence,
                        min_variance,
                        normalization,
                        transform,
                        analysis_method,
                        max_significance,
                        random_effects,
                        fixed_effects,
                        correction,
                        standardize,
                        reference) {
    res <- Maaslin2::Maaslin2(
      input_data = input_data,
      input_metadata = input_metadata,
      output = output,
      fixed_effects = fixed_effects,
      min_abundance = min_abundance,
      min_prevalence = min_prevalence,
      min_variance = min_variance,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      max_significance = max_significance,
      random_effects = random_effects,
      correction = correction,
      standardize = standardize,
      reference = reference,
      plot_heatmap = FALSE,
      plot_scatter = FALSE
    )
  }
  purrr::quietly(f_quietly)(
    input_data,
    input_metadata,
    output,
    var,
    rec,
    min_abundance,
    min_prevalence,
    min_variance,
    normalization,
    transform,
    analysis_method,
    max_significance,
    random_effects,
    fixed_effects,
    correction,
    standardize,
    reference
  )$result
}
