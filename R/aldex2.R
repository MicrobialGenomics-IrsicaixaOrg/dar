#' ALDEx2 analysis
#'
#' Differential abundance analysis using ALDEx2's Monte Carlo log-ratio
#' transformation and native statistical tests. With a centralized model,
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
#'   settings can change the response or link scale. Recompute previously
#'   prepared centralized ALDEx2 results to obtain these coefficients.
#'
#'   Deprecated model-free execution retains native two-group Welch tests
#'   and standardized effects from `ALDEx2::aldex.effect()`, identified as
#'   `"standardized_clr_effect"` in tidy results.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param max_significance Adjusted P value cutoff for the GLM, or the
#'   Benjamini-Hochberg corrected Welch P value in model-free execution.
#' @param mc.samples The number of Monte Carlo instances to use to estimate the
#'   underlying distributions. The default is 128.
#' @param denom An any variable (all, iqlr, zero, lvha, median, user) indicating
#'   features to use as the denominator for the Geometric Mean calculation The
#'   default "all" uses the geometric mean abundance of all features. Using
#'   "median" returns the median abundance of all features. Using "iqlr" uses
#'   the features that are between the first and third quartile of the variance
#'   of the clr values across all samples. Using "zero" uses the non-zero
#'   features in each grop as the denominator. This approach is an extreme case
#'   where there are many nonzero features in one condition but many zeros in
#'   another. Using "lvha" uses features that have low variance (bottom
#'   quartile) and high relative abundance (top quartile in every sample). It is
#'   also possible to supply a vector of row indices to use as the denominator.
#'   Here, the experimentalist is determining a-priori which rows are thought to
#'   be invariant. In the case of RNA-seq, this could include ribosomal protein
#'   genes and and other house-keeping genes. This should be used with caution
#'   because the offsets may be different in the original data and in the data
#'   used by the function because features that are 0 in all samples are removed
#'   by aldex.clr. Centralized models require `"all"` or explicit numeric
#'   denominator indices.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native `clr`
#'   and `fit` stages. The former `effect` stage is not supported: centralized
#'   effects now come from the GLM. Arguments managed by dar or exposed above
#'   cannot be overridden.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_aldex
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' 
#' test <-
#'  suppressWarnings(recipe(metaHIV_phy, "RiskGroup2", "Phylum")) |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_abundance() |> 
#'  step_aldex()
#' test <- suppressWarnings(prep(test, parallel = FALSE))
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_aldex(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_aldex(rec)
#'
#' # Increase the number of Monte Carlo instances.
#' rec_more_instances <- step_aldex(rec, mc.samples = 256)
step_aldex <- function(rec,
                       max_significance = 0.05,
                       mc.samples = 128,
                       denom = "all",
                       rarefy = FALSE,
                       id = rand_id("aldex"),
                       engine_args = list()) {
  
  check_recipe(rec)
  recipes_pkg_check(c("bioc::ALDEx2"), "step_aldex()")
  
  add_step(
    rec,
    step(
      subclass = "aldex",
      max_significance = max_significance,
      mc.samples = mc.samples,
      denom = denom,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("aldex", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_aldex <- function(rec,
                      max_significance,
                      mc.samples,
                      denom,
                      rarefy,
                      id,
                      engine_args = list()) {

  engine_args <- check_engine_args_execution(rec, "aldex", engine_args, id)

  if (!is.null(get_model(rec))) {
    return(run_aldex_model(
      rec, max_significance, mc.samples, denom, rarefy, engine_args
    ))
  }

  vars <- recipe_targets(rec)
  tax_level <- recipe_tax_level(rec)
  
  phy <- 
    get_phy(rec) %>% 
    use_rarefy(rarefy) %>%
    phyloseq::tax_glom(taxrank = tax_level, NArm = FALSE)
  
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
        purrr::map_dfr(function(comparison) {
          meta_vals <- phyloseq::get_variable(phy, var)
          s_phy <- phyloseq::prune_samples(meta_vals %in% comparison, phy)
          otu_table <- phyloseq::otu_table(s_phy)
          conds <- phyloseq::get_variable(s_phy, var) %>% as.character()

          clr <- suppressMessages(
            ALDEx2::aldex.clr(
              reads = otu_table,
              conds = conds,
              denom = denom,
              mc.samples = mc.samples
            )
          )

          dplyr::bind_cols(
            ALDEx2::aldex.ttest(clr, verbose = FALSE),
            ALDEx2::aldex.effect(clr, CI = TRUE, verbose = FALSE)
          ) %>%
            dplyr::mutate(
              comparison = stringr::str_c(comparison, collapse = "_"), 
              var = var
            ) %>%
            tibble::as_tibble(rownames = "taxa_id") %>%
            dplyr::left_join(analysis_tax_table(rec), by = "taxa_id") %>%
            dplyr::mutate(
              effect = effect, 
              signif = ifelse(we.eBH < max_significance, TRUE, FALSE),
            )
        })
    })
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_aldex <- function(x, ...) { c("bioc::ALDEx2") }
