#' ANCOM analysis
#'
#' Determine taxa whose absolute abundances, per unit volume, of the ecosystem
#' (e.g., gut) are significantly different with changes in the covariate of
#' interest (e.g., group). The current version of ancombc2 function implements
#' Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC2) in
#' cross-sectional and repeated measurements data. In addition to the two-group
#' comparison, ANCOM-BC2 also supports testing for continuous covariates and
#' multi-group comparisons, including the global test, pairwise directional
#' test, Dunnett's type of test, and trend test.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param fix_formula the character string expresses how the microbial absolute
#'   abundances for each taxon depend on the fixed effects in metadata. When
#'   specifying the fix_formula, make sure to include the group variable in the
#'   formula if it is not NULL. Defaults to `NULL`: use the centralized model,
#'   or the recipe targets during deprecated model-free execution. An explicit
#'   formula is overridden with a warning when a centralized model is defined.
#' @param rand_formula the character string expresses how the microbial absolute
#'   abundances for each taxon depend on the random effects in metadata.
#'   ANCOM-BC2 follows the lmerTest package in formulating the random effects.
#'   See ?lmerTest::lmer for more details. Default is NULL.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none". See ?stats::p.adjust for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#'   less than prv_cut will be excluded in the analysis. For instance, suppose
#'   there are 100 samples, if a taxon has nonzero counts presented in less than
#'   10 samples, it will not be further analyzed. Default is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#'   sizes. Samples with library sizes less than lib_cut will be excluded in the
#'   analysis. Default is 0, i.e. do not discard any sample.
#' @param s0_perc a numerical fraction between 0 and 1. Inspired by Significance
#'   Analysis of Microarrays (SAM) methodology, a small positive constant is
#'   added to the denominator of ANCOM-BC2 test statistic corresponding to each
#'   taxon to avoid the significance due to extremely small standard errors,
#'   especially for rare taxa. This small positive constant is chosen as
#'   s0_perc-th percentile of standard error values for each fixed effect.
#'   Default is 0.05 (5th percentile).
#' @param group character. The name of the group variable in metadata. group
#'   should be discrete. Specifying group is required for detecting structural
#'   zeros and performing multi-group comparisons (global test, pairwise
#'   directional test, Dunnett's type of test, and trend test). Default is NULL.
#'   If the group of interest contains only two categories, leave it as NULL.
#' @param struc_zero logical. Whether to detect structural zeros based on group.
#'   Default is FALSE. See Details for a more comprehensive discussion on
#'   structural zeros.
#' @param neg_lb logical. Whether to classify a taxon as a structural zero using
#'   its asymptotic lower bound. Default is FALSE.
#' @param alpha numeric. Level of significance. Default is 0.05.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#'   ?parallel::makeCluster. Default is 1 (no parallel computing).
#' @param verbose logical. Whether to generate verbose output during the
#'   ANCOM-BC2 fitting process. Default is FALSE.
#' @param global logical. Whether to perform the global test. Default is FALSE.
#' @param pairwise logical. Whether to perform the pairwise directional test.
#'   Default is FALSE.
#' @param dunnet logical. Whether to perform the Dunnett's type of test. Default
#'   is FALSE.
#' @param trend logical. Whether to perform trend test. Default is FALSE.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native `fit`
#'   stage. Arguments managed by dar or exposed above cannot be overridden.
#' @details Significance follows ANCOM-BC2's sensitivity-aware `diff_robust`
#'   decision when `pseudo_sens = TRUE` (the default), rather than the basic
#'   `diff` decision alone. With a centralized model,
#'   `engine_args = list(fit = list(pseudo_sens = FALSE))` explicitly disables
#'   sensitivity analysis and uses `diff`. Deprecated model-free execution
#'   uses the default sensitivity analysis. Adjusted p-values and effect
#'   estimates are not changed by this selection.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_ancom
#' @return An object of class `Recipe`
#' @autoglobal
#' @export
#' @tests
#' data(metaHIV_phy)
#' 
#' test <-
#'  suppressWarnings(recipe(metaHIV_phy, "RiskGroup2", "Order")) |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_prevalence() |> 
#'  step_ancom(id = "ancom__Utap")
#'  
#' res <- suppressWarnings(prep(test, parallel = FALSE))
#' expect_s4_class(res, "PrepRecipe")
#' expect_snapshot(res)
#' 
#' data(test_prep_rec)
#' expect_error(step_ancom(test_prep_rec, rarefy = TRUE))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus")
#'
#' rec_default <- step_ancom(rec)
#'
#' # Pass advanced convergence controls to ANCOM-BC2.
#' rec_custom <- step_ancom(
#'   rec,
#'   engine_args = list(
#'     fit = list(iter_control = list(tol = 1e-3, max_iter = 50))
#'   )
#' )
#'
#' # Rarefy counts only when this step runs.
#' rec_rarefied <- step_ancom(rec, rarefy = TRUE)
step_ancom <- function(rec,
                       fix_formula = NULL,
                       rand_formula = NULL,
                       p_adj_method = "holm",
                       prv_cut = 0.1,
                       lib_cut = 0,
                       s0_perc = 0.05,
                       group = NULL,
                       struc_zero = FALSE,
                       neg_lb = FALSE,
                       alpha = 0.05,
                       n_cl = 1,
                       verbose = FALSE,
                       global = FALSE,
                       pairwise = FALSE,
                       dunnet = FALSE,
                       trend = FALSE,
                       rarefy = FALSE,
                       id = rand_id("ancom"),
                       engine_args = list()) {

  check_recipe(rec)
  recipes_pkg_check(c("bioc::ANCOMBC"), "step_ancom()")
  add_step(
    rec,
    step(
      subclass = "ancom",
      fix_formula = fix_formula,
      rand_formula = rand_formula,
      p_adj_method = p_adj_method,
      prv_cut = prv_cut,
      lib_cut = lib_cut,
      s0_perc = s0_perc,
      group = group,
      struc_zero = struc_zero,
      neg_lb = neg_lb,
      alpha = alpha,
      n_cl = n_cl,
      verbose = verbose,
      global = global,
      pairwise = pairwise,
      dunnet = dunnet,
      trend = trend,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("ancom", engine_args)
    )
  )
}

#' @noRd
#' @keywords internal
run_ancom <- function(rec,
                      fix_formula,
                      rand_formula,
                      p_adj_method,
                      prv_cut,
                      lib_cut,
                      s0_perc,
                      group,
                      struc_zero,
                      neg_lb,
                      alpha,
                      n_cl,
                      verbose,
                      global,
                      pairwise,
                      dunnet,
                      trend,
                      rarefy,
                      id,
                      engine_args = list()) {

  engine_args <- check_engine_args_execution(
    rec, "ancom", engine_args, id
  )

  if (!is.null(get_model(rec))) {
    return(run_ancom_model(
      rec = rec,
      p_adj_method = p_adj_method,
      prv_cut = prv_cut,
      lib_cut = lib_cut,
      s0_perc = s0_perc,
      struc_zero = struc_zero,
      neg_lb = neg_lb,
      alpha = alpha,
      n_cl = n_cl,
      verbose = verbose,
      rarefy = rarefy,
      engine_args = engine_args
    ))
  }
  
  fix_formula <- fix_formula %||% recipe_targets(rec)
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
      to_remove <- rm_zero_variance(s_phy, var)
      if (length(to_remove) > 0) {
        s_phy <- phyloseq::prune_taxa(
          phyloseq::taxa_names(s_phy)[!phyloseq::taxa_names(s_phy) %in% to_remove],
          s_phy
        )
      }

      res <- ANCOMBC::ancombc2(
        data = mia::convertFromPhyloseq(s_phy),
        tax_level = tax_level,
        fix_formula = fix_formula,
        rand_formula = rand_formula,
        p_adj_method = p_adj_method,
        prv_cut = prv_cut,
        lib_cut = lib_cut,
        s0_perc = s0_perc,
        group = group,
        struc_zero = struc_zero,
        neg_lb = neg_lb,
        alpha = alpha,
        n_cl = n_cl,
        verbose = verbose,
        global = global,
        pairwise = pairwise,
        dunnet = dunnet,
        trend = trend,
        pseudo_sens = TRUE
      )

      ancom_stats_tbl(res$res, var, rec, comparison)
    })
  })
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_ancom <- function(x, ...) { c("bioc::ANCOMBC") }

#' @noRd
#' @keywords internal
#' @autoglobal
ancom_stats_tbl <- function(ancom_res, var, rec, comparison,
                            pseudo_sens = TRUE) {
  result <- ancom_res %>%
    dplyr::select(taxa_id = taxon, dplyr::contains(!!var)) %>%
    dplyr::right_join(analysis_tax_table(rec), ., by = "taxa_id") %>%
    stats::setNames(
      stringr::str_remove_all(names(.), stringr::str_c("_", var, ".*"))
    )
  signif_col <- ancom_significance_column(names(result), pseudo_sens)
  result %>%
    dplyr::mutate(
      comparison = stringr::str_c(comparison, collapse = "_"),
      var = var,
      effect = lfc,
      signif = .data[[signif_col]]
    )
}

#' @noRd
ancom_significance_column <- function(columns, pseudo_sens = TRUE,
                                      coefficient = NULL) {
  prefix <- if (pseudo_sens) "diff_robust" else "diff"
  column <- if (is.null(coefficient)) {
    if (prefix %in% columns) prefix else NA_character_
  } else {
    match_ancom_column(columns, prefix, coefficient)
  }
  if (is.na(column)) {
    cli::cli_abort(
      c(
        "ANCOM-BC2 did not return the expected {.val {prefix}} decision.",
        "i" = "Sensitivity analysis is {.val {pseudo_sens}}; significance cannot be inferred from adjusted p-values alone."
      ),
      class = "dar_error_invalid_result_contract"
    )
  }
  column
}

#' @noRd
#' @keywords internal
#' @autoglobal
rm_zero_variance <- function(phy, vars) {
  phyloseq::sample_data(phy) %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    dplyr::group_split(!!dplyr::sym(vars)) %>% 
    purrr::map_dfr(~ { 
      .x %>% 
        dplyr::pull(sample_id) %>% 
        phyloseq::prune_samples(phy) %>% 
        phyloseq::filter_taxa(function(x) var(x) > 0) %>% 
        tibble::as_tibble(rownames = "otu_id")
    }) %>% 
    dplyr::filter(value == FALSE) %>% 
    dplyr::pull(otu_id) 
}
