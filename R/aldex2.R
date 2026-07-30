#' ALDEx2 analysis
#'
#' A differential abundance analysis for the comparison of two or more
#' conditions. For example, single-organism and meta-RNA-seq high-throughput
#' sequencing assays, or of selected and unselected values from in-vitro
#' sequence selections. Uses a Dirichlet-multinomial model to infer abundance
#' from counts, that has been optimized for three or more experimental
#' replicates. Infers sampling variation and calculates the expected false
#' discovery rate given the biological and sampling variation using the Wilcox
#' rank test or Welches t-test (aldex.ttest) or the glm and Kruskal Wallis tests
#' (aldex.glm). Reports both P and fdr values calculated by the Benjamini
#' Hochberg correction (Not supported in `dar` package).
#'
#' @details The run_aldex function is a wrapper that performs log-ratio
#'   transformation and statistical testing in a single line of code.
#'   Specifically, this function: (a) generates Monte Carlo samples of the
#'   Dirichlet distribution for each sample, (b) converts each instance using a
#'   log-ratio transform, then (c) returns test results for two sample (Welch's
#'   t, Wilcoxon) test. This function also estimates effect size for two sample
#'   analyses.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param max_significance Benjamini-Hochberg corrected P value of Welch’s t
#'   test cutoff.
#' @param mc.samples The number of Monte Carlo instances to use to estimate the
#'   underlying distributions; since we are estimating central tendencies, 128
#'   is usually sufficient, but larger numbers may be .
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
#'   by aldex.clr.
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
#' @aliases step_aldex
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
#'  step_aldex() |> 
#'  prep()
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_aldex(test_prep_rec))
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
#' ## Define ALDEX step with default parameters and prep
#' rec <-
#'   step_aldex(rec) |>
#'   prep(parallel = FALSE)
#'
#' rec
#'
#' ## Wearing rarefaction only for this step
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_aldex(rarefy = TRUE)
#'
#' rec
step_aldex <- function(rec,
                       max_significance = 0.05,
                       mc.samples = 128,
                       denom = "all",
                       rarefy = FALSE,
                       id = rand_id("aldex")) {
  
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
      id = id
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
                      id) {

  if (!is.null(get_model(rec))) {
    return(run_aldex_model(
      rec, max_significance, mc.samples, denom, rarefy
    ))
  }

  vars <- get_var(rec)[[1]]
  tax_level <- get_tax(rec)[[1]]
  
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
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
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
