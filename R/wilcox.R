#' Wilcox analysis
#'
#' Performs a wilcox test to determine features (be it Operational Taxanomic
#' Unit (OTU), species, etc.) that are differentially abundant between two or
#' more groups of multiple samples.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param norm_method Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#'   'hellinger', 'identity', 'clr', 'alr', or any method from the
#'   vegan::decostand function.
#' @param p_adj_method Character. Specifying the method to adjust p-values for
#'   multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
#' @param max_significance The q-value threshold for significance.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. 
#' @param id A character string that is unique to this step to identify it.
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
#'  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_variance() |> 
#'  step_wilcox() |> 
#'  prep()
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_wilcox(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea"))
#'
#' ## Define step with default parameters
#' rec <- step_wilcox(rec)
#' rec
step_wilcox <- function(rec,
                        norm_method = "compositional",
                        max_significance = 0.05,
                        p_adj_method = "BH",
                        rarefy = FALSE,
                        id = rand_id("wilcox")) {
  
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
      id = id
    )
  )
}

#' @noRd
#' @keywords internal
required_pkgs_wilcox <- function(x, ...) {
  c("bioc::microbiome", "rstatix")
}

#' @noRd
#' @keywords internal
#' @autoglobal
run_wilcox <- function(rec,
                       norm_method,
                       max_significance,
                       p_adj_method,
                       rarefy,
                       id) {

  if (!is.null(get_model(rec))) {
    return(run_wilcox_model(
      rec, norm_method, max_significance, p_adj_method, rarefy
    ))
  }

  phy <- get_phy(rec) %>% use_rarefy(rarefy)
  tax_level <- get_tax(rec)

  prepro_df <- phy %>%
    phyloseq::tax_glom(taxrank = tax_level, NArm = FALSE) %>%
    microbiome::transform(transform = norm_method) %>%
    phyloseq::otu_table() %>%
    t() %>%
    to_tibble("sample_id") %>%
    dplyr::left_join(sample_data(rec), ., by = "sample_id")

  vars <- get_var(rec)
  
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      meta_cols <- c("sample_id", var, names(sample_data(rec)))
      features <- setdiff(names(prepro_df), meta_cols)
      
      features %>%
        purrr::map_dfr(function(feature) {
          prep <- prepro_df %>% dplyr::select(sample_id, !!var, !!feature)
          to_exclude <- prep %>%
            dplyr::group_by(!!dplyr::sym(var)) %>%
            dplyr::summarise(sum = sum(!!dplyr::sym(feature), na.rm = TRUE), .groups = "drop") %>%
            dplyr::filter(sum == 0) %>%
            dplyr::pull(!!var)

          comparisons <- get_comparisons(var, get_phy(rec), as_list = TRUE, n_cut = 1) %>% 
            purrr::map(~ c(.x[[2]], .x[[1]]))
          
          if (length(to_exclude) >= 2) {
            comparisons <- comparisons %>% 
              purrr::discard(~ all(.x %in% to_exclude))
          }
          
          if (length(comparisons) == 0) return(NULL)

          res <- rstatix::wilcox_test(
            data = prep,
            formula = stats::as.formula(paste0("`", feature, "` ~ `", var, "`")),
            comparisons = comparisons,
            p.adjust.method = p_adj_method,
            detailed = TRUE
          ) %>%
            tidyr::unite("comparison", group1:group2, sep = "_") %>% 
            dplyr::rename(taxa_id = .y.)
          
          p_col <- if ("p.adj" %in% names(res)) "p.adj" else "p"
          res <- 
            res %>% 
            dplyr::rename(padj = !!p_col) %>% 
            dplyr::left_join(tax_table(rec), by = "taxa_id")
        }) %>% 
        dplyr::mutate(
          estimate = -estimate, 
          effect = estimate,
          var = !!var,
          signif = padj < max_significance
        ) %>% 
        dplyr::relocate(taxa_id, .before = 1)
    })
}
