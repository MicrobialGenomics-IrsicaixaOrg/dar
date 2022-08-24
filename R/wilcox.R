#' Wilcox analysis
#'
#'
#' Performs a wilcox test to determine features (be it Operational Taxanomic
#' Unit (OTU), species, etc.) that are differentially abundant between two or
#' more groups of multiple samples.
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param norm_method Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#'   'hellinger', 'identity', 'clr', 'alr', or any method from the
#'   vegan::decostand function.
#' @param p_adj_method Character. Specifying the method to adjust p-values for
#'   multiple comparisons. Default is “BH” (Benjamini-Hochberg procedure).
#' @param max_significance The q-value threshold for significance.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_wilcox
#' @return An object of class `recipe`
#' @export
methods::setGeneric(
  name = "step_wilcox",
  def = function(rec,
                 norm_method = "compositional",
                 max_significance = 0.05,
                 p_adj_method = "BH",
                 rarefy = TRUE,
                 id = rand_id("wilcox")) {
    standardGeneric("step_wilcox")
  }
)

#' @rdname step_wilcox
#' @export
methods::setMethod(
  f = "step_wilcox",
  signature = c(rec = "recipe"),
  definition = function(rec,
                        norm_method,
                        max_significance,
                        p_adj_method,
                        rarefy,
                        id) {

    recipes_pkg_check(required_pkgs_wilcox(), "step_wilcox()")
    if (!rarefy & !contains_rarefaction(rec)) {
      rlang::inform(c(
        "!" = glue::glue(
          "Run wilcox without rarefaction is not recommended ", 
          "({crayon::blue(paste0('id = ', id))})"
        )
      ))
    }

    add_step(
      rec,
      step_wilcox_new(
        norm_method = norm_method,
        max_significance = max_significance,
        p_adj_method = p_adj_method,
        rarefy, rarefy,
        id = id
      )
    )
  }
)

#' @rdname step_wilcox
#' @export
methods::setMethod(
  f = "step_wilcox",
  signature = c(rec = "prep_recipe"),
  definition = function(rec,
                        norm_method,
                        max_significance,
                        p_adj_method,
                        rarefy,
                        id) {
    rlang::abort("This function needs a non-prep recipe!")
  }
)

#' @noRd
#' @keywords internal
step_wilcox_new <- function(rec,
                            norm_method,
                            max_significance,
                            p_adj_method,
                            rarefy,
                            id) {
  step(
    subclass = "wilcox",
    norm_method = norm_method,
    max_significance = max_significance,
    p_adj_method = p_adj_method,
    rarefy = rarefy,
    id = id
  )
}

#' @noRd
#' @keywords internal
required_pkgs_wilcox <-
  function(x, ...) {
    c("bioc::microbiome", "rstatix")
  }


#' @noRd
#' @keywords internal
run_wilcox <- function(rec,
                       norm_method,
                       max_significance,
                       p_adj_method,
                       rarefy) {

  phy <- get_phy(rec)
  if (rarefy) { 
    phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE) 
  }

  prepro_df <-
    phy %>%
    phyloseq::tax_glom(taxrank = get_tax(rec), NArm = FALSE) %>%
    microbiome::transform(transform = norm_method) %>%
    phyloseq::otu_table() %>%
    t() %>%
    to_tibble("sample_id") %>%
    dplyr::left_join(sample_data(rec), ., by = "sample_id")

  get_var(rec) %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      dplyr::select(prepro_df, -sample_id, -!!var) %>%
        names() %>%
        purrr::map_dfr( ~ {
          prep <-
            prepro_df %>%
            dplyr::select(sample_id, !!var, !!.x)

          to_exclude <-
            dplyr::group_by(prep, !!dplyr::sym(var)) %>%
            dplyr::summarise(sum = sum(!!dplyr::sym(.x))) %>%
            dplyr::filter(sum == 0) %>%
            dplyr::pull(!!var)

          comparisons <- 
            get_comparisons(var, get_phy(rec), as_list = TRUE, n_cut = 1)
          
          if (length(to_exclude) >= 2) {
            comparisons <-
              comparisons %>% purrr::discard(~ {
                all(.x %in% to_exclude)
              })
          }

          rstatix::wilcox_test(
            data = prep,
            formula = formula(paste(.x, "~", var)),
            comparisons = comparisons,
            p.adjust.method = p_adj_method,
            detailed = TRUE
          ) %>%
            dplyr::filter(p < max_significance) %>%
            tidyr::unite("comparison", group1:group2, sep = "_") %>%
            dplyr::rename(padj = p, taxa_id = .y.) %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id")
        })
    })
}
