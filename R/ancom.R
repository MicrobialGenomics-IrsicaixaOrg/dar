#' ANCOM analysis
#'
#' This function performs an analysis of composition of microbiomes (ANCOM).
#' ANCOM accounts for the underlying structure in the data and can be used for
#' comparing the composition of microbiomes in two or more populations. ANCOM
#' makes no distributional assumptions and can be implemented in a linear model
#' framework to adjust for covariates as well as model longitudinal data.
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param formula the character string expresses how the microbial absolute
#'   abundances for each taxon depend on the variables in metadata.
#' @param p_adj_method method to adjust p-values by. Default is "holm". Options
#'   include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#'   "none". See p.adjust for more details.
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of
#'   zeroes greater than zero_cut will be excluded in the analysis. Default is
#'   0.90.
#' @param lib_cut a numerical threshold for filtering samples based on library
#'   sizes. Samples with library sizes less than lib_cut will be excluded in the
#'   analysis. Default is 0, i.e. do not filter any sample.
#' @param group the name of the group variable in metadata. Specifying group is
#'   required for detecting structural zeros and performing global test.
#' @param struc_zero whether to detect structural zeros. Default is FALSE.
#' @param neg_lb whether to classify a taxon as a structural zero in the
#'   corresponding study group using its asymptotic lower bound. Default is
#'   FALSE.
#' @param tol the iteration convergence tolerance for the E-M algorithm. Default
#'   is 1e-05.
#' @param max_iter the maximum number of iterations for the E-M algorithm.
#'   Default is 100.
#' @param conserve whether to use a conservative variance estimate of the test
#'   statistic. It is recommended if the sample size is small and/or the number
#'   of differentially abundant taxa is believed to be large. Default is FALSE.
#' @param alpha level of significance. Default is 0.05.
#' @param global whether to perform global test. Default is FALSE.
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
#' @aliases step_ancom
#' @return An object of class `recipe`
#' @export
#' @examples 
#' data(metaHIV_phy)
#' 
#' ## Init recipe
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
#'   step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.4 * length(x))")
#' 
#' rec
#' 
#' ## Define step with default parameters and prep
#' rec <- 
#'   step_ancom(rec) %>% 
#'   prep(parallel = TRUE)
#'   
#' rec
#' 
#' ## Wearing rarefaction only for this step 
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
#'   step_ancom(rec, rarefy = TRUE)
#' 
#' rec
methods::setGeneric(
  name = "step_ancom",
  def = function(rec,
                 formula = get_var(rec),
                 p_adj_method = "holm",
                 zero_cut = 0.9,
                 lib_cut = 0,
                 group = NULL,
                 struc_zero = FALSE,
                 neg_lb = FALSE,
                 tol = 1e-05, 
                 max_iter = 100, 
                 conserve = FALSE, 
                 alpha = 0.05, 
                 global = FALSE, 
                 rarefy = FALSE, 
                 id = rand_id("ancom")) {
    standardGeneric("step_ancom")
  }
)

#' @rdname step_ancom
#' @export
methods::setMethod(
  f = "step_ancom",
  signature = c(rec = "recipe"),
  definition = function(rec,
                        formula,
                        p_adj_method,
                        zero_cut,
                        lib_cut,
                        group,
                        struc_zero,
                        neg_lb,
                        tol, 
                        max_iter, 
                        conserve, 
                        alpha, 
                        global, 
                        rarefy, 
                        id) {

    recipes_pkg_check(required_pkgs_ancom(), "step_ancom()")
    add_step(
      rec,
      step_ancom_new(
        formula = formula,
        p_adj_method = p_adj_method,
        zero_cut = zero_cut,
        lib_cut = lib_cut,
        group = group,
        struc_zero = struc_zero,
        neg_lb = neg_lb,
        tol = tol, 
        max_iter = max_iter, 
        conserve = conserve, 
        alpha = alpha, 
        global = global, 
        rarefy = rarefy, 
        id = id
      )
    )
  }
)

#' @noRd
#' @keywords internal
step_ancom_new <-
  function(formula,
           p_adj_method,
           zero_cut,
           lib_cut,
           group,
           struc_zero,
           neg_lb,
           tol, 
           max_iter, 
           conserve, 
           alpha, 
           global, 
           rarefy, 
           id) {
    step(
      subclass = "ancom",
      formula = formula,
      p_adj_method = p_adj_method,
      zero_cut = zero_cut,
      lib_cut = lib_cut,
      group = group,
      struc_zero = struc_zero,
      neg_lb = neg_lb,
      tol = tol, 
      max_iter = max_iter, 
      conserve = conserve, 
      alpha = alpha, 
      global = global, 
      rarefy = rarefy, 
      id = id
    )
  }

#' @noRd
#' @keywords internal
required_pkgs_ancom <- function(x, ...) { c("bioc::ANCOMBC") }

#' @noRd
#' @keywords internal
run_ancom <- function(rec,
                      formula,
                      p_adj_method,
                      zero_cut,
                      lib_cut,
                      group,
                      struc_zero,
                      neg_lb,
                      tol, 
                      max_iter, 
                      conserve, 
                      alpha, 
                      global, 
                      rarefy, 
                      id) {
  
  phy <- get_phy(rec)
  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  if (rarefy) {
    phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE)
  }
  
  phy <- phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
        purrr::map_dfr(function(comparison) {
          vct_comp <- glue::glue("c('{comparison[1]}', '{comparison[2]}')")
          s_phy <-
            glue::glue(
              "phyloseq::subset_samples(phy, {var} %in% {vct_comp})"
            ) %>%
            parse(text = .) %>%
            eval()
          
          res <- ANCOMBC::ancombc(
            phyloseq = s_phy,
            formula = var,
            p_adj_method = p_adj_method,
            zero_cut = zero_cut,
            lib_cut = lib_cut,
            group = group,
            struc_zero = struc_zero,
            neg_lb = neg_lb,
            tol = tol,
            max_iter = max_iter,
            conserve = conserve,
            alpha = alpha,
            global = global
          )
          
          res$res %>%
            purrr::map2_dfc(names(res$res), ~ {
              tibble::tibble(!!.y := .x[[1]])
            }) %>%
            dplyr::mutate(
              taxa_id = rownames(res$res$beta),
              .before = 1,
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = var
            ) %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::mutate(
              effect = .data$beta,
              signif = .data$diff_abn
            )
        })
    }) 
}
