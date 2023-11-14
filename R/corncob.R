#' corncob analysis
#'
#' Corncob an individual taxon regression model that uses abundance tables and
#' sample data. corncob is able to model differential abundance and differential
#' variability, and addresses each of the challenges presented below:
#' 
#' * different sequencing depth
#' * excessive zeros from unobserved taxa
#' * high variability of empirical relative abundances (overdispersion)
#' * within-taxon correlation
#' * hypothesis testing with categorical and continuous covariates
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param phi.formula An object of class formula without the response: a
#'   symbolic description of the model to be fitted to the dispersion.
#' @param formula_null Formula for mean under null, without response.
#' @param phi.formula_null Formula for overdispersion under null, without
#'   response.
#' @param link Link function for abundance covariates, defaults to "logit".
#' @param phi.link Link function for dispersion covariates, defaults to "logit".
#' @param test Character. Hypothesis testing procedure to use. One of "Wald" or
#'   "LRT" (likelihood ratio test).
#' @param boot Boolean. Defaults to FALSE. Indicator of whether or not to use
#'   parametric bootstrap algorithm. (See pbWald and pbLRT).
#' @param B Optional integer. Number of bootstrap iterations. Ignored if boot is
#'   FALSE. Otherwise, defaults to 1000.
#' @param filter_discriminant Boolean. Defaults to TRUE. If FALSE, discriminant
#'   taxa will not be filtered out.
#' @param fdr_cutoff Integer. Defaults to 0.05. Desired type 1 error rate.
#' @param fdr Character. Defaults to "fdr". False discovery rate control method,
#'   see p.adjust for more options.
#' @param log2FC log2FC cutoff.
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
#' @aliases step_corncob
#' @return An object of class `recipe`
#' @export
#' @autoglobal
#' @md
#' @tests testthat
#' data(test_prep_rec)
#' expect_error(step_corncob(test_prep_rec, rarefy = TRUE))
#' 
#' @examples 
#' data(metaHIV_phy)
#' 
#' ## Init recipe
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Class") %>% 
#'   step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))")
#' 
#' rec
#' 
#' ## Define step with default parameters and prep
#' rec <- 
#'   step_corncob(rec) %>% 
#'   prep(parallel = FALSE)
#'   
#' rec
#' 
#' ## Wearing rarefaction only for this step 
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
#'   step_corncob(rec, rarefy = TRUE)
#' 
#' rec
methods::setGeneric(
  name = "step_corncob",
  def = function(rec,
                 phi.formula = stats::formula(~ 1),
                 formula_null = stats::formula(~ 1),
                 phi.formula_null = stats::formula(~ 1),
                 link = "logit",
                 phi.link = "logit",
                 test = "Wald",
                 boot = FALSE,
                 B = 1000,
                 filter_discriminant = TRUE,
                 fdr_cutoff = 0.05,
                 fdr = "fdr",
                 log2FC = 0,
                 rarefy = FALSE,
                 id = rand_id("corncob")) {
    standardGeneric("step_corncob")
  }
)

#' @rdname step_corncob
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_corncob",
  signature = c(rec = "recipe"),
  definition = function(rec,
                        phi.formula,
                        formula_null,
                        phi.formula_null,
                        link,
                        phi.link,
                        test,
                        boot,
                        B,
                        filter_discriminant,
                        fdr_cutoff,
                        fdr,
                        log2FC,
                        rarefy,
                        id) {

    recipes_pkg_check(required_pkgs_corncob(), "step_croncob()")
    add_step(
      rec,
      step_corncob_new(
        phi.formula = phi.formula,
        formula_null = formula_null,
        phi.formula_null = phi.formula_null,
        link = link,
        phi.link = phi.link,
        test = test,
        boot = boot,
        B = B,
        filter_discriminant = filter_discriminant,
        fdr_cutoff = fdr_cutoff,
        fdr = fdr,
        log2FC = log2FC,
        rarefy = rarefy,
        id = id
      ))
  }
)

#' @rdname step_corncob
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_corncob",
  signature = c(rec = "prep_recipe"),
  definition = function(rec,
                        phi.formula,
                        formula_null,
                        phi.formula_null,
                        link,
                        phi.link,
                        test,
                        boot,
                        B,
                        filter_discriminant,
                        fdr_cutoff,
                        fdr,
                        log2FC,
                        rarefy,
                        id) {
    rlang::abort("This function needs a non-prep recipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_corncob_new <- function(phi.formula,
                             formula_null,
                             phi.formula_null,
                             link,
                             phi.link,
                             test,
                             boot,
                             B,
                             filter_discriminant,
                             fdr_cutoff,
                             fdr,
                             log2FC,
                             rarefy,
                             id) {
  step(
    subclass = "corncob",
    phi.formula = phi.formula,
    formula_null = formula_null,
    phi.formula_null = phi.formula_null,
    link = link,
    phi.link = phi.link,
    test = test,
    boot = boot,
    B = B,
    filter_discriminant = filter_discriminant,
    fdr_cutoff = fdr_cutoff,
    fdr = fdr,
    log2FC = log2FC,
    rarefy = rarefy,
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_corncob <- function(x, ...) { c("corncob") }

#' @noRd
#' @keywords internal
#' @autoglobal
run_corncob <- function(rec,
                        phi.formula,
                        formula_null,
                        phi.formula_null,
                        link,
                        phi.link,
                        test,
                        boot,
                        B,
                        filter_discriminant,
                        fdr_cutoff,
                        fdr,
                        log2FC,
                        rarefy) {
  
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

          ## Filter samples
          f_phy <- 
            phyloseq::sample_data(phy) %>% 
            to_tibble("sample_id") %>% 
            dplyr::filter(!!dplyr::sym(var) %in% comparison) %>% 
            dplyr::pull(sample_id) %>% 
            phyloseq::prune_samples(phy)
          
          corncob_res <- 
            tryCatch({
              corncob::differentialTest(
                formula = glue::glue("~ { var }") %>% stats::formula(),
                data = f_phy,
                phi.formula = phi.formula,
                formula_null = formula_null,
                phi.formula_null = phi.formula_null,
                link = link,
                phi.link = phi.link,
                test = test,
                boot = boot,
                B = 10e4,
                filter_discriminant = filter_discriminant,
                fdr_cutoff = Inf,
                fdr = fdr
              )
            }, error = function(e) { conditionMessage(e) })
        
          ## Skip error for no convergence
          check_non_convergence(corncob_res)
         
          signif_taxa <- corncob::otu_to_taxonomy(
            OTU = corncob_res$significant_taxa,
            data = corncob_res$data,
            level = tax_level %>% dplyr::pull()
          ) %>% stringr::str_c(" (", corncob_res$significant_taxa, ")")

          corncob_stats_tbl(
            corncob_res,
            signif_taxa,
            comparison,
            var,
            log2FC,
            fdr_cutoff
          )
        })
    })
}

#' @noRd
#' @keywords internal
#' @autoglobal
corncob_stats_tbl <- function(corncob_res,
                              signif_taxa,
                              comparison,
                              var,
                              log2FC,
                              fdr_cutoff) {
    
  qval <- stats::qnorm(0.975)
  corncob_res$significant_models %>%
    purrr::map2_dfr(stats::na.omit(corncob_res$p_fdr), function(m, p) {
      coefs_mu <- m$coefficients[seq_len(m$np.mu), , drop = FALSE]
      coefs_mu <- coefs_mu[2, , drop = FALSE]
      tibble::tibble(
        log2FC = coefs_mu[1, 1],
        log2FC_min = coefs_mu[1, 1] - qval * coefs_mu[1, 2],
        log2FC_xmax = coefs_mu[1, 1] + qval * coefs_mu[1, 2],
        padj = p
      )
    }) %>%
    dplyr::mutate(
      taxa = signif_taxa,
      comparison = stringr::str_c(comparison, collapse = "_"),
      var = var
    ) %>%
    tidyr::separate(taxa, c("taxa", "taxa_id"), sep = " ", remove = TRUE) %>%
    dplyr::mutate(taxa_id = stringr::str_remove_all(taxa_id, "[(]|[)]")) %>%
    dplyr::mutate(
      effect = log2FC,
      signif = ifelse(padj < fdr_cutoff & abs(log2FC) >= log2FC, TRUE, FALSE)
    )
}

#' @noRd
#' @keywords internal
#' @autoglobal
check_non_convergence <- function(corncob_res) {
  if (!methods::is(corncob_res, "differentialTest")) {
    if (stringr::str_detect(corncob_res, "failed to converge")) {
      rlang::abort(c(
        glue::glue(
          "{crayon::bgMagenta('corncob')}: All models failed to converge!"
        ),
        glue::glue(
          "{crayon::bgMagenta('corncob')}: If you are seeing this, it",
          " is likely that your model is overspecified. This occurs ", 
          "when your sample size is not large enough to estimate all ", 
          "the parameters of your model. This is most commonly due to", 
          " categorical variables that include many categories."
        ),
        glue::glue(
          "Please remove or edit the ", 
          "{crayon::bgMagenta('step_corncob()')} and rerun."
        )
      ))
    } else {
      rlang::abort(
        c(
          glue::glue("{crayon::bgMagenta('corncob')}: Internal error!"),
          glue::glue(
            "Please remove or edit the ", 
            "{crayon::bgMagenta('step_corncob()')} and rerun."
          ),
          glue::glue(
            "Please report this bug on GitHub: ", 
            "https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/", 
            "issues" 
          )
        )
      )
    }
  }
}
