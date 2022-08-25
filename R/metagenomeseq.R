#' MetagenomeSeq analysis
#'
#' metagenomeSeq is designed to determine features (be it Operational Taxanomic
#' Unit (OTU), species, etc.) that are differentially abundant between two or
#' more groups of multiple samples. metagenomeSeq is designed to address the
#' effects of both normalization and under-sampling of microbial communities on
#' disease association detection and the testing of feature correlations.
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param zeroMod The zero model, the model to account for the change in the
#'   number of OTUs observed as a linear effect of the depth of coverage.
#' @param useCSSoffset Boolean, whether to include the default scaling
#'   parameters in the model or not.
#' @param useMixedModel Estimate the correlation between duplicate features or
#'   replicates using duplicateCorrelation.
#' @param max_significance The q-value threshold for significance.
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
#' @aliases step_metagenomeseq
#' @return An object of class `recipe`
#' @export
#' @examples 
#' data(metaHIV_phy)
#' 
#' ## Init recipe
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#' 
#' ## Define metagenomeSeq step with default parameters
#' rec <- step_metagenomeseq(rec)
#' rec
methods::setGeneric(
  name = "step_metagenomeseq",
  def = function(rec,
                 zeroMod = NULL,
                 useCSSoffset = TRUE,
                 useMixedModel = FALSE,
                 max_significance = 0.05,
                 log2FC = 0,
                 rarefy = FALSE,
                 id = rand_id("metagenomeseq")) {
    standardGeneric("step_metagenomeseq")
  }
)

#' @rdname step_metagenomeseq
#' @export
methods::setMethod(
  f = "step_metagenomeseq",
  signature = c(rec = "recipe"),
  definition = function(rec,
                        zeroMod,
                        useCSSoffset,
                        useMixedModel,
                        max_significance,
                        log2FC,
                        rarefy,
                        id) {

    recipes_pkg_check(required_pkgs_metagenomeseq(), "step_matagenomeseq()")
    add_step(
      rec,
      step_metagenomeseq_new(
        zeroMod = zeroMod,
        useCSSoffset = useCSSoffset,
        useMixedModel = useMixedModel,
        max_significance = max_significance,
        log2FC = log2FC,
        rarefy = rarefy,
        id = id
      )
    )
  }
)

#' @rdname step_metagenomeseq
#' @export
methods::setMethod(
  f = "step_metagenomeseq",
  signature = c(rec = "prep_recipe"),
  definition = function(rec,
                        zeroMod,
                        useCSSoffset,
                        useMixedModel,
                        max_significance,
                        log2FC,
                        rarefy,
                        id) {
    rlang::abort("This function needs a non-prep recipe!")
  }
)


#' @noRd
#' @keywords internal
step_metagenomeseq_new <- function(rec,
                                   zeroMod,
                                   useCSSoffset,
                                   useMixedModel,
                                   max_significance,
                                   log2FC,
                                   rarefy,
                                   id) {
  step(
    subclass = "metagenomeseq",
    zeroMod = zeroMod,
    useCSSoffset = useCSSoffset,
    useMixedModel = useMixedModel,
    max_significance = max_significance,
    log2FC = log2FC,
    rarefy = rarefy,
    id = id
  )
}

#' @noRd
#' @keywords internal
required_pkgs_metagenomeseq <-
  function(x, ...) {
    c("bioc::metagenomeSeq", "bioc::limma", "bioc::Biobase")
  }

#' @noRd
#' @keywords internal
run_metagenomeseq <- function(rec,
                              zeroMod,
                              useCSSoffset,
                              useMixedModel,
                              max_significance,
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
      vct_var <- phyloseq::sample_data(phy)[[var]]
      model <- stats::model.matrix(~ vct_var)
      colnames(model) <- levels(factor(vct_var))
      mr_obj <-
        phyloseq_to_MRexperiment(phy) %>%
        metagenomeSeq::cumNorm(p = metagenomeSeq::cumNormStat(.)) %>%
        metagenomeSeq::fitZig(
          mod = model,
          zeroMod = zeroMod,
          useCSSoffset = useCSSoffset,
          useMixedModel = useMixedModel,
          control = metagenomeSeq::zigControl(verbose = FALSE)
        )

      contrasts_df <-
        get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
        purrr::map_chr(~ stringr::str_c(.x, collapse = " - "))

      f_fit <-
        limma::makeContrasts(
          contrasts = contrasts_df,
          levels = methods::slot(mr_obj, "fit") %>% purrr::pluck("design")
        ) %>%
        limma::contrasts.fit(fit = methods::slot(mr_obj, "fit"), contrasts = .) %>%
        limma::eBayes()

      contrasts_df %>%
        purrr::imap_dfr(~ {
          limma::topTable(f_fit, coef = .y, number = Inf) %>%
            tibble::as_tibble(rownames = "taxa_id") %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::rename(pvalue = P.Value, padj = adj.P.Val) %>%
            dplyr::mutate(
              comparison = stringr::str_replace_all(.x, " - ", "_"), 
              var = var
            ) %>%
            dplyr::filter(
              abs(.data$logFC) >= log2FC & .data$padj < max_significance
            )
        })
    })
}

#' Phyloseq to MRexpereiment object
#'
#' Wrapper to convert Phyloseq objects to MRexpereiment objects.
#'
#' @param phy Phyloseq object
#'
#' @noRd
#' @keywords internal
phyloseq_to_MRexperiment <- function(phy) {
  counts <- phyloseq::otu_table(phy) %>% data.frame(check.names = FALSE)

  phenoData <-
    phyloseq::sample_data(phy) %>%
    Biobase::AnnotatedDataFrame()

  otuData <-
    phyloseq::tax_table(phy) %>%
    data.frame() %>%
    Biobase::AnnotatedDataFrame()

  metagenomeSeq::newMRexperiment(counts, phenoData, otuData)
}
