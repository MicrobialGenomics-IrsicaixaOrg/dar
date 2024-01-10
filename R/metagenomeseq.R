#' MetagenomeSeq analysis
#'
#' metagenomeSeq is designed to determine features (be it Operational Taxanomic
#' Unit (OTU), species, etc.) that are differentially abundant between two or
#' more groups of multiple samples. metagenomeSeq is designed to address the
#' effects of both normalization and under-sampling of microbial communities on
#' disease association detection and the testing of feature correlations.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
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
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @param rm_zeros Proportion of samples of the same categorical level with more
#'   than 0 counts.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_metagenomeseq
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(test_prep_rec)
#' expect_error(step_metagenomeseq(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.02 * length(x))")
#'
#' rec
#'
#' ## Define step with default parameters and prep
#' rec <-
#'   step_metagenomeseq(rec, rm_zeros = 0.01) |>
#'   prep(parallel = FALSE)
#'
#' rec
#'
#' ## Wearing rarefaction only for this step
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_metagenomeseq(rarefy = TRUE)
#'
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
                 rm_zeros = 0,
                 id = rand_id("metagenomeseq")) {
    standardGeneric("step_metagenomeseq")
  }
)

#' @rdname step_metagenomeseq
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_metagenomeseq",
  signature = c(rec = "Recipe"),
  definition = function(rec,
                        zeroMod,
                        useCSSoffset,
                        useMixedModel,
                        max_significance,
                        log2FC,
                        rarefy,
                        rm_zeros, 
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
        rm_zeros = rm_zeros, 
        id = id
      )
    )
  }
)

#' @rdname step_metagenomeseq
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_metagenomeseq",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec,
                        zeroMod,
                        useCSSoffset,
                        useMixedModel,
                        max_significance,
                        log2FC,
                        rarefy,
                        rm_zeros,
                        id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)


#' @noRd
#' @keywords internal
#' @autoglobal
step_metagenomeseq_new <- function(rec,
                                   zeroMod,
                                   useCSSoffset,
                                   useMixedModel,
                                   max_significance,
                                   log2FC,
                                   rarefy,
                                   rm_zeros, 
                                   id) {
  step(
    subclass = "metagenomeseq",
    zeroMod = zeroMod,
    useCSSoffset = useCSSoffset,
    useMixedModel = useMixedModel,
    max_significance = max_significance,
    log2FC = log2FC,
    rarefy = rarefy,
    rm_zeros = rm_zeros, 
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_metagenomeseq <-
  function(x, ...) {
    c("bioc::metagenomeSeq", "bioc::limma", "bioc::Biobase")
  }
#' @noRd
#' @keywords internal
#' @autoglobal
run_metagenomeseq <- function(rec,
                              zeroMod,
                              useCSSoffset,
                              useMixedModel,
                              max_significance,
                              log2FC,
                              rarefy, 
                              rm_zeros) {
  
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
          
          ## Filter samples
          phy2 <- 
            phyloseq::sample_data(phy) %>% 
            to_tibble("sample_id") %>% 
            dplyr::filter(!!dplyr::sym(var) %in% comparison) %>% 
            dplyr::pull(sample_id) %>% 
            phyloseq::prune_samples(phy)
          
          ## Filter taxa
          to_prune <-  
            zero_otu(phy2, var = var, pct_cutoff = rm_zeros) %>%
            dplyr::pull(taxa_id) %>%
            unique() 
          
          if (length(to_prune) > 0) {
            phy2 <- phyloseq::prune_taxa(to_prune, phy2)
          }
      
          mr_exp <- phyloseq_to_MRexperiment(phy2)
          p <- suppressMessages(metagenomeSeq::cumNormStat(mr_exp))
          mr_obj <- metagenomeSeq::cumNorm(mr_exp, p = p)
          
          vct_var <- phyloseq::sample_data(phy2)[[var]]
          norm_factor <- metagenomeSeq::normFactors(mr_obj)
          norm_factor <- log2(norm_factor / stats::median(norm_factor) + 1)
          
          model <- stats::model.matrix( ~ 1 + vct_var + norm_factor)

          metagenomeSeq::fitZig(
            mr_obj,
            mod = model,
            useCSSoffset = useCSSoffset,
            zeroMod = zeroMod,
            useMixedModel = useMixedModel,
            control = metagenomeSeq::zigControl(
              verbose = FALSE,
              maxit = 100,
              dfMethod = "modified"
            )
          ) %>%
            metagenomeSeq::MRfulltable(number = Inf, by = 2) %>%
            tibble::as_tibble(rownames = "taxa_id") %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::rename(pvalue = pvalues, padj = adjPvalues) %>%
            stats::setNames(
              names(.) %>% stringr::str_replace_all("vct_var.*", "vct_var")
            ) %>%
            dplyr::mutate(
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = var
            ) %>% 
            dplyr::mutate(
              effect = vct_var,
              signif = ifelse(padj < max_significance, TRUE, FALSE)
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
#' @autoglobal
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
