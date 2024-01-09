#' DESeq2 analysis
#'
#' Differential expression analysis based on the Negative Binomial (a.k.a.
#' Gamma-Poisson) distribution. This function performs a default analysis
#' through the steps: 1) estimation of size factors: `estimateSizeFactors`. 2)
#' estimation of dispersion: `estimateDispersions`. 3) Negative Binomial GLM
#' fitting and Wald statistics: nbinomWaldTest. For complete details on each
#' step, see the manual pages of the respective functions. After the DESeq
#' function returns a DESeqDataSet object, results tables (log2 fold changes and
#' p-values) can be generated using the results function. Shrunken LFC can then
#' be generated using the lfcShrink function.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param test Either "Wald" or "LRT", which will then use either Wald
#'   significance tests (defined by nbinomWaldTest), or the likelihood ratio
#'   test on the difference in deviance between a full and reduced model formula
#'   (defined by nbinomLRT).
#' @param fitType either "parametric", "local", "mean", or "glmGamPoi" for the
#'   type of fitting of dispersions to the mean intensity. See
#'   estimateDispersions for description.
#' @param betaPrior whether or not to put a zero-mean normal prior on the
#'   non-intercept coefficients See nbinomWaldTest for description of the
#'   calculation of the beta prior. In versions >=1.16, the default is set to
#'   FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
#' @param type "apeglm" is the adaptive Student's t prior shrinkage estimator
#'   from the 'apeglm' package; "ashr" is the adaptive shrinkage estimator from
#'   the 'ashr' package, using a fitted mixture of normals prior - see the
#'   Stephens (2016) reference below for citation; "normal" is the 2014 DESeq2
#'   shrinkage estimator using a Normal prior.
#' @param max_significance The q-value threshold for significance.
#' @param log2FC log2FC cutoff.
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
#' @aliases step_deseq
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' 
#' test <-
#'  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_rarity(0.95) |> 
#'  step_deseq() |> 
#'  prep()
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_deseq(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Init Recipe
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Class") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.4 * length(x))")
#'
#' rec
#'
#' ## Define step with default parameters and prep
#' rec <-
#'   step_deseq(rec) |>
#'   prep(parallel = FALSE)
#'
#' rec
#'
#' ## Wearing rarefaction only for this step
#' rec <-
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_deseq(rarefy = TRUE)
#'
#' rec
methods::setGeneric(
  name = "step_deseq",
  def = function(rec,
                 test = "Wald",
                 fitType = "local",
                 betaPrior = FALSE,
                 type = "ashr",
                 max_significance = 0.05,
                 log2FC = 0,
                 rarefy = FALSE,
                 id = rand_id("deseq")) {
    standardGeneric("step_deseq")
  }
)

#' @rdname step_deseq
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_deseq",
  signature = c(rec = "Recipe"),
  definition = function(rec,
                        test,
                        fitType,
                        betaPrior,
                        type,
                        max_significance,
                        log2FC,
                        rarefy,
                        id) {

    if (type == "ashr") {
      recipes_pkg_check(required_pkgs_deseq()[-3], "step_deseq()")
    } else {
      recipes_pkg_check(required_pkgs_deseq()[-2], "step_deseq()")
    }

    add_step(
      rec,
      step_deseq_new(
        test = test,
        fitType = fitType,
        betaPrior = betaPrior,
        type = type,
        max_significance = max_significance,
        log2FC = log2FC,
        rarefy = rarefy,
        id = id
      )
    )
  }
)

#' @rdname step_deseq
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_deseq",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec,
                        test,
                        fitType,
                        betaPrior,
                        type,
                        max_significance,
                        log2FC,
                        rarefy,
                        id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_deseq_new <-
  function(rec,
           test,
           fitType,
           betaPrior,
           type,
           max_significance,
           log2FC,
           rarefy,
           id) {

  step(
    subclass = "deseq",
    test = test,
    fitType = fitType,
    betaPrior = betaPrior,
    type = type,
    max_significance = max_significance,
    log2FC = log2FC,
    rarefy = rarefy,
    id = id
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_deseq <- function(x, ...) {
  c("bioc::DESeq2", "bioc::apeglm", "ashr")
}

#' @noRd
#' @keywords internal
#' @autoglobal
run_deseq <- function(rec,
                      test,
                      fitType,
                      betaPrior,
                      type,
                      max_significance,
                      log2FC,
                      rarefy) {

  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  phy <- 
    get_phy(rec) %>% 
    use_rarefy(rarefy)
  
  phy <- phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)
  phyloseq::sample_data(phy) <-  
    to_tibble(phyloseq::sample_data(phy), "sample_id") %>% 
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>% 
    data.frame(row.names = 1) %>% 
    phyloseq::sample_data() 
    
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      dds <- suppressMessages(
        phy %>%
          phyloseq::phyloseq_to_deseq2(
            design = stats::as.formula(stringr::str_c("~", var))
          )
      )
      
      dds <-
        dds %>%
        DESeq2::estimateSizeFactors(geoMeans = apply(
          X = DESeq2::counts(.),
          MARGIN = 1,
          FUN = function(x) {
            exp(sum(log(x[x > 0]), na.rm = TRUE) / length(x))
          }
        )) %>%
        DESeq2::DESeq(
          fitType = fitType,
          test = test,
          betaPrior = betaPrior,
          quiet = TRUE
        )

      contrasts_df <- get_comparisons(var, phy, as_list = FALSE, n_cut = 1)
      seq_len(nrow(contrasts_df)) %>%
        purrr::map_dfr(function(it) {
          x <- contrasts_df[[it, "x"]]
          y <- contrasts_df[[it, "y"]]
          DESeq2::lfcShrink(
            dds = dds,
            contrast = c(var, y, x),
            type = type,
            quiet = TRUE
          ) %>%
            tibble::as_tibble(rownames = "taxa_id") %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::mutate(
              comparison = stringr::str_c(x, "_", y), var = !!var
            ) %>%
            dplyr::mutate(
              effect = log2FoldChange,
              signif = ifelse(
                padj < max_significance & abs(log2FoldChange) >= log2FC,
                TRUE,
                FALSE
              )
            )
        })
    })
}
