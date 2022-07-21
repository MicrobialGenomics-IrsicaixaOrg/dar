#' DESeq2 analysis
#'
#' Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson)
#' distribution. This function performs a default analysis through the steps: 1)
#' estimation of size factors: `estimateSizeFactors`. 2) estimation of dispersion:
#' `estimateDispersions`. 3) Negative Binomial GLM fitting and Wald statistics:
#' nbinomWaldTest. For complete details on each step, see the manual pages of the
#' respective functions. After the DESeq function returns a DESeqDataSet object, results
#' tables (log2 fold changes and p-values) can be generated using the results function.
#' Shrunken LFC can then be generated using the lfcShrink function.
#'
#' @param rec A recipe object. The step will be added to the sequence of operations for
#'   this recipe.
#' @param test Either "Wald" or "LRT", which will then use either Wald significance tests
#'   (defined by nbinomWaldTest), or the likelihood ratio test on the difference in
#'   deviance between a full and reduced model formula (defined by nbinomLRT).
#' @param fitType either "parametric", "local", "mean", or "glmGamPoi" for the type of
#'   fitting of dispersions to the mean intensity. See estimateDispersions for
#'   description.
#' @param betaPrior whether or not to put a zero-mean normal prior on the non-intercept
#'   coefficients See nbinomWaldTest for description of the calculation of the beta prior.
#'   In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained
#'   afterwards using lfcShrink.
#' @param type "apeglm" is the adaptive Student's t prior shrinkage estimator from the
#'   'apeglm' package; "ashr" is the adaptive shrinkage estimator from the 'ashr' package,
#'   using a fitted mixture of normals prior - see the Stephens (2016) reference below for
#'   citation; "normal" is the 2014 DESeq2 shrinkage estimator using a Normal prior.
#' @param max_significance The q-value threshold for significance.
#' @param log2FC log2FC cutoff.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This rarefaction uses
#'   the standard R sample function to resample from the abundance values in the otu_table
#'   component of the first argument, physeq. Often one of the major goals of this
#'   procedure is to achieve parity in total number of counts between samples, as an
#'   alternative to other formal normalization procedures, which is why a single value for
#'   the sample.size is expected.
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_deseq
#' @return An object of class `recipe`
#' @export
methods::setGeneric(
  name = "step_deseq",
  def = function(rec,
                 test = "Wald",
                 fitType = "parametric",
                 betaPrior = FALSE,
                 type = "ashr",
                 max_significance = 0.05,
                 log2FC = 1,
                 rarefy = FALSE,
                 id = rand_id("deseq")) {
    standardGeneric("step_deseq")
  }
)

#' @rdname step_deseq
#' @export
methods::setMethod(
  f = "step_deseq",
  signature = c(rec = "recipe"),
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
#' @keywords internal
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
required_pkgs_deseq <- function(x, ...) { c("bioc::DESeq2", "bioc::apeglm", "ashr") }

#' @rdname step_deseq
#' @keywords internal
run_deseq <- function(rec, test, fitType, betaPrior, type, max_significance, log2FC, rarefy) {

  phy <- get_phy(rec)
  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  if (rarefy) { phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE) }
  phy <-  phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)

  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      dds <-
        suppressWarnings(
          suppressMessages(
            phy %>%
              phyloseq::phyloseq_to_deseq2(design = stats::as.formula(stringr::str_c("~", var))) %>%
              DESeq2::estimateSizeFactors(geoMeans = apply(
                X = DESeq2::counts(.),
                MARGIN = 1,
                FUN = function(x) {
                  exp(sum(log(x[x > 0]), na.rm = T) / length(x))
                }
              )) %>%
              DESeq2::DESeq(
                fitType = fitType,
                test = test,
                betaPrior = betaPrior,
                quiet = TRUE
              )
          )
        )

      contrasts_df <- get_comparisons(var, phy, as_list = FALSE, n_cut = 1)
      seq_len(nrow(contrasts_df)) %>%
        purrr::map_dfr(function(it) {
          x <- contrasts_df[[it, "x"]]
          y <- contrasts_df[[it, "y"]]
          DESeq2::lfcShrink(
            dds = dds,
            contrast = c(var, x, y),
            type = type,
            quiet = TRUE
          ) %>%
            tibble::as_tibble(rownames = "taxa_id") %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::mutate(comparison = stringr::str_c(x, "_vs_", y), var = !!var) %>%
            dplyr::filter(abs(.data$log2FoldChange) >= log2FC & .data$padj < max_significance)
        })
    })
}
