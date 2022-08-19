#' ANCOM analysis
#'
#' This function performs an analysis of composition of microbiomes (ANCOM). ANCOM
#' accounts for the underlying structure in the data and can be used for comparing the
#' composition of microbiomes in two or more populations. ANCOM makes no distributional
#' assumptions and can be implemented in a linear model framework to adjust for covariates
#' as well as model longitudinal data.
#'
#' @param rec A recipe object. The step will be added to the sequence of operations for
#'   this recipe.
#' @param out_cut Numerical fraction between 0 and 1. For each taxon, observations with
#'   proportion of mixture distribution less than out_cut will be detected as outlier
#'   zeros; while observations with proportion of mixture distribution greater than 1 -
#'   out_cut will be detected as outlier values.
#' @param zero_cut Numerical fraction between 0 and 1. Taxa with proportion of zeroes
#'   greater than zero_cut are not included in the analysis.
#' @param lib_cut Numeric. Samples with library size less than lib_cut are not included in
#'   the analysis.
#' @param neg_lb Logical. TRUE indicates a taxon would be classified as a structural zero
#'   in the corresponding experimental group using its asymptotic lower bound. More
#'   specifically, neg_lb = TRUE indicates you are using both criteria stated in section
#'   3.2 of ANCOM-II to detect structural zeros; Otherwise, neg_lb = FALSE will only use
#'   the equation 1 in section 3.2 of ANCOM-II for declaring structural zeros.
#' @param p_adj_method Character. Specifying the method to adjust p-values for multiple
#'   comparisons. Default is “BH” (Benjamini-Hochberg procedure).
#' @param alpha Level of significance. Default is 0.05.
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
#' @aliases step_ancom
#' @return An object of class `recipe`
#' @export
methods::setGeneric(
  name = "step_ancom",
  def = function(rec,
                 out_cut = 0.05,
                 zero_cut = 0.9,
                 lib_cut = 0,
                 neg_lb = FALSE,
                 p_adj_method = "BH",
                 alpha = 0.05,
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
                        out_cut,
                        zero_cut,
                        lib_cut,
                        neg_lb,
                        p_adj_method,
                        alpha,
                        rarefy,
                        id) {

    recipes_pkg_check(required_pkgs_ancom(), "step_ancom()")
    add_step(
      rec,
      step_ancom_new(
        out_cut = out_cut,
        zero_cut = zero_cut,
        lib_cut = lib_cut,
        neg_lb = neg_lb,
        p_adj_method = p_adj_method,
        alpha = alpha,
        rarefy = rarefy,
        id = id
      )
    )
  }
)

#' @noRd
#' @keywords internal
step_ancom_new <-
  function(out_cut,
           zero_cut,
           lib_cut,
           neg_lb,
           p_adj_method,
           alpha,
           rarefy,
           id) {
    step(
      subclass = "ancom",
      out_cut = out_cut,
      zero_cut = zero_cut,
      lib_cut = lib_cut,
      neg_lb = neg_lb,
      p_adj_method = p_adj_method,
      alpha = alpha,
      rarefy = rarefy,
      id = id
    )
  }

#' @noRd
#' @keywords internal
required_pkgs_ancom <- function(x, ...) { c("xec-cm/ANCOM") }

#' @noRd
#' @keywords internal
run_ancom <- function(rec, out_cut, zero_cut, lib_cut, neg_lb, p_adj_method, alpha, rarefy) {

  phy <- get_phy(rec)
  vars <- get_var(rec)
  tax_level <- get_tax(rec)
  if (rarefy) { phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE) }
  phy <- phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)

  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
        purrr::map_dfr(function(comparison) {
          sample_data <- dplyr::filter(sample_data(rec), !!dplyr::sym(var) %in% comparison)
          otu_table <-
            otu_table(rec) %>%
            dplyr::select(taxa_id, dplyr::all_of(sample_data$sample_id)) %>%
            data.frame(row.names = 1) %>%
            as.matrix()

          prepro <- ANCOM::feature_table_pre_process(
            feature_table = otu_table,
            meta_data = sample_data,
            sample_var = 'sample_id',
            lib_cut = lib_cut,
            out_cut = out_cut,
            zero_cut = zero_cut
          )

          ancom_res <- ANCOM::ANCOM(
            feature_table = prepro$feature_table,
            meta_data = prepro$meta_data,
            struc_zero = prepro$structure_zeros,
            main_var = var,
            p_adj_method = p_adj_method,
            alpha = alpha
          )

          n_taxa <- ifelse(
            is.null(prepro$struc_zero),
            nrow(prepro$feature_table),
            sum(apply(prepro$struc_zero, 1, sum) == 0)
          )

          ancom_res$out %>%
            tibble::as_tibble() %>%
            dplyr::mutate(
              clr_mean_diff = ancom_res$fig$data$x,
              structural_zero = ancom_res$fig$data$zero_ind,
              x_label = min(ancom_res$fig$data$x),
              y_label = 0.7 * (n_taxa - 1),
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = var
            ) %>%
            dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
            dplyr::filter(.data$detected_0.7 == TRUE)
        })
    })
}
