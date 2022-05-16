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
#' @param struc_zero A matrix consists of 0 and 1s with 1 indicating the taxon is
#'   identified as a structural zero in the corresponding group. Can be the output value
#'   from feature_table_pre_process.
#' @param p_adj_method Character. Specifying the method to adjust p-values for multiple
#'   comparisons. Default is “BH” (Benjamini-Hochberg procedure).
#' @param alpha Level of significance. Default is 0.05.
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
                 out_cut = 0,
                 zero_cut = 1,
                 lib_cut = 0,
                 neg_lb = FALSE,
                 struc_zero = NULL,
                 p_adj_method = "BH",
                 alpha = 0.05,
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
                        struc_zero,
                        p_adj_method,
                        alpha,
                        id) {

    recipes_pkg_check(required_pkgs_ancom())

    add_step(
      rec,
      step_ancom_new(
        out_cut = out_cut,
        zero_cut = zero_cut,
        lib_cut = lib_cut,
        neg_lb = neg_lb,
        struc_zero = struc_zero,
        p_adj_method = p_adj_method,
        alpha = alpha,
        id = id
      )
    )

  }
)

#' @rdname step_ancom
#' @keywords internal
step_ancom_new <-
  function(out_cut,
           zero_cut,
           lib_cut,
           neg_lb,
           struc_zero,
           p_adj_method,
           alpha,
           id) {
    step(
      subclass = "ancom",
      out_cut = out_cut,
      zero_cut = zero_cut,
      lib_cut = lib_cut,
      neg_lb = neg_lb,
      struc_zero = struc_zero,
      p_adj_method = p_adj_method,
      alpha = alpha,
      id = id
    )
  }

#' @noRd
#' @keywords internal
required_pkgs_ancom <- function(x, ...) { c("xec-cm/ANCOM") }

run_ancom <- function(phy,
                      tax_level,
                      out_cut,
                      zero_cut,
                      lib_cut,
                      neg_lb,
                      struc_zero,
                      p_adj_method,
                      alpha) {

  phy <- phyloseq::tax_glom(phy, taxrank = tax_level, NArm = FALSE)
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, phy, as_list = TRUE, n_cut = 1)



    })

}

