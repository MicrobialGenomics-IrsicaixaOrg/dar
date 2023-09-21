#' Phyloseq Quality Control Metrics
#'
#' phy_qc() returns a tibble. It will have information about some important
#' metrics about the sparsity of the count matrix. The content of the table is
#' as follows:
#' * var_levels: levels of the categorical variable of interest. "all" refers 
#' to all rows of the dataset (without splitting by categorical levels).
#' * n: total number of values in the count matrix.
#' * n_zero: number of zeros in the count matrix.
#' * pct_zero: percentage of zeros in the count matrix.
#' * pct_all_zero: percentage of taxa with zero counts in all samples.
#' * pct_singletons: percentage of taxa with counts in a single sample.
#' * pct_doubletons: percentage of taxa with counts in two samples.
#' * count_mean: average of the mean counts per sample.
#' * count_min: average of the min counts per sample.
#' * count_max: average of the max counts per sample.
#'
#' @param rec A recipe or recipe step.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' phy_qc(rec)
methods::setGeneric("phy_qc", function(rec) standardGeneric("phy_qc"))

#' @rdname phy_qc
#' @export
#' @autoglobal
methods::setMethod(
  f = "phy_qc",
  signature = "recipe",
  definition = function(rec) {
    .zero_stats(rec) %>% 
      dplyr::left_join(.zero_groups(rec), by = get_var(rec)[[1]]) %>% 
      dplyr::left_join(.count_summary(rec), by = get_var(rec)[[1]]) %>% 
      dplyr::rename(var_levels = !!dplyr::sym(get_var(rec)[[1]]))
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
.zero_prepro <- function(rec) {
  var <- get_var(rec)[[1]]
  otu_table(rec) %>%
    tidyr::pivot_longer(-1, names_to = "sample_id") %>%
    dplyr::left_join(sample_data(rec), by = "sample_id") %>% 
    dplyr::bind_rows(., dplyr::mutate(., !!var := "all"))
}

#' @noRd
#' @keywords internal
#' @autoglobal
.zero_stats <- function(rec) {
  var <- get_var(rec)[[1]]
  .zero_prepro(rec) %>%
    dplyr::group_by(!!dplyr::sym(var)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_zero = sum(value == 0),
      pct_zero = n_zero / n * 100
    )
}
  
#' @noRd
#' @keywords internal
#' @autoglobal
.zero_groups <- function(rec) {
  var <- get_var(rec)[[1]]
  .zero_prepro(rec) %>% 
    dplyr::mutate(no_zero = dplyr::if_else(value == 0, FALSE, TRUE)) %>% 
    dplyr::group_by(!!dplyr::sym(var), taxa_id) %>% 
    dplyr::summarise(
      all = sum(no_zero) == 0,
      single = sum(no_zero) == 1,
      double = sum(no_zero) == 2, 
      .groups = "drop_last"
    ) %>% 
    dplyr::summarise(
      pct_all_zero = sum(all) / dplyr::n() * 100, 
      pct_singletons = sum(single) / dplyr::n() * 100, 
      pct_doubletons = sum(double) / dplyr::n() * 100, 
      .groups = "drop"
    )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.count_summary <- function(rec) {
  .zero_prepro(rec) %>% 
    dplyr::group_by(!!dplyr::sym(get_var(rec)[[1]]), sample_id) %>% 
    dplyr::summarise(
      count_mean = mean(value), 
      count_min = min(value),
      count_max = max(value),
      .groups = "drop_last", 
    ) %>% 
    dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x)))
}
