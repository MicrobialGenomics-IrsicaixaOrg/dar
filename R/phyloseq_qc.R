#' Phyloseq Quality Control Metrics
#'
#' `phy_qc()` returns a tibble containing important metrics about the sparsity 
#' and sequencing depth of the count matrix. These metrics are calculated both 
#' globally and split by the levels of the categorical variable of interest.
#' The content of the table is as follows:
#' 
#' * `var_levels`: levels of the categorical variable of interest. "all" refers 
#'   to all rows of the dataset (without splitting by categorical levels).
#' * `n`: total number of values in the count matrix.
#' * `n_zero`: number of zeros in the count matrix.
#' * `pct_zero`: percentage of zeros in the count matrix.
#' * `pct_all_zero`: percentage of taxa with zero counts in all samples.
#' * `pct_singletons`: percentage of taxa with counts in a single sample.
#' * `pct_doubletons`: percentage of taxa with counts in exactly two samples.
#' * `n_samples`: total number of samples in the group.
#' * `lib_size_min`: minimum library size (total counts) across samples in the group.
#' * `lib_size_max`: maximum library size (total counts) across samples in the group.
#' * `count_mean`: average of the mean counts per sample.
#' * `count_max`: average of the max counts per sample.
#'
#' @param rec A `Recipe` object.
#'
#' @return A tibble with QC metrics.
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' test <- phy_qc(test_rec)
#' 
#' ## Check columns
#' expected_cols <- c(
#'   "var_levels", "n", "n_zero", "pct_zero", "pct_all_zero", "pct_singletons", 
#'   "pct_doubletons", "n_samples", "lib_size_min", "lib_size_max", 
#'   "count_mean", "count_max"
#' )
#' expect_equal(colnames(test), expected_cols)
#' 
#' ## Check column types
#' expect_true(all(sapply(test[,1], is.character)))
#' expect_true(all(sapply(test[,-1], is.numeric)))
#' 
#' @examples
#' data(metaHIV_phy)
#'
#' ## 1. Init Recipe
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' 
#' ## 2. Get QC metrics
#' phy_qc(rec)
phy_qc <- function(rec) {  
  check_recipe(rec)
  var_name <- get_var(rec)[[1]]
  zero_groups <- .zero_groups(rec)
  count_summary <- .count_summary(rec)
  .zero_stats(rec) %>% 
    dplyr::left_join(.zero_groups(rec), by = var_name) %>% 
    dplyr::left_join(.count_summary(rec), by = var_name) %>% 
    dplyr::rename(var_levels = !!dplyr::sym(var_name))
}

#' Pre-process OTU table for QC metrics
#' 
#' @description Pivots the count matrix to a long format, joins sample metadata, 
#' and duplicates the dataset to include an "all" grouping level.
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

#' Calculate basic zero-inflation statistics
#' 
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
  
#' Calculate prevalence-based sparsity metrics
#' 
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

#' Calculate library size and count distributions
#' 
#' @noRd
#' @keywords internal
#' @autoglobal
.count_summary <- function(rec) {
  .zero_prepro(rec) %>% 
    dplyr::group_by(!!dplyr::sym(get_var(rec)[[1]]), sample_id) %>% 
    dplyr::summarise(
      lib_size = sum(value),     
      count_mean = mean(value), 
      count_max = max(value),
      .groups = "drop_last"
    ) %>% 
    dplyr::summarise(
      n_samples = dplyr::n(),     
      lib_size_min = min(lib_size),
      lib_size_max = max(lib_size),
      count_mean = mean(count_mean),
      count_max = mean(count_max),
      .groups = "drop"
    )
}
