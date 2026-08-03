#' Recipe quality-control metrics
#'
#' `recipe_qc()` summarizes sparsity and sequencing-depth metrics from the count
#' matrix stored in a [Recipe-class] or [PrepRecipe-class] object. The summary
#' is independent of the centralized statistical model and can therefore be
#' used before model configuration and after recipe preprocessing.
#'
#' With `group_by = NULL`, the result contains one overall row. When a sample
#' metadata column is supplied, the overall row is followed by one row for each
#' observed group. The `group_by` and `group` columns are `NA` in the overall
#' row, so a real group named `"all"` remains unambiguous.
#'
#' The returned metrics are:
#'
#' * `n`: number of values in the count matrix.
#' * `n_zero`: number of zero values.
#' * `pct_zero`: percentage of zero values.
#' * `pct_all_zero`: percentage of taxa absent from every sample.
#' * `pct_singletons`: percentage of taxa detected in exactly one sample.
#' * `pct_doubletons`: percentage of taxa detected in exactly two samples.
#' * `n_samples`: number of samples.
#' * `lib_size_min`: minimum sample library size.
#' * `lib_size_max`: maximum sample library size.
#' * `count_mean`: mean of the per-sample mean counts.
#' * `count_max`: mean of the per-sample maximum counts.
#'
#' @param rec A [Recipe-class] or [PrepRecipe-class] object.
#' @param group_by `NULL` or the name of one sample metadata column used to
#'   stratify the QC summary.
#'
#' @return A tibble with grouping identifiers and QC metrics.
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' result <- recipe_qc(test_rec)
#' expect_named(result, c(
#'   "group_by", "group", "n", "n_zero", "pct_zero", "pct_all_zero",
#'   "pct_singletons", "pct_doubletons", "n_samples", "lib_size_min",
#'   "lib_size_max", "count_mean", "count_max"
#' ))
#' expect_equal(nrow(result), 1L)
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy)
#'
#' ## Overall metrics do not require a statistical model
#' recipe_qc(rec)
#'
#' ## Add explicit sample-level stratification when useful
#' recipe_qc(rec, group_by = "RiskGroup2")
recipe_qc <- function(rec, group_by = NULL) {
  check_any_recipe(rec)
  grouping <- validate_qc_group(rec, group_by)
  group_by <- grouping$name
  counts <- qc_count_data(rec)

  overall <- dplyr::bind_cols(
    tibble::tibble(group_by = NA_character_, group = NA_character_),
    qc_metrics(counts)
  )
  if (is.null(group_by)) {
    return(overall)
  }

  metadata <- grouping$metadata
  group_values <- metadata[[group_by]]
  sample_groups <- as.character(group_values)
  names(sample_groups) <- rownames(metadata)
  counts$.qc_group <- unname(sample_groups[counts$sample_id])

  grouped <- purrr::map_dfr(qc_group_levels(group_values), function(level) {
    keep <- if (is.na(level)) {
      is.na(counts$.qc_group)
    } else {
      !is.na(counts$.qc_group) & counts$.qc_group == level
    }
    dplyr::bind_cols(
      tibble::tibble(group_by = group_by, group = level),
      qc_metrics(counts[keep, , drop = FALSE])
    )
  })

  dplyr::bind_rows(overall, grouped)
}

#' Defunct recipe QC interface
#'
#' `phy_qc()` is defunct. Use [recipe_qc()] with an explicit `group_by`
#' argument when grouped metrics are required.
#'
#' @param rec A recipe object. This argument is no longer evaluated.
#' @param ... Additional arguments. These arguments are no longer evaluated.
#'
#' @return This function does not return a value; it raises a defunct error.
#' @export
#' @keywords internal
phy_qc <- function(rec, ...) {
  .Defunct(new = "recipe_qc", package = "dar")
}

#' Validate the recipe QC grouping specification
#' @noRd
validate_qc_group <- function(rec, group_by) {
  if (is.null(group_by)) {
    return(list(name = NULL, metadata = NULL))
  }
  if (!is.character(group_by) || length(group_by) != 1L ||
      is.na(group_by) || !nzchar(group_by)) {
    cli::cli_abort(
      c(
        "x" = "{.arg group_by} must be {.code NULL} or one non-empty sample metadata column name.",
        "i" = "Supply a single character value such as {.val condition}."
      ),
      class = "dar_error_invalid_qc_group"
    )
  }

  metadata <- qc_sample_metadata(rec)
  if (!group_by %in% colnames(metadata)) {
    cli::cli_abort(
      c(
        "x" = "Column {.val {group_by}} does not exist in the sample metadata.",
        "i" = "Available columns: {.val {colnames(metadata)}}."
      ),
      class = "dar_error_invalid_qc_group"
    )
  }
  values <- metadata[[group_by]]
  if (!is.atomic(values) || !is.null(dim(values))) {
    cli::cli_abort(
      "Sample metadata column {.val {group_by}} cannot be used for grouping.",
      class = "dar_error_invalid_qc_group"
    )
  }
  list(name = group_by, metadata = metadata)
}

#' Extract complete sample metadata in phyloseq sample order
#' @noRd
qc_sample_metadata <- function(rec) {
  phy <- get_phy(rec)
  metadata <- tryCatch(
    as(phyloseq::sample_data(phy), "data.frame"),
    error = function(cnd) NULL
  )
  if (is.null(metadata)) {
    return(data.frame(row.names = phyloseq::sample_names(phy)))
  }
  metadata[phyloseq::sample_names(phy), , drop = FALSE]
}

#' Reshape a recipe count matrix once for all QC summaries
#' @noRd
qc_count_data <- function(rec) {
  phy <- get_phy(rec)
  counts <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    counts <- t(counts)
  }
  counts <- counts[, phyloseq::sample_names(phy), drop = FALSE]
  tibble::as_tibble(counts, rownames = "taxa_id") |>
    tidyr::pivot_longer(
      cols = -"taxa_id", names_to = "sample_id", values_to = "value"
    )
}

#' Return deterministic observed QC groups
#' @noRd
qc_group_levels <- function(values) {
  observed <- values[!is.na(values)]
  group_levels <- if (is.factor(values)) {
    factor_levels <- levels(values)
    factor_levels[factor_levels %in% as.character(observed)]
  } else {
    as.character(sort(unique(observed), na.last = TRUE))
  }
  c(group_levels, if (anyNA(values)) NA_character_)
}

#' Calculate QC metrics for one sample cohort
#' @noRd
qc_metrics <- function(counts) {
  detections <- tapply(counts$value != 0, counts$taxa_id, sum)
  library_sizes <- tapply(counts$value, counts$sample_id, sum)
  sample_means <- tapply(counts$value, counts$sample_id, mean)
  sample_maxima <- tapply(counts$value, counts$sample_id, max)
  n <- as.numeric(nrow(counts))
  n_zero <- as.numeric(sum(counts$value == 0))

  tibble::tibble(
    n = n,
    n_zero = n_zero,
    pct_zero = n_zero / n * 100,
    pct_all_zero = mean(detections == 0) * 100,
    pct_singletons = mean(detections == 1) * 100,
    pct_doubletons = mean(detections == 2) * 100,
    n_samples = as.numeric(length(library_sizes)),
    lib_size_min = as.numeric(min(library_sizes)),
    lib_size_max = as.numeric(max(library_sizes)),
    count_mean = as.numeric(mean(sample_means)),
    count_max = as.numeric(mean(sample_maxima))
  )
}
