#' Extract abundances in the shared taxonomic analysis view
#'
#' `tidy_abundance()` returns one row per taxon and sample. When a centralized
#' model is defined, taxa are aggregated through its selected taxonomic rank,
#' using the same identities as differential-abundance results. Without a
#' centralized model, the stored features are returned without aggregation.
#'
#' This is a descriptive view of the data currently stored in the recipe. A
#' prepared recipe includes common preprocessing and the shared model cohort;
#' an unprepared recipe does not execute its configured steps. Neither view
#' reproduces method-specific filtering, rarefaction or transformations, or
#' adjusted model predictions.
#'
#' Taxonomic columns are prefixed with `taxonomy_` and sample metadata columns
#' with `metadata_`, so annotations cannot overwrite identifiers or abundances.
#' Full lineage through the selected rank is retained; lower ranks removed by
#' aggregation remain missing. Taxa sharing a label but not a lineage retain
#' distinct `taxa_id` values.
#'
#' @param rec A [Recipe-class] or [PrepRecipe-class] object.
#' @param taxa_ids Optional character vector of taxon identifiers to retain.
#'   IDs refer to the shared analysis view, not absorbed features. `NULL`
#'   retains every analysis taxon; `character(0)` returns no rows. Output keeps
#'   analysis taxon order rather than selector order.
#' @param transform Either `"identity"` for stored abundances or
#'   `"compositional"` for within-sample fractions. Fractions use all shared
#'   analysis taxa as their denominator before `taxa_ids` selection. Zero-depth
#'   samples have missing fractions and produce a warning.
#' @param metadata Character vector of sample metadata columns to include.
#'   `NULL` includes every column; `character(0)` includes none. Supply the
#'   original, unprefixed column names.
#'
#' @return A tibble with `taxa_id`, `taxa`, `sample_id`, `abundance`, complete
#'   namespaced taxonomy and selected namespaced sample metadata. `taxa` uses
#'   the configured taxonomic label, or feature IDs when no single rank is
#'   configured. Sample metadata types are preserved. Rows follow analysis
#'   taxon order, then sample order.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus")
#' abundances <- tidy_abundance(rec, metadata = "RiskGroup2")
#' head(abundances)
#'
#' fractions <- tidy_abundance(
#'   rec, transform = "compositional", metadata = character(0)
#' )
tidy_abundance <- function(rec, taxa_ids = NULL, transform = "identity", metadata = NULL) {
  check_any_recipe(rec)
  if (!is.null(taxa_ids)) {
    checkmate::assert_character(taxa_ids, any.missing = FALSE, unique = TRUE)
  }
  if (!is.null(metadata)) {
    checkmate::assert_character(metadata, any.missing = FALSE, unique = TRUE)
  }
  transform <- match.arg(transform, c("identity", "compositional"))

  phy <- analysis_phyloseq(rec)
  counts <- recipe_count_matrix(rec, phy = phy)
  available_ids <- rownames(counts)
  unavailable <- setdiff(taxa_ids, available_ids)
  if (length(unavailable) > 0L) {
    cli::cli_abort(
      "Taxon ID{?s} {.val {unavailable}} {?is/are} outside the shared analysis view.",
      class = "dar_error_invalid_taxa_id"
    )
  }

  if (identical(transform, "compositional")) {
    if (any(!is.finite(counts)) || any(counts < 0)) {
      cli::cli_abort(
        "Compositional transformation requires finite, non-negative abundances.",
        class = "dar_error_invalid_input"
      )
    }
    library_sizes <- colSums(counts)
    zero_samples <- colnames(counts)[library_sizes == 0]
    counts <- sweep(counts, 2L, library_sizes, "/")
    if (length(zero_samples) > 0L) {
      counts[, zero_samples] <- NA_real_
      cli::cli_warn(
        "Compositional abundances are missing for zero-depth sample{?s}: {.val {zero_samples}}.",
        class = "dar_warning_zero_library"
      )
    }
  }

  if (!is.null(taxa_ids)) {
    counts <- counts[available_ids %in% taxa_ids, , drop = FALSE]
  }
  result <- tibble::tibble(
    taxa_id = rep(as.character(rownames(counts)), each = ncol(counts)),
    sample_id = rep(colnames(counts), times = nrow(counts)),
    abundance = as.numeric(t(counts))
  )

  taxonomy <- recipe_taxonomy_data(rec, phy = phy)
  tax_level <- recipe_tax_level(rec)
  result$taxa <- if (length(tax_level) == 1L &&
                     tax_level %in% names(taxonomy)) {
    taxonomy[[tax_level]][match(result$taxa_id, taxonomy$taxa_id)]
  } else {
    result$taxa_id
  }
  result <- result[c("taxa_id", "taxa", "sample_id", "abundance")]
  taxonomy <- taxonomy[
    match(result$taxa_id, taxonomy$taxa_id),
    setdiff(names(taxonomy), "taxa_id"),
    drop = FALSE
  ]
  if (ncol(taxonomy) > 0L) {
    names(taxonomy) <- paste0("taxonomy_", names(taxonomy))
  }

  sample_metadata <- tryCatch(
    as(phyloseq::sample_data(phy), "data.frame"),
    error = function(cnd) NULL
  )
  metadata_names <- if (is.null(metadata)) names(sample_metadata) else metadata
  unavailable <- setdiff(metadata_names, names(sample_metadata))
  if (length(unavailable) > 0L) {
    cli::cli_abort(
      "Sample metadata column{?s} {.val {unavailable}} {?does/do} not exist.",
      class = "dar_error_invalid_input"
    )
  }
  if (length(metadata_names) == 0L) {
    return(dplyr::bind_cols(result, taxonomy))
  }
  sample_metadata <- sample_metadata[
    match(result$sample_id, rownames(sample_metadata)),
    metadata_names,
    drop = FALSE
  ]
  names(sample_metadata) <- paste0("metadata_", names(sample_metadata))
  dplyr::bind_cols(result, taxonomy, tibble::as_tibble(sample_metadata))
}

#' Inspect the recorded execution of a prepared recipe
#'
#' `get_execution()` exposes the manifest recorded by [prep()] without
#' recomputing contrasts, modifying results or rerunning any analysis. Samples
#' listed in `dropped_samples` were removed by the centralized model's
#' missing-value policy, not by method-specific filtering.
#'
#' Older prepared recipes may lack all or part of this manifest. Unavailable
#' fields are returned as `NULL`: missing provenance does not imply that no
#' samples were removed or no steps were skipped.
#'
#' @param rec A [PrepRecipe-class] object.
#'
#' @return A list containing the recorded `contrasts`, `dropped_samples`,
#'   `executed_steps` and `skipped_steps`. The contrasts and skipped-step tables
#'   retain their existing columns, including incompatibility reasons. No
#'   centralized model also yields `NULL` contrasts in current recipes.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' prepared <- recipe(metaHIV_phy) |>
#'   step_filter_by_prevalence(threshold = 0.1) |>
#'   prep(parallel = FALSE)
#' get_execution(prepared)
get_execution <- function(rec) {
  check_prep_recipe(rec)
  fields <- c("contrasts", "dropped_samples", "executed_steps", "skipped_steps")
  stats::setNames(lapply(fields, function(field) rec@execution[[field]]), fields)
}
