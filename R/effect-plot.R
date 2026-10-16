#' Plot method-specific differential-abundance effects
#'
#' Displays point estimates returned by the selected differential-abundance
#' steps for one contrast. Each configured step and effect metric has its own
#' panel and axis scales; effects are not pooled or standardized across methods.
#' Colours retain the originating step's significance decision rather than
#' recomputing it from a common p-value threshold. Confidence intervals are not
#' inferred from the common result contract.
#'
#' @param rec A [PrepRecipe-class] object.
#' @param taxa_ids Optional character vector of taxon identifiers to display.
#'   When supplied, all finite returned effects for these identifiers are
#'   included and `top_n` is not applied. Taxa with identical names remain
#'   separate units and are labelled with their identifiers when needed.
#' @param steps Character vector of executed differential-abundance step IDs.
#'   By default, all executed DA steps are included.
#' @param target Optional modeled target to include.
#' @param contrast_id Optional contrast identifier. This is required when more
#'   than one contrast remains after selecting `target`.
#' @param top_n Maximum number of taxa selected independently in each step and
#'   effect-metric panel when `taxa_ids` is `NULL`. Taxa are ranked by decreasing
#'   absolute finite effect within that panel, with taxon ID as a tie-breaker.
#'   This is a display selection, not a significance or consensus rule.
#'
#' @return An editable ggplot2 object. Returns invisible `NULL` with an
#'   informational message when no finite effects are available. Returned
#'   nonfinite effects are omitted; absent result rows are never
#'   interpreted as zero effects or nonsignificant findings.
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' results <- tidy_results(test_prep_rec)
#' effect_plt(
#'   test_prep_rec,
#'   steps = results$step_id[[1]],
#'   contrast_id = results$contrast_id[[1]],
#'   top_n = 10
#' )
effect_plt <- function(rec,
                       taxa_ids = NULL,
                       steps = steps_ids(rec, "da"),
                       target = NULL,
                       contrast_id = NULL,
                       top_n = 20) {
  check_prep_recipe(rec)
  checkmate::assert_count(top_n, positive = TRUE)
  if (!is.null(taxa_ids)) {
    checkmate::assert_character(
      taxa_ids,
      any.missing = FALSE,
      unique = TRUE,
      pattern = ".+"
    )
    unknown <- setdiff(taxa_ids, analysis_tax_table(rec)$taxa_id)
    if (length(unknown) > 0L) {
      cli::cli_abort(
        "Unknown analysis taxon identifier{?s}: {.val {unknown}}.",
        class = "dar_error_invalid_taxa_ids"
      )
    }
  }

  results <- tidy_results(rec, steps = steps)
  results <- select_effect_plot_contrast(
    results,
    rec,
    target = target,
    contrast_id = contrast_id
  )
  if (!is.null(taxa_ids)) {
    results <- dplyr::filter(results, .data$taxa_id %in% .env$taxa_ids)
  }
  results <- dplyr::filter(results, is.finite(.data$effect_size))
  if (nrow(results) == 0L) {
    cli::cli_inform("No finite effects are available for the selected hypotheses and taxa.")
    return(invisible(NULL))
  }

  results$.step_order <- match(results$step_id, steps)
  results <- results |>
    dplyr::group_by(.data$step_id, .data$effect_metric) |>
    dplyr::arrange(.data$.step_order, dplyr::desc(abs(.data$effect_size)), .data$taxa_id)
  if (is.null(taxa_ids)) { results <- dplyr::slice_head(results, n = top_n) }
  results <- dplyr::ungroup(results) |>
    dplyr::arrange(
      .data$.step_order,
      dplyr::desc(abs(.data$effect_size)),
      .data$taxa_id
    ) |>
    dplyr::select(-".step_order")
  panel_names <- paste(results$step_id, results$effect_metric, sep = "\n")
  results$panel <- factor(panel_names, levels = unique(panel_names))
  results$decision <- factor(
    ifelse(
      is.na(results$significant),
      "Not available",
      ifelse(results$significant, "Significant", "Not significant")
    ),
    levels = c("Significant", "Not significant", "Not available")
  )
  labels <- effect_plot_taxon_labels(results)
  results$.plot_taxon <- factor(
    paste(as.integer(results$panel), results$taxa_id, sep = "\r"),
    levels = rev(unique(paste(as.integer(results$panel), results$taxa_id, sep = "\r")))
  )
  plot_labels <- labels[results$taxa_id]
  names(plot_labels) <- as.character(results$.plot_taxon)

  ggplot2::ggplot(
    results,
    ggplot2::aes(x = .data$effect_size, y = .data$.plot_taxon, colour = .data$decision)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free") +
    ggplot2::scale_y_discrete(labels = plot_labels) +
    ggplot2::scale_colour_manual(
      values = c(
        "Significant" = "#0072B2", "Not significant" = "grey55",
        "Not available" = "#E69F00"
      ),
      drop = TRUE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Effect estimate (panel-specific metric)",
      y = NULL,
      colour = "Method decision",
      title = unique(results$contrast_id),
      caption = "Point estimates on method-specific scales; no confidence intervals are inferred."
    )
}

#' @noRd
select_effect_plot_contrast <- function(results, rec, target = NULL, contrast_id = NULL) {
  if (!is.null(get_model(rec))) {
    plan <- resolve_model(rec)$contrast_plan
    results <- filter_model_hypotheses(
      results,
      rec,
      target = target,
      contrast_id = contrast_id
    )
  } else {
    plan <- dplyr::distinct(results, .data$contrast_id, .data$var)
    available_targets <- unique(plan$var[!is.na(plan$var)])
    if (!is.null(target)) {
      check_effect_plot_selector(target, "target", available_targets)
      results <- dplyr::filter(results, .data$var %in% .env$target)
    }
  }
  if (!is.null(target)) { plan <- dplyr::filter(plan, .data$var %in% .env$target) }
  if (!is.null(contrast_id)) {
    check_effect_plot_selector(contrast_id, "contrast_id", plan$contrast_id)
    plan <- dplyr::filter(plan, .data$contrast_id %in% .env$contrast_id)
    results <- dplyr::filter(results, .data$contrast_id %in% .env$contrast_id)
  }
  if (dplyr::n_distinct(plan$contrast_id) > 1L) {
    cli::cli_abort(
      "{.arg contrast_id} is required when multiple contrasts are available.",
      class = "dar_error_ambiguous_plot_contrast"
    )
  }
  results
}

#' @noRd
check_effect_plot_selector <- function(value, argument, available) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || !value %in% available) {
    cli::cli_abort(
      "Invalid {.arg {argument}}. Available values: {.val {available}}.",
      class = "dar_error_invalid_contrast_selector"
    )
  }
  invisible(value)
}

#' @noRd
effect_plot_taxon_labels <- function(results) {
  taxa <- dplyr::distinct(results, .data$taxa_id, .data$taxa)
  labels <- taxa$taxa
  labels[is.na(labels) | !nzchar(labels)] <- taxa$taxa_id[
    is.na(labels) | !nzchar(labels)
  ]
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_labels] <- paste0(
    labels[duplicated_labels], " [", taxa$taxa_id[duplicated_labels], "]"
  )
  stats::setNames(labels, taxa$taxa_id)
}
