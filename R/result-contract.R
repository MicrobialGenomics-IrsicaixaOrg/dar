#' Collect differential-abundance results in a tidy contract
#'
#' `tidy_results()` combines the results of executed differential-abundance
#' steps without changing the engine-specific tables stored in `rec@results`.
#' Effect sizes are labelled with their engine-specific metric and must not be
#' compared directly across incompatible metrics.
#'
#' @param rec A [PrepRecipe-class] object.
#' @param steps Character vector of executed differential-abundance step IDs to
#'   collect. By default, all executed DA steps are included.
#' @param significant_only If `TRUE`, retain only rows classified as
#'   significant by their originating step.
#'
#' @return A deterministic tibble with one row per taxon, contrast and
#'   configured step. It contains taxon identifiers (`taxa_id`, `taxa`),
#'   contrast metadata (`contrast_id`, `comparison`, `contrast_type`, `var`,
#'   `numerator`, `denominator`, `at_var`, `at_level`), statistics
#'   (`effect_size`, `effect_metric`, `p_value`, `adj_p_value`, `significant`),
#'   the stable engine in `method`, and the configured instance in `step_id`.
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' results <- tidy_results(test_prep_rec)
#' results
#'
#' tidy_results(test_prep_rec, significant_only = TRUE)
tidy_results <- function(rec,
                         steps = steps_ids(rec, "da"),
                         significant_only = FALSE) {
  check_prep_recipe(rec)
  checkmate::assert_character(steps, any.missing = FALSE, unique = TRUE)
  checkmate::assert_flag(significant_only)

  if (length(steps) == 0L) {
    return(empty_tidy_results())
  }

  executed <- names(rec@results)
  unavailable <- setdiff(steps, executed)
  if (length(unavailable) > 0L) {
    skipped <- rec@execution$skipped_steps
    skipped_reasons <- character()
    if (inherits(skipped, "data.frame") && nrow(skipped) > 0L) {
      matched <- skipped[skipped$step_id %in% unavailable, , drop = FALSE]
      if (nrow(matched) > 0L) {
        skipped_reasons <- paste0(matched$step_id, ": ", matched$reason)
      }
    }
    information <- if (length(skipped_reasons) > 0L) {
      "Skipped step{?s}: {skipped_reasons}."
    } else {
      "Use {.code steps_ids(rec, 'da')} to list executed DA steps."
    }
    cli::cli_abort(
      c(
        "x" = "Requested DA step{?s} {?is/are} not available in this {.cls PrepRecipe}: {.val {unavailable}}.",
        "i" = information
      ),
      class = "dar_error_invalid_step_id"
    )
  }

  result <- purrr::map_dfr(
    steps,
    ~ standardize_da_result(rec@results[[.x]], rec, .x)
  )
  validate_tidy_result_contract(result)

  if (significant_only) {
    result <- dplyr::filter(result, .data$significant %in% TRUE)
  }

  order_tidy_results(result, rec, steps)
}

#' @noRd
tidy_result_columns <- function() {
  c(
    "taxa_id", "taxa", "contrast_id", "comparison", "contrast_type",
    "var", "numerator", "denominator", "at_var", "at_level",
    "effect_size", "effect_metric", "p_value", "adj_p_value",
    "significant", "method", "step_id"
  )
}

#' @noRd
empty_tidy_results <- function() {
  tibble::tibble(
    taxa_id = character(), taxa = character(), contrast_id = character(),
    comparison = character(), contrast_type = character(), var = character(),
    numerator = character(), denominator = character(), at_var = character(),
    at_level = character(), effect_size = numeric(), effect_metric = character(),
    p_value = numeric(), adj_p_value = numeric(), significant = logical(),
    method = character(), step_id = character()
  )
}

#' @noRd
da_result_registry <- function() {
  list(
    deseq = list(
      effect_metric = "log2_fold_change",
      p_value = c("pvalue"), adj_p_value = c("padj")
    ),
    aldex = list(
      effect_metric = "standardized_clr_effect",
      p_value = c("pval", "we.ep"), adj_p_value = c("padj", "we.eBH")
    ),
    ancom = list(
      effect_metric = "bias_corrected_log_fold_change",
      p_value = c("pval", "p"), adj_p_value = c("padj", "q")
    ),
    corncob = list(
      effect_metric = "logit_mean_abundance_coefficient",
      p_value = c("pval"), adj_p_value = c("padj")
    ),
    maaslin = list(
      effect_metric = "transformed_abundance_coefficient",
      p_value = c("pval"), adj_p_value = c("padj", "qval")
    ),
    wilcox = list(
      effect_metric = "median_difference",
      p_value = c("pval", "p"), adj_p_value = c("padj", "p.adj")
    ),
    lefse = list(
      effect_metric = "lda_score",
      p_value = c("pvalue"), adj_p_value = c("padj", "adjp")
    )
  )
}

#' @noRd
result_step <- function(rec, step_id) {
  matches <- purrr::keep(rec@steps, ~ identical(.x$id, step_id))
  if (length(matches) != 1L) {
    cli::cli_abort(
      "Could not resolve DA step {.val {step_id}} for result normalization.",
      class = "dar_error_invalid_result_contract"
    )
  }
  matches[[1L]]
}

#' @noRd
result_method <- function(step) {
  sub("^step_", "", class(step)[[1L]])
}

#' @noRd
result_effect_metric <- function(method, step) {
  registry <- da_result_registry()
  if (!method %in% names(registry)) {
    cli::cli_abort(
      "DA method {.val {method}} has no registered result contract.",
      class = "dar_error_invalid_result_contract"
    )
  }
  if (identical(method, "wilcox")) {
    norm_method <- step$norm_method %||% "compositional"
    return(paste0("median_", norm_method, "_difference"))
  }
  registry[[method]]$effect_metric
}

#' @noRd
result_numeric_column <- function(result, candidates, context) {
  selected <- candidates[candidates %in% names(result)][1L]
  if (length(selected) == 0L || is.na(selected)) {
    return(rep(NA_real_, nrow(result)))
  }
  value <- result[[selected]]
  if (!is.numeric(value)) {
    cli::cli_abort(
      "Result column {.field {selected}} for {context} must be numeric.",
      class = "dar_error_invalid_result_contract"
    )
  }
  as.numeric(value)
}

#' Read ALDEx2 coefficient p-values without modifying its native result table
#' @noRd
aldex_model_p_value <- function(result, rec, context) {
  if (is.null(get_model(rec)) || !"contrast_id" %in% names(result)) {
    return(rep(NA_real_, nrow(result)))
  }
  compiled <- compile_model(rec, "aldex")
  columns <- purrr::map_chr(
    seq_len(nrow(compiled$contrasts)),
    function(index) {
      contrast <- compiled$contrasts[index, , drop = FALSE]
      coefficient <- reparameterize_model_contrast(compiled, contrast)$coefficient
      paste0(coefficient, ":pval")
    }
  )
  names(columns) <- compiled$contrasts$contrast_id
  values <- rep(NA_real_, nrow(result))
  for (contrast_id in intersect(unique(result$contrast_id), names(columns))) {
    column <- columns[[contrast_id]]
    if (!column %in% names(result)) {
      next
    }
    if (!is.numeric(result[[column]])) {
      cli::cli_abort(
        "Result column {.field {column}} for {context} must be numeric.",
        class = "dar_error_invalid_result_contract"
      )
    }
    rows <- result$contrast_id == contrast_id
    values[rows] <- as.numeric(result[[column]][rows])
  }
  values
}

#' @noRd
result_p_value <- function(result, method, rec, candidates, context) {
  selected <- candidates[candidates %in% names(result)][1L]
  if (length(selected) > 0L && !is.na(selected)) {
    return(result_numeric_column(result, candidates, context))
  }
  if (identical(method, "aldex")) {
    return(aldex_model_p_value(result, rec, context))
  }
  rep(NA_real_, nrow(result))
}

#' @noRd
result_character_column <- function(result, name, default = NA_character_) {
  if (!name %in% names(result)) {
    return(rep(default, nrow(result)))
  }
  as.character(result[[name]])
}

#' @noRd
standardize_da_result <- function(result, rec, step_id) {
  step <- result_step(rec, step_id)
  method <- result_method(step)
  registry <- da_result_registry()
  if (!method %in% names(registry)) {
    cli::cli_abort(
      "DA method {.val {method}} has no registered result contract.",
      class = "dar_error_invalid_result_contract"
    )
  }

  result <- flatten_model_output(result)
  if (nrow(result) == 0L) {
    return(empty_tidy_results())
  }
  result <- harmonize_model_result(result, rec)

  required <- c("taxa_id", "effect", "signif")
  missing <- setdiff(required, names(result))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "x" = "DA step {.val {step_id}} returned an invalid result contract.",
        "i" = "Missing required column{?s}: {.field {missing}}."
      ),
      class = "dar_error_invalid_result_contract"
    )
  }
  if (!is.numeric(result$effect) || !is.logical(result$signif)) {
    cli::cli_abort(
      "DA step {.val {step_id}} must return numeric {.field effect} and logical {.field signif} columns.",
      class = "dar_error_invalid_result_contract"
    )
  }
  taxa_id <- as.character(result$taxa_id)
  if (anyNA(taxa_id) || any(!nzchar(taxa_id))) {
    cli::cli_abort(
      "DA step {.val {step_id}} returned missing or empty taxon identifiers.",
      class = "dar_error_invalid_result_contract"
    )
  }

  comparison <- result_character_column(result, "comparison")
  var <- result_character_column(result, "var")
  contrast_id <- result_character_column(result, "contrast_id")
  legacy <- is.na(contrast_id) | !nzchar(contrast_id)
  legacy_id <- ifelse(
    !is.na(var) & nzchar(var),
    paste0(var, "[", comparison, "]"),
    comparison
  )
  contrast_id[legacy] <- legacy_id[legacy]
  if (anyNA(contrast_id) || any(!nzchar(contrast_id))) {
    cli::cli_abort(
      "DA step {.val {step_id}} could not derive a contrast identifier.",
      class = "dar_error_invalid_result_contract"
    )
  }

  config <- registry[[method]]
  standardized <- tibble::tibble(
    taxa_id = taxa_id,
    taxa = result_character_column(result, "taxa"),
    contrast_id = contrast_id,
    comparison = comparison,
    contrast_type = result_character_column(
      result, "contrast_type", default = "legacy_pairwise"
    ),
    var = var,
    numerator = result_character_column(result, "numerator"),
    denominator = result_character_column(result, "denominator"),
    at_var = result_character_column(result, "at_var"),
    at_level = result_character_column(result, "at_level"),
    effect_size = as.numeric(result$effect),
    effect_metric = rep(result_effect_metric(method, step), nrow(result)),
    p_value = result_p_value(
      result, method, rec, config$p_value,
      paste0(method, " step ", step_id)
    ),
    adj_p_value = result_numeric_column(
      result, config$adj_p_value, paste0(method, " step ", step_id)
    ),
    significant = as.logical(result$signif),
    method = rep(method, nrow(result)),
    step_id = rep(step_id, nrow(result))
  )
  validate_tidy_result_contract(standardized)
  standardized
}

#' @noRd
validate_tidy_result_contract <- function(result) {
  expected <- tidy_result_columns()
  if (!identical(names(result), expected)) {
    cli::cli_abort(
      c(
        "x" = "The tidy DA result contract has invalid columns.",
        "i" = "Expected, in order: {.field {expected}}."
      ),
      class = "dar_error_invalid_result_contract"
    )
  }
  numeric_columns <- c("effect_size", "p_value", "adj_p_value")
  logical_columns <- "significant"
  character_columns <- setdiff(expected, c(numeric_columns, logical_columns))
  if (!all(vapply(result[character_columns], is.character, logical(1))) ||
      !all(vapply(result[numeric_columns], is.numeric, logical(1))) ||
      !is.logical(result$significant)) {
    cli::cli_abort(
      "The tidy DA result contract contains invalid column types.",
      class = "dar_error_invalid_result_contract"
    )
  }
  identifier_columns <- c(
    "taxa_id", "contrast_id", "effect_metric", "method", "step_id"
  )
  invalid_identifier <- vapply(
    result[identifier_columns],
    function(value) anyNA(value) || any(!nzchar(value)),
    logical(1)
  )
  if (any(invalid_identifier)) {
    cli::cli_abort(
      "The tidy DA result contract contains missing or empty identifiers in {.field {names(invalid_identifier)[invalid_identifier]}}.",
      class = "dar_error_invalid_result_contract"
    )
  }
  for (column in c("p_value", "adj_p_value")) {
    values <- result[[column]]
    invalid <- !is.na(values) & (!is.finite(values) | values < 0 | values > 1)
    if (any(invalid)) {
      cli::cli_abort(
        "The tidy DA result column {.field {column}} must contain values between 0 and 1 or {.val NA}.",
        class = "dar_error_invalid_result_contract"
      )
    }
  }
  keys <- c("taxa_id", "contrast_id", "step_id")
  if (anyDuplicated(result[keys])) {
    cli::cli_abort(
      "The tidy DA result contract contains duplicated taxon-contrast-step keys.",
      class = "dar_error_invalid_result_contract"
    )
  }
  invisible(TRUE)
}

#' @noRd
order_tidy_results <- function(result, rec, steps) {
  if (nrow(result) == 0L) {
    return(empty_tidy_results())
  }
  contrast_order <- rec@execution$contrasts$contrast_id
  contrast_order <- unique(c(contrast_order, result$contrast_id))
  taxa_order <- rownames(rec@phyloseq@otu_table)
  taxa_order <- unique(c(taxa_order, result$taxa_id))
  result |>
    dplyr::mutate(
      .step_order = match(.data$step_id, steps),
      .contrast_order = match(.data$contrast_id, contrast_order),
      .taxa_order = match(.data$taxa_id, taxa_order)
    ) |>
    dplyr::arrange(
      .data$.step_order, .data$.contrast_order, .data$.taxa_order,
      .data$taxa_id
    ) |>
    dplyr::select(-dplyr::starts_with("."))
}
