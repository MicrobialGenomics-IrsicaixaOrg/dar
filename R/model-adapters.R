# ENGINE ADAPTERS FOR CENTRAL MODELS -----------------------------------------

#' @noRd
model_phyloseq <- function(rec, rarefy = FALSE) {
  phy <- get_phy(rec) |>
    use_rarefy(rarefy) |>
    phyloseq::tax_glom(taxrank = recipe_tax_level(rec), NArm = FALSE)
  metadata <- resolve_model(rec)$data |>
    tibble::column_to_rownames("sample_id") |>
    data.frame(check.names = FALSE) |>
    phyloseq::sample_data()
  phyloseq::sample_data(phy) <- metadata
  phy
}

#' Reparameterize a fixed model so a planned contrast is one coefficient
#' @noRd
reparameterize_model_contrast <- function(compiled, contrast) {
  data <- compiled$data
  resolved <- compiled$resolved
  type <- contrast$contrast_type[[1]]
  target <- contrast$var[[1]]

  relevel_var <- function(var, ref) {
    data[[var]] <<- stats::relevel(factor(data[[var]]), ref = ref)
  }

  if (identical(type, "main")) {
    relevel_var(target, contrast$denominator[[1]])
  } else if (identical(type, "condition_within_time")) {
    relevel_var(target, contrast$denominator[[1]])
    relevel_var(resolved$time, contrast$at_level[[1]])
  } else if (identical(type, "time_within_condition")) {
    relevel_var(resolved$time, contrast$denominator[[1]])
    relevel_var(contrast$at_var[[1]], contrast$at_level[[1]])
  } else if (identical(type, "difference_in_differences")) {
    numerator <- strsplit(contrast$numerator[[1]], ":", fixed = TRUE)[[1]]
    denominator <- strsplit(contrast$denominator[[1]], ":", fixed = TRUE)[[1]]
    relevel_var(target, denominator[[1]])
    relevel_var(resolved$time, denominator[[2]])
    if (length(numerator) != 2L) {
      cli::cli_abort(
        "Invalid difference-in-differences contrast specification.",
        class = "dar_error_invalid_result_contract"
      )
    }
  }

  new_matrix <- stats::model.matrix(resolved$fixed_formula, data = data)
  transformation <- qr.solve(compiled$matrix, new_matrix)
  new_weights <- as.numeric(t(transformation) %*% contrast$weights[[1]])
  active <- which(abs(new_weights) > sqrt(.Machine$double.eps))
  if (length(active) != 1L || !isTRUE(all.equal(abs(new_weights[[active]]), 1))) {
    cli::cli_abort(
      c(
        "x" = "Contrast {.val {contrast$contrast_id}} could not be represented by one model coefficient.",
        "i" = "The engine requires coefficient-wise reparameterization."
      ),
      class = "dar_error_invalid_result_contract"
    )
  }
  list(
    data = data,
    matrix = new_matrix,
    coefficient = colnames(new_matrix)[[active]],
    coefficient_index = active,
    sign = sign(new_weights[[active]])
  )
}

#' @noRd
bind_model_result <- function(result, contrast, rec, effect, padj, signif) {
  result <- tibble::as_tibble(result, rownames = "taxa_id") |>
    dplyr::left_join(tax_table(rec), by = "taxa_id")
  result$contrast_id <- contrast$contrast_id[[1]]
  result$comparison <- contrast$comparison[[1]]
  result$contrast_type <- contrast$contrast_type[[1]]
  result$var <- contrast$var[[1]]
  result$effect <- effect
  result$padj <- padj
  result$signif <- signif
  result
}

# DESEQ2 ----------------------------------------------------------------------

#' @noRd
run_deseq_model <- function(rec, test, fitType, betaPrior, type,
                            max_significance, log2FC, rarefy,
                            engine_args = list()) {
  compiled <- compile_model(rec, "deseq")
  phy <- model_phyloseq(rec, rarefy)
  result_stage <- if (identical(type, "normal") &&
                      length(compiled$resolved$interaction_targets) > 0L) {
    "results"
  } else {
    "shrink"
  }
  check_unused_engine_args(
    "deseq", engine_args, c("size_factors", "fit", result_stage)
  )

  fit_dds <- function(phy, design) {
    dds <- suppressMessages(phyloseq::phyloseq_to_deseq2(phy, design = design))
    geo_means <- apply(DESeq2::counts(dds), 1, function(x) {
      positive <- x[x > 0]
      if (length(positive) == 0L) 1 else exp(mean(log(positive)))
    })
    dds <- exec_engine_stage(
      "deseq", "size_factors", engine_args,
      fixed = list(object = dds, geoMeans = geo_means)
    )
    tryCatch(
      exec_engine_stage(
        "deseq", "fit", engine_args,
        fixed = list(
          object = dds, fitType = fitType, test = test,
          betaPrior = betaPrior
        ),
        defaults = list(quiet = TRUE)
      ),
      error = function(cnd) {
        if (!grepl("all gene-wise dispersion estimates", conditionMessage(cnd),
                   fixed = TRUE) || !identical(test, "Wald")) {
          stop(cnd)
        }
        dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
        DESeq2::dispersions(dds) <-
          SummarizedExperiment::rowData(dds)$dispGeneEst
        DESeq2::nbinomWaldTest(
          dds, betaPrior = betaPrior, modelMatrixType = "standard",
          quiet = TRUE
        )
      }
    )
  }

  out <- purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    parameterized <- reparameterize_model_contrast(compiled, contrast)
    contrast_phy <- phy
    phyloseq::sample_data(contrast_phy) <- parameterized$data |>
      tibble::column_to_rownames("sample_id") |>
      data.frame(check.names = FALSE) |>
      phyloseq::sample_data()
    contrast_dds <- fit_dds(contrast_phy, compiled$formula)
    result_names <- DESeq2::resultsNames(contrast_dds)
    normalize <- function(x) gsub("[^[:alnum:]]", "", tolower(x))
    coefficient <- normalize(parameterized$coefficient)
    normalized_names <- normalize(result_names)
    coefficient_index <- which(normalized_names == coefficient)
    if (length(coefficient_index) == 0L) {
      coefficient_index <- which(startsWith(normalized_names, paste0(coefficient, "vs")))
    }
    if (length(coefficient_index) == 0L) {
      candidates <- which(startsWith(normalized_names, coefficient))
      if (length(candidates) == 1L) {
        coefficient_index <- candidates
      }
    }
    if (length(coefficient_index) != 1L) {
      cli::cli_abort(
        "DESeq2 did not expose coefficient {.val {parameterized$coefficient}} for contrast {.val {contrast$contrast_id}}.",
        class = "dar_error_invalid_result_contract"
      )
    }
    res <- if (identical(result_stage, "results")) {
      exec_engine_stage(
        "deseq", "results", engine_args,
        fixed = list(
          object = contrast_dds,
          name = result_names[[coefficient_index]]
        )
      )
    } else {
      exec_engine_stage(
        "deseq", "shrink", engine_args,
        fixed = list(
          dds = contrast_dds, coef = coefficient_index, type = type
        ),
        defaults = list(quiet = TRUE)
      )
    }
    effect_sign <- parameterized$sign

    res_tbl <- tibble::as_tibble(res, rownames = "taxa_id") |>
      dplyr::left_join(tax_table(rec), by = "taxa_id") |>
      dplyr::mutate(
        contrast_id = contrast$contrast_id[[1]],
        comparison = contrast$comparison[[1]],
        contrast_type = contrast$contrast_type[[1]],
        var = contrast$var[[1]],
        effect = .data$log2FoldChange * effect_sign,
        signif = !is.na(.data$padj) & .data$padj < max_significance &
          abs(.data$effect) >= log2FC
      )
    res_tbl
  })
  list(model = out)
}

# ALDEX2 ----------------------------------------------------------------------

#' @noRd
run_aldex_model <- function(rec, max_significance, mc.samples, denom, rarefy,
                            engine_args = list()) {
  if (!identical(denom, "all") && !is.numeric(denom)) {
    cli::cli_abort(
      c(
        "x" = "ALDEx2 model matrices require {.arg denom = 'all'} or explicit offsets.",
        "i" = "The requested denominator cannot be shared safely across model contrasts."
      ),
      class = "dar_error_invalid_model"
    )
  }
  compiled <- compile_model(rec, "aldex")
  phy <- model_phyloseq(rec, rarefy)
  check_unused_engine_args("aldex", engine_args, c("clr", "fit", "effect"))
  reads <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    reads <- t(reads)
  }
  reads <- reads[, compiled$data$sample_id, drop = FALSE]

  out <- purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    parameterized <- reparameterize_model_contrast(compiled, contrast)
    clr <- suppressMessages(exec_engine_stage(
      "aldex", "clr", engine_args,
      fixed = list(
        reads = reads, conds = parameterized$matrix,
        denom = denom, mc.samples = mc.samples
      ),
      defaults = list(verbose = FALSE)
    ))
    glm <- exec_engine_stage(
      "aldex", "fit", engine_args, fixed = list(clr = clr),
      defaults = list(verbose = FALSE, fdr.method = "BH")
    )
    effects <- exec_engine_stage(
      "aldex", "effect", engine_args, fixed = list(clr = clr),
      defaults = list(verbose = FALSE, CI = TRUE)
    )
    coefficient <- parameterized$coefficient
    estimate_col <- paste0(coefficient, ":Est")
    padj_col <- paste0(coefficient, ":pval.padj")
    if (!all(c(estimate_col, padj_col) %in% names(glm)) ||
        !coefficient %in% names(effects)) {
      cli::cli_abort(
        "ALDEx2 did not return coefficient {.val {coefficient}}.",
        class = "dar_error_invalid_result_contract"
      )
    }
    effect <- effects[[coefficient]]$effect * parameterized$sign
    padj <- glm[[padj_col]]
    bind_model_result(
      glm,
      contrast,
      rec,
      effect = effect,
      padj = padj,
      signif = !is.na(padj) & padj < max_significance
    )
  })
  list(model = out)
}

# ANCOM-BC2 -------------------------------------------------------------------

#' @noRd
match_ancom_column <- function(columns, prefix, coefficient) {
  normalize <- function(x) gsub("[^[:alnum:]]", "", tolower(x))
  candidates <- columns[startsWith(columns, paste0(prefix, "_"))]
  suffix <- sub(paste0("^", prefix, "_"), "", candidates)
  hit <- candidates[normalize(suffix) == normalize(coefficient)]
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

#' @noRd
run_ancom_model <- function(rec, p_adj_method, prv_cut, lib_cut, s0_perc,
                            struc_zero, neg_lb, alpha, n_cl, verbose,
                            rarefy, engine_args = list()) {
  compiled <- compile_model(rec, "ancom")
  phy <- model_phyloseq(rec, rarefy)
  check_unused_engine_args("ancom", engine_args, "fit")

  out <- purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    parameterized <- reparameterize_model_contrast(compiled, contrast)
    contrast_phy <- phy
    phyloseq::sample_data(contrast_phy) <- parameterized$data |>
      tibble::column_to_rownames("sample_id") |>
      data.frame(check.names = FALSE) |>
      phyloseq::sample_data()
    fit <- exec_engine_stage(
      "ancom", "fit", engine_args,
      fixed = list(
        data = mia::convertFromPhyloseq(contrast_phy),
        tax_level = recipe_tax_level(rec),
        fix_formula = compiled$fix_formula,
        rand_formula = compiled$rand_formula,
        p_adj_method = p_adj_method,
        prv_cut = prv_cut,
        lib_cut = lib_cut,
        s0_perc = s0_perc,
        group = NULL,
        struc_zero = struc_zero,
        neg_lb = neg_lb,
        alpha = alpha,
        n_cl = n_cl,
        verbose = verbose,
        global = FALSE,
        pairwise = FALSE,
        dunnet = FALSE,
        trend = FALSE
      )
    )
    stats <- fit$res
    lfc_col <- match_ancom_column(names(stats), "lfc", parameterized$coefficient)
    q_col <- match_ancom_column(names(stats), "q", parameterized$coefficient)
    diff_col <- match_ancom_column(names(stats), "diff", parameterized$coefficient)
    if (anyNA(c(lfc_col, q_col))) {
      cli::cli_abort(
        "ANCOM-BC2 did not return coefficient {.val {parameterized$coefficient}}.",
        class = "dar_error_invalid_result_contract"
      )
    }
    tibble::as_tibble(stats) |>
      dplyr::transmute(
        taxa_id = .data$taxon,
        padj = .data[[q_col]],
        effect = .data[[lfc_col]] * parameterized$sign,
        signif = if (!is.na(diff_col)) .data[[diff_col]] else
          !is.na(.data[[q_col]]) & .data[[q_col]] < alpha,
        contrast_id = contrast$contrast_id[[1]],
        comparison = contrast$comparison[[1]],
        contrast_type = contrast$contrast_type[[1]],
        var = contrast$var[[1]]
      ) |>
      dplyr::left_join(tax_table(rec), by = "taxa_id") |>
      dplyr::relocate("taxa_id", "taxa")
  })
  list(model = out)
}

# CORNCOB ---------------------------------------------------------------------

#' @noRd
run_corncob_model <- function(rec, phi.formula, link, phi.link,
                              filter_discriminant, fdr_cutoff, fdr,
                              log2FC, rarefy, engine_args = list()) {
  compiled <- compile_model(rec, "corncob")
  phy <- model_phyloseq(rec, rarefy)
  check_unused_engine_args("corncob", engine_args, "fit")
  reads <- as(phyloseq::otu_table(phy), "matrix")
  if (phyloseq::taxa_are_rows(phy)) {
    reads <- t(reads)
  }
  reads <- reads[compiled$data$sample_id, , drop = FALSE]
  metadata <- data.frame(compiled$data, row.names = compiled$data$sample_id)
  total <- rowSums(reads)

  per_taxon <- purrr::map_dfr(seq_len(ncol(reads)), function(taxon_index) {
    fit_data <- metadata
    fit_data$W <- reads[, taxon_index]
    fit_data$M <- total
    fit_formula <- stats::update(compiled$formula, cbind(W, M - W) ~ .)
    fit <- tryCatch(
      exec_engine_stage(
        "corncob", "fit", engine_args,
        fixed = list(
          formula = fit_formula, phi.formula = phi.formula,
          data = fit_data, link = link, phi.link = phi.link
        )
      ),
      error = function(cnd) NULL
    )
    if (is.null(fit)) {
      return(tibble::tibble())
    }
    covariance <- tryCatch(
      {
        hessian <- utils::getFromNamespace("hessian", "corncob")(fit)
        solve(hessian)[seq_len(fit$np.mu), seq_len(fit$np.mu), drop = FALSE]
      },
      error = function(cnd) NULL
    )
    if (is.null(covariance)) {
      return(tibble::tibble())
    }
    purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
      contrast <- compiled$contrasts[index, , drop = FALSE]
      weights <- contrast$weights[[1]]
      estimate <- sum(weights * fit$b.mu)
      standard_error <- sqrt(drop(t(weights) %*% covariance %*% weights))
      pvalue <- if (is.finite(standard_error) && standard_error > 0) {
        2 * stats::pnorm(abs(estimate / standard_error), lower.tail = FALSE)
      } else {
        NA_real_
      }
      tibble::tibble(
        taxa_id = colnames(reads)[[taxon_index]],
        contrast_id = contrast$contrast_id[[1]],
        comparison = contrast$comparison[[1]],
        contrast_type = contrast$contrast_type[[1]],
        var = contrast$var[[1]],
        effect = estimate,
        stderr = standard_error,
        pval = pvalue,
        separated = isTRUE(fit$sep_da)
      )
    })
  })

  if (nrow(per_taxon) == 0L) {
    cli::cli_abort(
      "corncob could not fit the centralized model for any taxon.",
      class = "dar_error_model_engine"
    )
  }

  out <- per_taxon |>
    dplyr::group_by(.data$contrast_id) |>
    dplyr::mutate(
      padj = stats::p.adjust(.data$pval, method = fdr),
      padj = ifelse(filter_discriminant & .data$separated, NA_real_, .data$padj),
      signif = !is.na(.data$padj) & .data$padj < fdr_cutoff &
        abs(.data$effect) >= log2FC
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(tax_table(rec), by = "taxa_id") |>
    dplyr::relocate("taxa_id", "taxa")
  list(model = out)
}

# MAASLIN3 --------------------------------------------------------------------

#' @noRd
run_maaslin_model <- function(rec, min_abundance, min_prevalence, min_variance,
                              normalization, transform, max_significance,
                              correction, standardize,
                              median_comparison_abundance, rarefy,
                              engine_args = list()) {
  compiled <- compile_model(rec, "maaslin")
  phy <- model_phyloseq(rec, rarefy)
  check_unused_engine_args("maaslin", engine_args, c("fit", "contrast"))
  input_metadata <- data.frame(compiled$data, row.names = compiled$data$sample_id)
  input_metadata$sample_id <- NULL
  input_data <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    input_data <- t(input_data)
  }
  input_data <- input_data[, rownames(input_metadata), drop = FALSE]
  output_dir <- tempfile("dar_maaslin3_")
  dir.create(output_dir)

  fit <- purrr::quietly(function() {
    exec_engine_stage(
      "maaslin", "fit", engine_args,
      fixed = list(
        input_data = input_data,
        input_metadata = input_metadata,
        output = output_dir,
        formula = compiled$formula_text,
        min_abundance = min_abundance,
        min_prevalence = min_prevalence,
        min_variance = min_variance,
        normalization = normalization,
        transform = transform,
        max_significance = max_significance,
        correction = correction,
        standardize = standardize,
        median_comparison_abundance = median_comparison_abundance
      ),
      defaults = list(
        verbosity = "ERROR", plot_summary_plot = FALSE,
        plot_associations = FALSE, max_pngs = 0
      )
    )
  })()$result

  contrast_matrix <- do.call(rbind, compiled$contrasts$weights)
  colnames(contrast_matrix) <- colnames(compiled$matrix)
  rownames(contrast_matrix) <- compiled$contrasts$contrast_id
  tested <- exec_engine_stage(
    "maaslin", "contrast", engine_args,
    fixed = list(
      maaslin3_fit = fit,
      contrast_mat = contrast_matrix,
      max_significance = max_significance,
      correction = correction,
      median_comparison_abundance = median_comparison_abundance,
      evaluate_only = "abundance"
    )
  )
  out <- tested$fit_data_abundance$results |>
    tibble::as_tibble() |>
    dplyr::transmute(
      taxa_id = .data$feature,
      contrast_id = as.character(.data$test),
      effect = .data$coef,
      stderr = .data$stderr,
      pval = .data$pval_individual,
      padj = .data$qval_individual,
      signif = !is.na(.data$qval_individual) & .data$qval_individual < max_significance
    ) |>
    dplyr::left_join(
      dplyr::select(
        compiled$contrasts, "contrast_id", "comparison", "contrast_type", "var"
      ),
      by = "contrast_id"
    ) |>
    dplyr::left_join(tax_table(rec), by = "taxa_id") |>
    dplyr::relocate("taxa_id", "taxa")
  list(model = out)
}

# SIMPLE-MODEL ENGINES --------------------------------------------------------

#' @noRd
run_wilcox_model <- function(rec, norm_method, max_significance,
                             p_adj_method, rarefy, engine_args = list()) {
  compiled <- compile_model(rec, "wilcox")
  phy <- model_phyloseq(rec, rarefy)
  check_unused_engine_args("wilcox", engine_args, "test")
  transformed <- phy |>
    microbiome::transform(transform = norm_method) |>
    phyloseq::otu_table() |>
    as("matrix")
  if (phyloseq::taxa_are_rows(phy)) {
    transformed <- t(transformed)
  }
  data <- compiled$data

  out <- purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    target <- contrast$var[[1]]
    numerator <- contrast$numerator[[1]]
    denominator <- contrast$denominator[[1]]
    keep <- data[[target]] %in% c(numerator, denominator)
    pvalues <- apply(transformed[keep, , drop = FALSE], 2, function(values) {
      groups <- data[[target]][keep]
      suppressWarnings(exec_engine_stage(
        "wilcox", "test", engine_args,
        fixed = list(
          x = values[groups == numerator],
          y = values[groups == denominator]
        ),
        defaults = list(exact = FALSE)
      )$p.value)
    })
    effects <- apply(transformed[keep, , drop = FALSE], 2, function(values) {
      groups <- data[[target]][keep]
      stats::median(values[groups == numerator]) -
        stats::median(values[groups == denominator])
    })
    padj <- stats::p.adjust(pvalues, method = p_adj_method)
    tibble::tibble(
      taxa_id = colnames(transformed),
      contrast_id = contrast$contrast_id[[1]],
      comparison = contrast$comparison[[1]],
      contrast_type = contrast$contrast_type[[1]],
      var = target,
      effect = as.numeric(effects),
      pval = as.numeric(pvalues),
      padj = as.numeric(padj),
      signif = !is.na(padj) & padj < max_significance
    )
  }) |>
    dplyr::left_join(tax_table(rec), by = "taxa_id") |>
    dplyr::relocate("taxa_id", "taxa")
  list(model = out)
}

#' @noRd
run_lefse_model <- function(rec, kruskal.threshold, wilcox.threshold,
                            lda.threshold, assay, trim.names, rarefy,
                            engine_args = list()) {
  compiled <- compile_model(rec, "lefse")
  lefse_mat <- prepro_lefse(rec, rarefy)
  check_unused_engine_args("lefse", engine_args, c("transform", "fit"))
  out <- purrr::map_dfr(seq_len(nrow(compiled$contrasts)), function(index) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    target <- contrast$var[[1]]
    numerator <- contrast$numerator[[1]]
    denominator <- contrast$denominator[[1]]
    metadata <- compiled$data %>%
      dplyr::filter(.data[[target]] %in% c(.env$denominator, .env$numerator))
    metadata[[target]] <- factor(
      metadata[[target]],
      levels = c(denominator, numerator)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = lefse_mat[, metadata$sample_id, drop = FALSE]),
      colData = data.frame(metadata, row.names = metadata$sample_id)
    )
    se <- exec_engine_stage(
      "lefse", "transform", engine_args,
      fixed = list(se = se)
    )
    raw <- exec_engine_stage(
      "lefse", "fit", engine_args,
      fixed = list(
        relab = se,
        classCol = target,
        kruskal.threshold = 1,
        wilcox.threshold = 1,
        lda.threshold = 0,
        subclassCol = NULL,
        assay = assay,
        trim.names = trim.names
      )
    ) |>
      tibble::as_tibble() |>
      dplyr::rename(lefse_id = "features") |>
      dplyr::mutate(lefse_id = stringr::str_remove_all(.data$lefse_id, "`"))
    tests <- kruskal_test(se, metadata[[target]])
    abundance <- SummarizedExperiment::assay(se, i = assay)
    direction <- apply(abundance, 1, function(values) {
      sign(stats::median(values[metadata[[target]] == numerator]) -
           stats::median(values[metadata[[target]] == denominator]))
    })

    raw |>
      dplyr::left_join(tests, by = c("lefse_id" = "otu")) |>
      dplyr::mutate(
        taxa = stringr::str_remove_all(.data$lefse_id, ".*[|]"),
        effect = abs(.data$scores) * direction[.data$lefse_id],
        padj = .data$adjp,
        signif = !is.na(.data$padj) & .data$padj < kruskal.threshold &
          abs(.data$effect) >= lda.threshold,
        contrast_id = contrast$contrast_id[[1]],
        comparison = contrast$comparison[[1]],
        contrast_type = contrast$contrast_type[[1]],
        var = target
      ) |>
      dplyr::left_join(tax_table(rec), by = "taxa") |>
      dplyr::relocate("taxa_id", "taxa")
  })
  list(model = out)
}
