make_maaslin_contrast_fits <- function(partial = FALSE, mixed = FALSE) {
  withr::local_seed(130)
  metadata <- data.frame(
    condition = factor(rep(c("A", "B", "C"), each = 12)),
    subject = factor(rep(seq_len(18), each = 2)),
    visit = rep(0:1, 18)
  )
  fits <- lapply(seq_len(8), function(index) {
    metadata$value <- 3 + index / 4 +
      c(0, index / 5, -index / 7)[metadata$condition] +
      sin(as.numeric(metadata$subject)) + stats::rnorm(nrow(metadata), sd = 0.4)
    keep <- if (partial && index == 1L) metadata$condition != "C" else TRUE
    if (mixed) {
      lmerTest::lmer(value ~ condition + visit + (1 | subject),
                     data = metadata[keep, , drop = FALSE])
    } else {
      stats::lm(value ~ condition, data = metadata[keep, , drop = FALSE])
    }
  })
  names(fits) <- paste0("taxon_", seq_along(fits))
  list(fit_data_abundance = list(
    fits = fits,
    results = data.frame(
      feature = names(fits), N = nrow(metadata),
      N_not_zero = if (partial) c(24L, rep(36L, 7)) else 36L
    )
  ))
}

maaslin_contrast_matrix <- function(mixed = FALSE) {
  weights <- rbind(B_A = c(0, 1, 0), C_A = c(0, 0, 1), C_B = c(0, -1, 1))
  colnames(weights) <- c("(Intercept)", "conditionB", "conditionC")
  if (mixed) {
    weights <- cbind(weights, visit = 0)
  }
  weights
}

sort_maaslin_contrasts <- function(results) {
  results <- as.data.frame(results)
  results <- results[order(results$feature, results$test), , drop = FALSE]
  rownames(results) <- NULL
  results
}

test_that("MaAsLin3 contrast-wise evaluation preserves native complete results", {
  skip_if_not_installed("maaslin3")
  fit <- make_maaslin_contrast_fits()
  weights <- maaslin_contrast_matrix()
  for (median_comparison in c(FALSE, TRUE)) {
    for (subtract_median in c(FALSE, TRUE)) {
      withr::local_seed(130)
      native <- maaslin3::maaslin_contrast_test(
        fit, weights, correction = "BH", evaluate_only = "abundance",
        median_comparison_abundance = median_comparison,
        subtract_median = subtract_median
      )$fit_data_abundance$results
      native_seed <- .Random.seed
      set.seed(130)
      actual <- dar:::maaslin_model_contrasts(
        fit, weights, 0.05, "BH", median_comparison,
        engine_args = list(contrast = list(subtract_median = subtract_median))
      )
      expect_equal(sort_maaslin_contrasts(actual),
                   sort_maaslin_contrasts(native), tolerance = 1e-12)
      expect_identical(.Random.seed, native_seed)
    }
  }
})

test_that("MaAsLin3 retains estimable contrasts when another group is absent", {
  skip_if_not_installed("maaslin3")
  fit <- make_maaslin_contrast_fits(partial = TRUE)
  weights <- maaslin_contrast_matrix()
  for (median_comparison in c(FALSE, TRUE)) {
    withr::local_seed(130)
    expected <- lapply(seq_len(nrow(weights)), function(index) {
      maaslin3::maaslin_contrast_test(
        fit, weights[index, , drop = FALSE], correction = "BH",
        evaluate_only = "abundance",
        median_comparison_abundance = median_comparison
      )$fit_data_abundance$results
    }) |>
      dplyr::bind_rows()
    pvalues <- expected$pval_individual
    pvalues[!is.na(expected$error)] <- NA_real_
    expected$qval_individual <- stats::p.adjust(pvalues, method = "BH")
    expected$qval_joint <- expected$qval_individual
    set.seed(130)
    actual <- dar:::maaslin_model_contrasts(
      fit, weights, 0.05, "BH", median_comparison
    )
    expect_equal(sort_maaslin_contrasts(actual),
                 sort_maaslin_contrasts(expected), tolerance = 1e-12)
    selected <- actual[actual$feature == "taxon_1", ]
    expect_true(is.finite(selected$coef[selected$test == "B_A"]))
    expect_true(all(is.na(selected$coef[selected$test != "B_A"])))
    expect_true(all(is.na(selected$qval_individual[selected$test != "B_A"])))
    expect_identical(anyDuplicated(actual[c("feature", "test")]), 0L)
    expect_equal(nrow(actual), 24L)
  }
})

test_that("MaAsLin3 joint correction masks native errors without changing estimates", {
  skip_if_not_installed("maaslin3")
  testthat::local_mocked_bindings(
    exec_engine_stage = function(engine, stage, engine_args, fixed) {
      data <- data.frame(
        feature = c("one", "two"), test = rownames(fixed$contrast_mat),
        coef = c(1, 2), pval_individual = c(0.01, 1e-8),
        qval_individual = c(0.02, 2e-8), qval_joint = c(0.02, 2e-8),
        error = c(NA_character_, "Native error")
      )
      list(fit_data_abundance = list(results = data))
    },
    .package = "dar"
  )
  actual <- dar:::maaslin_model_contrasts(
    NULL, maaslin_contrast_matrix(), 0.05, "bonferroni", FALSE
  )
  expect_equal(actual$coef, rep(c(1, 2), 3))
  expect_equal(actual$pval_individual, rep(c(0.01, 1e-8), 3))
  expect_equal(actual$qval_individual, rep(c(0.03, NA_real_), 3))
  expect_identical(actual$qval_joint, actual$qval_individual)
})

test_that("MaAsLin3 contrast-wise evaluation preserves native mixed-model tests", {
  skip_if_not_installed("maaslin3")
  skip_if_not_installed("lmerTest")
  fit <- make_maaslin_contrast_fits(mixed = TRUE)
  weights <- maaslin_contrast_matrix(mixed = TRUE)
  for (median_comparison in c(FALSE, TRUE)) {
    withr::local_seed(130)
    native <- maaslin3::maaslin_contrast_test(
      fit, weights, correction = "BH", evaluate_only = "abundance",
      median_comparison_abundance = median_comparison
    )$fit_data_abundance$results
    native_seed <- .Random.seed
    set.seed(130)
    actual <- dar:::maaslin_model_contrasts(
      fit, weights, 0.05, "BH", median_comparison
    )
    expect_equal(sort_maaslin_contrasts(actual),
                 sort_maaslin_contrasts(native), tolerance = 1e-12)
    expect_identical(.Random.seed, native_seed)
  }
})

test_that("MaAsLin3 recovers only estimable mixed-model contrasts", {
  skip_if_not_installed("maaslin3")
  skip_if_not_installed("lmerTest")
  fit <- make_maaslin_contrast_fits(partial = TRUE, mixed = TRUE)
  weights <- maaslin_contrast_matrix(mixed = TRUE)
  withr::local_seed(130)
  actual <- dar:::maaslin_model_contrasts(
    fit, weights, 0.05, "BH", FALSE
  )
  native <- maaslin3::maaslin_contrast_test(
    fit, weights["B_A", , drop = FALSE], correction = "BH",
    evaluate_only = "abundance", median_comparison_abundance = FALSE
  )$fit_data_abundance$results
  columns <- c("feature", "test", "coef", "stderr", "pval_individual")
  expect_equal(
    sort_maaslin_contrasts(actual[actual$test == "B_A", columns]),
    sort_maaslin_contrasts(native[, columns]), tolerance = 1e-12
  )
  selected <- actual[actual$feature == "taxon_1", ]
  expect_true(is.finite(selected$coef[selected$test == "B_A"]))
  expect_true(all(is.na(selected$coef[selected$test != "B_A"])))
})

test_that("MaAsLin3 fits once and routes each planned contrast independently", {
  skip_if_not_installed("maaslin3")
  phy <- make_multilevel_phy()
  counts <- as(phyloseq::otu_table(phy), "matrix")
  metadata <- data.frame(phyloseq::sample_data(phy))
  counts[1, metadata$condition == "treated"] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  calls <- character()
  native_stage <- dar:::engine_stage_function
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      native <- native_stage(engine, stage)
      function(...) {
        calls <<- c(calls, stage)
        native(...)
      }
    },
    .package = "dar"
  )
  withr::local_seed(130)
  result <- dar:::run_maaslin_model(
    rec, 0, 0, 0, "TSS", "LOG", 0.05, "BH", FALSE, FALSE, FALSE,
    engine_args = dar:::normalize_engine_args("maaslin", list(
      fit = list(cores = 1), contrast = list(subtract_median = FALSE)
    ))
  )
  expect_equal(sum(calls == "fit"), 1L)
  expect_equal(sum(calls == "contrast"), 3L)
  expect_complete_model_result(result, rec)
  table <- dar:::flatten_model_output(result)
  plan <- dar:::compile_model(rec, "maaslin")$contrasts
  estimable_id <- plan$contrast_id[
    plan$numerator == "rescue" & plan$denominator == "control"
  ]
  selected <- table[table$taxa_id == rownames(counts)[1], ]
  expect_length(estimable_id, 1L)
  expect_true(is.finite(selected$effect[selected$contrast_id == estimable_id]))
  expect_true(all(is.na(selected$effect[selected$contrast_id != estimable_id])))
})
