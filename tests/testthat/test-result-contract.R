contract_recipe <- function() {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  engines <- c("deseq", "aldex", "ancom", "corncob", "maaslin", "wilcox", "lefse")
  for (engine in engines) {
    arguments <- list(subclass = engine, id = paste0(engine, "__contract"))
    if (identical(engine, "wilcox")) {
      arguments$norm_method <- "clr"
    }
    if (identical(engine, "corncob")) {
      arguments$test <- "Wald"
      arguments$boot <- FALSE
    }
    rec <- dar:::add_step(rec, do.call(dar:::step, arguments))
  }
  rec
}

contract_results <- function(rec = contract_recipe()) {
  contrast <- dar:::resolve_model(rec)$contrast_plan[1, , drop = FALSE]
  common <- tibble::tibble(
    taxa_id = "taxon_1",
    taxa = "Species_1",
    contrast_id = contrast$contrast_id,
    comparison = contrast$comparison,
    contrast_type = contrast$contrast_type,
    var = contrast$var,
    effect = 1,
    signif = TRUE
  )
  list(
    deseq__contract = list(model = dplyr::mutate(
      common, pvalue = 0.01, padj = 0.02
    )),
    aldex__contract = list(model = dplyr::mutate(
      common, we.ep = 0.011, we.eBH = 0.021
    )),
    ancom__contract = list(model = dplyr::mutate(
      common, p = 0.012, q = 0.022
    )),
    corncob__contract = list(model = dplyr::mutate(
      common, pval = 0.013, padj = 0.023
    )),
    maaslin__contract = list(model = dplyr::mutate(
      common, pval = 0.014, qval = 0.024
    )),
    wilcox__contract = list(model = dplyr::mutate(
      common, p = 0.015, padj = 0.025
    )),
    lefse__contract = list(model = dplyr::mutate(
      common, pvalue = 0.016, adjp = 0.026
    ))
  )
}

contract_prepared <- function() {
  rec <- contract_recipe()
  results <- contract_results(rec)
  dar:::prep_recipe(
    rec,
    results,
    list(),
    execution = list(
      contrasts = dar:::resolve_model(rec)$contrast_plan,
      executed_steps = names(results),
      skipped_steps = tibble::tibble(
        step_id = character(), engine = character(), compatible = logical(),
        reason = character()
      )
    )
  )
}

test_that("tidy_results exposes the canonical contract for every engine", {
  prepared <- contract_prepared()
  raw <- prepared@results
  result <- tidy_results(prepared)

  expect_named(result, dar:::tidy_result_columns())
  expect_s3_class(result, "tbl_df")
  expect_true(all(vapply(
    result[c("effect_size", "p_value", "adj_p_value")],
    is.numeric,
    logical(1)
  )))
  expect_true(is.logical(result$significant))
  expect_true(all(vapply(
    result[setdiff(
      dar:::tidy_result_columns(),
      c("effect_size", "p_value", "adj_p_value", "significant")
    )],
    is.character,
    logical(1)
  )))
  expect_identical(prepared@results, raw)
  expect_equal(nrow(result), 7L)
  expect_equal(
    result$method,
    c("deseq", "aldex", "ancom", "corncob", "maaslin", "wilcox", "lefse")
  )
  expect_equal(result$step_id, names(raw))
  expect_equal(result$p_value, seq(0.01, 0.016, by = 0.001))
  expect_equal(result$adj_p_value, seq(0.02, 0.026, by = 0.001))
  expect_equal(result$effect_metric, c(
    "log2_fold_change",
    "standardized_clr_effect",
    "bias_corrected_log_fold_change",
    "logit_mean_abundance_coefficient",
    "transformed_abundance_coefficient",
    "median_clr_difference",
    "lda_score"
  ))
})

test_that("tidy_results supports filtering, step selection and empty recipes", {
  prepared <- contract_prepared()
  prepared@results$aldex__contract$model$signif <- FALSE

  expect_equal(nrow(tidy_results(prepared, significant_only = TRUE)), 6L)
  expect_equal(
    tidy_results(prepared, steps = "aldex__contract")$method,
    "aldex"
  )

  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  empty <- dar:::prep_recipe(rec, list(), list())
  expect_identical(tidy_results(empty), dar:::empty_tidy_results())
})

test_that("modeled ALDEx2 p-values are read without changing native columns", {
  rec <- contract_recipe()
  results <- contract_results(rec)
  compiled <- dar:::compile_model(rec, "aldex")
  contrast <- compiled$contrasts[1, , drop = FALSE]
  coefficient <- dar:::reparameterize_model_contrast(
    compiled, contrast
  )$coefficient
  native <- results$aldex__contract$model |>
    dplyr::select(-"we.ep", -"we.eBH")
  native[[paste0(coefficient, ":pval")]] <- 0.017
  native[[paste0(coefficient, ":pval.padj")]] <- 0.027
  native$padj <- 0.027
  results$aldex__contract$model <- native
  prepared <- dar:::prep_recipe(rec, results, list())
  before <- prepared@results$aldex__contract

  result <- tidy_results(prepared, steps = "aldex__contract")

  expect_equal(result$p_value, 0.017)
  expect_equal(result$adj_p_value, 0.027)
  expect_identical(prepared@results$aldex__contract, before)
})

test_that("tidy_results distinguishes multiple configurations of one engine", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(
    rec, dar:::step(
      "deseq", test = "Wald", fitType = "mean", id = "deseq__one"
    )
  )
  rec <- dar:::add_step(
    rec, dar:::step(
      "deseq", test = "Wald", fitType = "local", id = "deseq__two"
    )
  )
  contrast <- dar:::resolve_model(rec)$contrast_plan[1, , drop = FALSE]
  row <- tibble::tibble(
    taxa_id = "taxon_1", taxa = "Species_1",
    contrast_id = contrast$contrast_id, comparison = contrast$comparison,
    contrast_type = contrast$contrast_type, var = contrast$var,
    effect = 1, pvalue = 0.01, padj = 0.02, signif = TRUE
  )
  prepared <- dar:::prep_recipe(
    rec,
    list(deseq__one = list(model = row), deseq__two = list(model = row)),
    list()
  )
  result <- tidy_results(prepared)

  expect_equal(result$method, c("deseq", "deseq"))
  expect_equal(result$step_id, c("deseq__one", "deseq__two"))
  expect_silent(dar:::validate_tidy_result_contract(result))
})

test_that("tidy_results ordering follows requested steps, contrasts and taxa", {
  rec <- contract_recipe()
  results <- contract_results(rec)
  first <- results$deseq__contract$model
  second <- dplyr::mutate(
    first, taxa_id = "taxon_2", taxa = "Species_2", effect = -1
  )
  results$deseq__contract$model <- dplyr::bind_rows(second, first)
  prepared <- dar:::prep_recipe(rec, results, list())

  result <- tidy_results(
    prepared,
    steps = c("lefse__contract", "deseq__contract")
  )

  expect_identical(
    unique(result$step_id),
    c("lefse__contract", "deseq__contract")
  )
  deseq <- dplyr::filter(result, .data$step_id == "deseq__contract")
  expect_identical(deseq$taxa_id, c("taxon_1", "taxon_2"))
  expect_equal(deseq$effect_size, c(1, -1))
})

test_that("the result contract rejects invalid probabilities and duplicate keys", {
  prepared <- contract_prepared()
  invalid <- tidy_results(prepared)
  invalid$p_value[[1]] <- 2
  expect_error(
    dar:::validate_tidy_result_contract(invalid),
    class = "dar_error_invalid_result_contract"
  )

  duplicated <- dplyr::bind_rows(invalid[2, ], invalid[2, ])
  expect_error(
    dar:::validate_tidy_result_contract(duplicated),
    class = "dar_error_invalid_result_contract"
  )

  invalid_id <- tidy_results(prepared)
  invalid_id$step_id[[1]] <- ""
  expect_error(
    dar:::validate_tidy_result_contract(invalid_id),
    class = "dar_error_invalid_result_contract"
  )
})

test_that("unavailable and malformed engine results fail contextually", {
  prepared <- contract_prepared()
  prepared@execution$skipped_steps <- tibble::tibble(
    step_id = "wilcox__skipped", engine = "wilcox", compatible = FALSE,
    reason = "does not support adjusted confounders"
  )
  expect_error(
    tidy_results(prepared, steps = "wilcox__skipped"),
    "Skipped step",
    class = "dar_error_invalid_step_id"
  )

  prepared@results$deseq__contract$model$effect <- NULL
  expect_error(
    tidy_results(prepared, steps = "deseq__contract"),
    class = "dar_error_invalid_result_contract"
  )
})

test_that("legacy results receive deterministic contrast identifiers", {
  data(test_prep_rec)
  result <- tidy_results(test_prep_rec)

  expect_true(all(startsWith(result$contrast_id, paste0(result$var, "["))))
  expect_true(all(result$contrast_type == "legacy_pairwise"))
  expect_true(all(is.na(result$numerator)))
  expect_true(all(is.na(result$denominator)))
  expect_identical(
    anyDuplicated(result[c("taxa_id", "contrast_id", "step_id")]), 0L
  )
})
