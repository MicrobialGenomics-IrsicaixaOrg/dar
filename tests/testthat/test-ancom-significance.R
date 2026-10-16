ancom_significance_recipe <- function() {
  recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
}

ancom_significance_fixture <- function(coefficient = "conditiontreated") {
  result <- tibble::tibble(taxon = paste0("taxon_", seq_len(6)))
  result[[paste0("lfc_", coefficient)]] <- c(1, -1, 0.5, -0.5, 2, -2)
  result[[paste0("q_", coefficient)]] <- rep(0.01, 6)
  result[[paste0("diff_", coefficient)]] <- c(TRUE, TRUE, FALSE, NA, TRUE, FALSE)
  result[[paste0("passed_ss_", coefficient)]] <- c(
    FALSE, TRUE, TRUE, TRUE, NA, TRUE
  )
  result[[paste0("diff_robust_", coefficient)]] <- c(
    FALSE, TRUE, FALSE, NA, NA, FALSE
  )
  result
}

mock_ancom_significance_fit <- function(rec, robust = TRUE) {
  compiled <- dar:::compile_model(rec, "ancom")
  coefficient <- dar:::reparameterize_model_contrast(
    compiled, compiled$contrasts[1, , drop = FALSE]
  )$coefficient
  result <- ancom_significance_fixture(coefficient)
  if (!robust) {
    result[[paste0("diff_robust_", coefficient)]] <- NULL
  }
  observed <- logical()
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      stopifnot(identical(engine, "ancom"), identical(stage, "fit"))
      function(pseudo_sens = FALSE, ...) {
        observed <<- c(observed, pseudo_sens)
        list(res = result)
      }
    },
    .package = "dar",
    .env = parent.frame()
  )
  function() observed
}

test_that("centralized ANCOM-BC2 uses native robust decisions by default", {
  skip_if_not_installed("ANCOMBC")
  rec <- ancom_significance_recipe()
  observed <- mock_ancom_significance_fit(rec)
  prepared <- prep(
    step_ancom(rec, id = "ancom__robust"), parallel = FALSE
  )
  result <- tidy_results(prepared)

  expect_identical(observed(), TRUE)
  expect_identical(
    result$significant, c(FALSE, TRUE, FALSE, NA, NA, FALSE)
  )
  expect_equal(result$adj_p_value, rep(0.01, 6))
  expect_equal(result$effect_size, c(1, -1, 0.5, -0.5, 2, -2))
  expect_identical(
    tidy_results(prepared, significant_only = TRUE)$taxa_id, "taxon_2"
  )
  expect_identical(
    cool(bake(prepared, count_cutoff = 1))$taxa_id, "taxon_2"
  )
})

test_that("explicit sensitivity analysis also selects native robust decisions", {
  skip_if_not_installed("ANCOMBC")
  rec <- ancom_significance_recipe()
  observed <- mock_ancom_significance_fit(rec)
  prepared <- prep(
    step_ancom(
      rec, engine_args = list(fit = list(pseudo_sens = TRUE)),
      id = "ancom__explicit_robust"
    ),
    parallel = FALSE
  )

  expect_identical(observed(), TRUE)
  expect_identical(
    tidy_results(prepared)$significant,
    c(FALSE, TRUE, FALSE, NA, NA, FALSE)
  )
})

test_that("disabling sensitivity analysis selects native basic decisions", {
  skip_if_not_installed("ANCOMBC")
  rec <- ancom_significance_recipe()
  observed <- mock_ancom_significance_fit(rec)
  prepared <- prep(
    step_ancom(
      rec, engine_args = list(fit = list(pseudo_sens = FALSE)),
      id = "ancom__basic"
    ),
    parallel = FALSE
  )

  expect_identical(observed(), FALSE)
  expect_identical(
    tidy_results(prepared)$significant,
    c(TRUE, TRUE, FALSE, NA, TRUE, FALSE)
  )
  expect_identical(
    tidy_results(prepared, significant_only = TRUE)$taxa_id,
    c("taxon_1", "taxon_2", "taxon_5")
  )
})

test_that("basic ANCOM-BC2 does not require robust columns", {
  skip_if_not_installed("ANCOMBC")
  rec <- ancom_significance_recipe()
  observed <- mock_ancom_significance_fit(rec, robust = FALSE)
  prepared <- prep(
    step_ancom(
      rec, engine_args = list(fit = list(pseudo_sens = FALSE)),
      id = "ancom__basic_only"
    ),
    parallel = FALSE
  )

  expect_identical(observed(), FALSE)
  expect_identical(
    tidy_results(prepared)$significant,
    c(TRUE, TRUE, FALSE, NA, TRUE, FALSE)
  )
})

test_that("robust ANCOM-BC2 does not silently fall back to basic calls or q", {
  skip_if_not_installed("ANCOMBC")
  rec <- ancom_significance_recipe()
  mock_ancom_significance_fit(rec, robust = FALSE)

  expect_error(
    prep(step_ancom(rec, id = "ancom__missing_robust"), parallel = FALSE),
    class = "dar_error_invalid_result_contract"
  )
})



test_that("ANCOM-BC2 significance columns must match the requested mode", {
  columns <- c(
    "diff_(Intercept)", "diff_conditiontreated",
    "diff_robust_(Intercept)", "diff_robust_conditiontreated"
  )
  expect_identical(
    dar:::ancom_significance_column(columns, coefficient = "conditiontreated"),
    "diff_robust_conditiontreated"
  )
  expect_identical(
    dar:::ancom_significance_column(
      columns, pseudo_sens = FALSE, coefficient = "conditiontreated"
    ),
    "diff_conditiontreated"
  )
  expect_error(
    dar:::ancom_significance_column(
      "diff_robust_conditiontreated", pseudo_sens = FALSE,
      coefficient = "conditiontreated"
    ),
    class = "dar_error_invalid_result_contract"
  )
})
