benchmark_fixture_truth <- function() {
  data.frame(
    taxa_id = paste0("t", 1:4), contrast_id = "condition_treatment_vs_control",
    beta_latent = c(1, -1, 0, 0), beta_clr_U = c(1, -1, 0, 0),
    beta_median_U = c(1, -1, 0, 0), retained = rep(TRUE, 4)
  )
}

benchmark_fixture_result <- function() {
  data.frame(
    taxa_id = paste0("t", 1:4), contrast_id = "condition_treatment_vs_control",
    effect_size = c(1, 1, 0.1, NA_real_),
    significant = c(TRUE, TRUE, TRUE, FALSE), padj = c(0.9, 0.9, 0.9, NA_real_)
  )
}

testthat::test_that("selection scoring honors native calls and direction", {
  score <- benchmark_score(benchmark_fixture_result(), benchmark_fixture_truth(),
                           "linda_pilot", "P2", 1)
  testthat::expect_equal(nrow(score), 3L)
  testthat::expect_equal(score$R, rep(3L, 3))
  testthat::expect_equal(score$fdp, rep(1 / 3, 3))
  testthat::expect_equal(score$power_endtoend, rep(1, 3))
  testthat::expect_equal(score$directional_power, rep(0.5, 3))
  testthat::expect_equal(score$error_directional, rep(0.5, 3))
  testthat::expect_equal(score$directional_fdp, rep(2 / 3, 3))
  testthat::expect_equal(score$n_q_finite, rep(3L, 3))
})

testthat::test_that("availability reads the canonical tidy-results adjusted p column", {
  result <- benchmark_fixture_result()
  result$adj_p_value <- c(0.01, NA_real_, 0.05, Inf)
  result$padj <- rep(0.02, 4)
  score <- benchmark_score(result, benchmark_fixture_truth(), "linda_pilot", "P2", 1)
  testthat::expect_equal(score$n_q_finite, rep(2L, 3))
  testthat::expect_equal(score$R, rep(3L, 3))
})

testthat::test_that("removed signals and missing results are not recovered", {
  truth <- benchmark_fixture_truth()
  truth$retained[[2L]] <- FALSE
  truth$beta_clr_original <- truth$beta_clr_U
  truth$beta_median_original <- truth$beta_median_U
  truth$beta_clr_U <- c(2 / 3, NA_real_, -1 / 3, -1 / 3)
  truth$beta_median_U <- c(1, NA_real_, 0, 0)
  result <- benchmark_fixture_result()[c(1L, 3L), ]
  score <- benchmark_score(result, truth, "linda_pilot", "P2", 1)
  testthat::expect_equal(score$power_endtoend[[1L]], 0.5)
  testthat::expect_equal(score$power_U[[1L]], 1)
  testthat::expect_equal(score$signal_retention[[1L]], 0.5)
  testthat::expect_equal(score$n_missing, rep(1L, 3))
  testthat::expect_equal(score$fdp[[2L]], 0.5)
  testthat::expect_equal(score$fdp_U[[2L]], 0)
  testthat::expect_equal(score$estimand_shift[[2L]], 1 / 3)
  testthat::expect_error(
    benchmark_score(result, truth[c(setdiff(names(truth), "beta_clr_original"))],
                    "linda_pilot", "P2", 1), "Original truth"
  )
})

testthat::test_that("native missing calls are distinct from valid negatives", {
  result <- benchmark_fixture_result()
  result$significant[[4L]] <- NA
  score <- benchmark_score(result, benchmark_fixture_truth(), "linda_pilot", "P2", 1)
  testthat::expect_true(all(score$available))
  testthat::expect_equal(score$n_significance_missing, rep(1L, 3))
  testthat::expect_equal(score$n_effect_nonfinite, rep(1L, 3))
  result$significant <- rep(FALSE, 4)
  score <- benchmark_score(result, benchmark_fixture_truth(), "linda_pilot", "P2", 1)
  testthat::expect_equal(score$fdp, rep(0, 3))
  testthat::expect_true(all(is.na(score$error_directional)))
})

testthat::test_that("contract-invalid calls and failures cannot become FDP zero", {
  for (value in c(NA_real_, Inf, 0)) {
    result <- benchmark_fixture_result()
    result$effect_size[[1L]] <- value
    score <- benchmark_score(result, benchmark_fixture_truth(), "ancom_pilot", "P2", 1)
    testthat::expect_equal(score$status, rep("contract_invalid", 3))
    testthat::expect_true(all(is.na(score$fdp)))
    testthat::expect_equal(score$conservative_power, rep(0, 3))
    testthat::expect_equal(score$n_returned, rep(4L, 3))
  }
  for (result in list(
    rbind(benchmark_fixture_result(), benchmark_fixture_result()[1L, ]),
    transform(benchmark_fixture_result(), taxa_id = paste0("unknown", 1:4))
  )) {
    score <- benchmark_score(result, benchmark_fixture_truth(), "ancom_pilot", "P2", 1)
    testthat::expect_equal(score$status, rep("contract_invalid", 3))
  }
  score <- benchmark_score(NULL, benchmark_fixture_truth(), "majority_3of4", "P2", 1,
                           status = "unavailable")
  testthat::expect_true(all(is.na(score$fdp)))
  testthat::expect_true(all(is.na(score$power_endtoend)))
  testthat::expect_equal(score$conservative_directional_power, rep(0, 3))
})

testthat::test_that("null FDR is the probability of any discovery", {
  truth <- benchmark_fixture_truth()
  truth$beta_latent <- truth$beta_clr_U <- truth$beta_median_U <- rep(0, 4)
  metrics <- do.call(rbind, lapply(1:20, function(replicate) {
    result <- benchmark_fixture_result()
    result$significant <- rep(replicate <= 4L, 4)
    result$effect_size <- rep(1, 4)
    benchmark_score(result, truth, "linda_pilot", "P1", replicate)
  }))
  summary <- benchmark_summarize(metrics, bootstrap = 10)
  testthat::expect_equal(summary$fdr, rep(0.2, 3))
  testthat::expect_equal(summary$null_events, rep(4L, 3))
  testthat::expect_true(all(is.na(summary$directional_power)))
  ci <- benchmark_null_interval(rep(FALSE, 20))
  testthat::expect_equal(unname(ci[["upper_one_sided"]]), 1 - 0.05^(1 / 20))
})

testthat::test_that("summaries use paired datasets, preserve failures and RNG", {
  truth <- benchmark_fixture_truth()
  result <- benchmark_fixture_result()
  metrics <- do.call(rbind, lapply(1:3, function(replicate) {
    reference <- benchmark_score(result, truth, "linda_pilot", "P2", replicate)
    status <- if (replicate == 3L) "unavailable" else "success"
    consensus <- benchmark_score(result[1L, ], truth, "majority_3of4", "P2", replicate, status)
    rbind(reference, consensus)
  }))
  set.seed(123)
  before <- .Random.seed
  summary <- benchmark_summarize(metrics, bootstrap = 20, planned_replicates = 3)
  testthat::expect_identical(.Random.seed, before)
  row <- summary[summary$procedure == "majority_3of4" & summary$truth == "truth_latent", ]
  testthat::expect_equal(row$n_valid, 2L)
  testthat::expect_equal(row$n_paired_fdp, 2L)
  testthat::expect_equal(row$n_paired_power, 2L)
  testthat::expect_equal(row$delta_fdp, -1 / 3)
  testthat::expect_equal(row$delta_directional_power, 0)
  testthat::expect_equal(row$fdr_identification_upper, 1 / 3)
  testthat::expect_equal(row$conservative_power, 1 / 3)
  testthat::expect_identical(summary, benchmark_summarize(metrics, bootstrap = 20,
                                                        planned_replicates = 3))
  testthat::expect_error(benchmark_summarize(rbind(metrics, metrics[1L, ])), "Invalid")
})
