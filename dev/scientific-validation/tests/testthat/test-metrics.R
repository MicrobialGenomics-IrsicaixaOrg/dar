validation_metric_fixture <- function(reverse = FALSE) {
  truth <- expand.grid(
    taxa_id = sprintf("taxon_%02d", 1:6),
    contrast_id = "condition[treated-control]",
    stringsAsFactors = FALSE
  )
  truth$scenario <- "cross_sectional_signal"
  truth$replicate <- 1L
  truth$true_effect <- c(1, -1, 0, 0, 0, 0)
  truth$is_null <- truth$true_effect == 0
  truth$expected_direction <- sign(truth$true_effect)
  truth$truth_scale <- "centered_log_abundance"
  effects <- c(1.4, -1.4, 0.1, 0.1, 0.1, 0.1)
  if (reverse) effects[1:2] <- -effects[1:2]
  results <- data.frame(
    engine = "deseq",
    scenario = truth$scenario,
    replicate = truth$replicate,
    taxa_id = truth$taxa_id,
    contrast_id = truth$contrast_id,
    effect = effects,
    padj = c(0.01, 0.02, 0.01, 0.5, 0.5, 0.5),
    std_error = 0.2,
    effect_metric = "log2_fold_change",
    truth_multiplier = 1 / log(2),
    comparable = TRUE
  )
  list(results = results, truth = truth)
}

testthat::test_that("scientific metrics are computed without mixing hypotheses", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth)

  testthat::expect_equal(metrics$type1_error, 0.25)
  testthat::expect_equal(metrics$fdr, 1 / 3)
  testthat::expect_equal(metrics$power, 1)
  testthat::expect_equal(metrics$direction_recovery, 1)
  testthat::expect_equal(metrics$finite_rate, 1)
  testthat::expect_true(is.finite(metrics$abs_relative_bias))
})

testthat::test_that("a reversed contrast direction fails its gate", {
  fixture <- validation_metric_fixture(reverse = TRUE)
  metrics <- score_validation_results(fixture$results, fixture$truth)
  repeated <- metrics[rep(1L, 20L), , drop = FALSE]
  repeated$replicate <- seq_len(20L)
  summary <- summarize_validation_metrics(repeated)
  gates <- evaluate_validation_gates(summary)
  direction <- gates[gates$metric == "direction_recovery", , drop = FALSE]

  testthat::expect_equal(metrics$direction_recovery, 0)
  testthat::expect_equal(direction$status, "fail")
  testthat::expect_true(direction_reversal_self_test())
})

testthat::test_that("insufficient Monte Carlo replication is informational", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth)
  gates <- evaluate_validation_gates(summarize_validation_metrics(metrics))

  testthat::expect_false(any(gates$status == "fail"))
  testthat::expect_true(any(gates$status == "informational"))
})

testthat::test_that("validation artifacts are machine readable", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth)
  summary <- summarize_validation_metrics(metrics)
  gates <- evaluate_validation_gates(summary)
  output <- tempfile("dar-validation-")
  write_validation_artifacts(
    output, fixture$results, data.frame(), fixture$truth, data.frame(),
    metrics, summary, gates
  )

  expected <- c(
    "results.csv", "runs.csv", "truth.csv", "manifests.csv", "metrics.csv",
    "summary.csv", "gates.csv", "validation-results.rds", "session-info.txt"
  )
  testthat::expect_setequal(list.files(output), expected)
  testthat::expect_equal(nrow(utils::read.csv(file.path(output, "results.csv"))), 6L)
})
