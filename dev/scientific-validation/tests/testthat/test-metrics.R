testthat::test_that("scientific metrics are computed without mixing hypotheses", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth)

  testthat::expect_equal(metrics$type1_error, 0.25)
  testthat::expect_equal(metrics$fdr, 1 / 3)
  testthat::expect_equal(metrics$power, 1)
  testthat::expect_equal(metrics$direction_recovery, 1)
  testthat::expect_equal(metrics$finite_rate, 1)
  testthat::expect_true(is.finite(metrics$abs_relative_bias))

  second_engine <- fixture$results
  second_engine$engine <- "ancom"
  second_engine$effect_metric <- "bias_corrected_log_fold_change"
  second_engine$truth_multiplier <- 1
  combined <- score_validation_results(
    rbind(fixture$results, second_engine), fixture$truth
  )
  testthat::expect_equal(nrow(combined), 2L)
  testthat::expect_setequal(combined$engine, c("deseq", "ancom"))
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
