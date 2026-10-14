testthat::test_that("scientific metrics are computed without mixing hypotheses", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)

  testthat::expect_equal(metrics$type1_error, 0.25)
  testthat::expect_equal(metrics$fdr, 1 / 3)
  testthat::expect_equal(metrics$power, 1)
  testthat::expect_equal(metrics$direction_recovery, 1)
  testthat::expect_equal(metrics$finite_rate, 1)
  testthat::expect_true(is.finite(metrics$relative_mae))
  testthat::expect_equal(metrics$n_expected, 6L)
  testthat::expect_equal(metrics$n_returned, 6L)
  testthat::expect_equal(metrics$n_missing, 0L)
  testthat::expect_equal(metrics$n_signals, 2L)
  testthat::expect_equal(metrics$n_tested, 6L)
  testthat::expect_equal(metrics$n_signal_finite, 2L)
  testthat::expect_equal(metrics$n_intervals, 2L)

  second_engine <- fixture$results
  second_engine$engine <- "ancom"
  second_engine$effect_metric <- "bias_corrected_log_fold_change"
  second_engine$truth_multiplier <- 1
  second_run <- fixture$runs
  second_run$engine <- "ancom"
  combined <- score_validation_results(
    rbind(fixture$results, second_engine), fixture$truth,
    rbind(fixture$runs, second_run)
  )
  testthat::expect_equal(nrow(combined), 2L)
  testthat::expect_setequal(combined$engine, c("deseq", "ancom"))
})

testthat::test_that("a reversed contrast direction fails its gate", {
  fixture <- validation_metric_fixture(reverse = TRUE)
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
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
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
  gates <- evaluate_validation_gates(summarize_validation_metrics(metrics))

  testthat::expect_false(any(gates$status == "fail"))
  testthat::expect_true(any(gates$status == "informational"))
})

testthat::test_that("validation artifacts are machine readable", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
  summary <- summarize_validation_metrics(metrics)
  gates <- evaluate_validation_gates(summary)
  output <- tempfile("dar-validation-")
  write_validation_artifacts(
    output, fixture$results, fixture$runs, fixture$truth, data.frame(),
    metrics, summary, gates
  )

  expected <- c(
    "results.csv", "runs.csv", "truth.csv", "manifests.csv", "metrics.csv",
    "summary.csv", "gates.csv", "versions.csv", "validation-results.rds",
    "session-info.txt"
  )
  testthat::expect_setequal(list.files(output), expected)
  testthat::expect_equal(nrow(utils::read.csv(file.path(output, "results.csv"))), 6L)
  restored <- score_validation_results(
    utils::read.csv(file.path(output, "results.csv")),
    utils::read.csv(file.path(output, "truth.csv")),
    utils::read.csv(file.path(output, "runs.csv"))
  )
  testthat::expect_equal(restored, metrics, tolerance = 1e-12)
})

testthat::test_that("longitudinal power gates target only strong contrasts", {
  summary <- data.frame(
    engine = "linda",
    scenario = "longitudinal_interaction",
    contrast_id = c(
      "condition[treated-control]@time[1]",
      "condition[treated-control]:time[1-0]"
    ),
    metric = "power",
    estimate = c(0.9, 0.1),
    lower = c(0.8, 0.05),
    upper = c(1, 0.2),
    n_replicates = 30L,
    stringsAsFactors = FALSE
  )
  gates <- evaluate_validation_gates(summary)

  testthat::expect_equal(gates$status, c("pass", "not_applicable"))
})

testthat::test_that("LinDA comparable metrics have active gates", {
  summary <- data.frame(
    engine = "linda",
    scenario = "cross_sectional_signal",
    contrast_id = "condition[treated-control]",
    metric = c("coverage", "relative_mae"),
    estimate = c(0.93, 0.08),
    lower = c(0.90, 0.05),
    upper = c(0.96, 0.11),
    n_replicates = 50L,
    stringsAsFactors = FALSE
  )

  testthat::expect_equal(evaluate_validation_gates(summary)$status, c("pass", "pass"))
})

validation_four_signal_fixture <- function() {
  fixture <- validation_metric_fixture()
  fixture$truth <- fixture$truth[1:4, ]
  fixture$truth$true_effect <- c(1, 1, -1, -1)
  fixture$truth$is_null <- FALSE
  fixture$truth$expected_direction <- sign(fixture$truth$true_effect)
  fixture$results <- fixture$results[1:4, ]
  fixture$results$effect <- fixture$truth$true_effect / log(2)
  fixture$results$padj <- 0.01
  fixture
}

testthat::test_that("absent taxa remain in scientific denominators", {
  fixture <- validation_four_signal_fixture()
  partial <- fixture$results[1, , drop = FALSE]
  universe <- validation_result_universe(partial, fixture$truth, fixture$runs)
  metrics <- score_validation_results(partial, fixture$truth, fixture$runs)

  testthat::expect_equal(metrics$power, 0.25)
  testthat::expect_equal(metrics$direction_recovery, 0.25)
  testthat::expect_equal(metrics$finite_rate, 0.25)
  testthat::expect_equal(metrics$n_expected, 4L)
  testthat::expect_equal(metrics$n_returned, 1L)
  testthat::expect_equal(metrics$n_missing, 3L)
  testthat::expect_equal(metrics$n_signal_finite, 1L)
  testthat::expect_equal(metrics$n_intervals, 1L)
  testthat::expect_true(all(is.na(universe$effect[!universe$returned])))
  testthat::expect_true(all(is.na(universe$padj[!universe$returned])))
  testthat::expect_equal(nrow(partial), 1L)
})

testthat::test_that("absent contrasts and replicates are not silently dropped", {
  fixture <- validation_four_signal_fixture()
  second_truth <- fixture$truth
  second_truth$contrast_id <- "condition[treated-control]@time[1]"
  metrics <- score_validation_results(
    fixture$results, rbind(fixture$truth, second_truth), fixture$runs
  )
  absent <- metrics[metrics$contrast_id == second_truth$contrast_id[1], ]

  testthat::expect_equal(nrow(metrics), 2L)
  testthat::expect_equal(absent$n_expected, 4L)
  testthat::expect_equal(absent$n_missing, 4L)
  testthat::expect_equal(absent$power, 0)
  testthat::expect_equal(absent$finite_rate, 0)
  testthat::expect_equal(absent$direction_recovery, 0)
  testthat::expect_true(is.na(absent$relative_mae))
  testthat::expect_true(is.na(absent$coverage))

  second_truth <- fixture$truth
  second_truth$replicate <- 2L
  second_run <- fixture$runs
  second_run$replicate <- 2L
  metrics <- score_validation_results(
    fixture$results, rbind(fixture$truth, second_truth),
    rbind(fixture$runs, second_run)
  )
  testthat::expect_equal(metrics$power, c(1, 0))
  summary <- summarize_validation_metrics(metrics)
  power <- summary[summary$metric == "power", ]
  testthat::expect_equal(power$n_replicates, 2L)
  testthat::expect_equal(power$estimate, 0.5)
})

testthat::test_that("failed fits count as unavailable but capability skips do not", {
  fixture <- validation_four_signal_fixture()
  for (status in c("success", "engine_error")) {
    fixture$runs$status <- status
    metrics <- score_validation_results(data.frame(), fixture$truth, fixture$runs)
    testthat::expect_equal(nrow(metrics), 1L)
    testthat::expect_equal(metrics$power, 0)
    testthat::expect_equal(metrics$finite_rate, 0)
    testthat::expect_equal(metrics$n_missing, 4L)
  }
  for (status in c("expected_skip", "dependency_missing")) {
    fixture$runs$status <- status
    testthat::expect_equal(
      nrow(score_validation_results(data.frame(), fixture$truth, fixture$runs)),
      0L
    )
  }
})

testthat::test_that("non-finite estimates cannot improve recovery by exclusion", {
  fixture <- validation_four_signal_fixture()
  fixture$results$effect[2:4] <- c(NA, Inf, NaN)
  fixture$results$padj[2:4] <- c(NA, Inf, NaN)
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)

  testthat::expect_equal(metrics$n_returned, 4L)
  testthat::expect_equal(metrics$n_missing, 0L)
  testthat::expect_equal(metrics$n_tested, 1L)
  testthat::expect_equal(metrics$n_finite, 1L)
  testthat::expect_equal(metrics$power, 0.25)
  testthat::expect_equal(metrics$direction_recovery, 0.25)
  testthat::expect_equal(metrics$finite_rate, 0.25)
})

testthat::test_that("balanced attenuation does not cancel absolute estimation error", {
  fixture <- validation_four_signal_fixture()
  fixture$results$effect <- 0
  fixture$results$padj <- 0.9
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)

  testthat::expect_equal(metrics$bias, 0)
  testthat::expect_equal(metrics$relative_mae, 1)
  testthat::expect_equal(metrics$power, 0)
  testthat::expect_equal(metrics$finite_rate, 1)
  repeated <- metrics[rep(1L, 30L), ]
  repeated$replicate <- seq_len(30L)
  gates <- evaluate_validation_gates(summarize_validation_metrics(repeated))
  testthat::expect_equal(gates$status[gates$metric == "relative_mae"], "fail")

  fixture$results$effect <- -fixture$truth$true_effect / log(2)
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
  summary <- summarize_validation_metrics(metrics)
  error <- summary[summary$metric == "relative_mae", ]
  testthat::expect_equal(error$estimate, 2)
  testthat::expect_equal(error$upper, 2)
})

testthat::test_that("null and incomparable effects have no relative error or coverage", {
  fixture <- validation_metric_fixture()
  fixture$truth$true_effect <- 0
  fixture$truth$is_null <- TRUE
  fixture$truth$expected_direction <- 0
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
  testthat::expect_true(is.na(metrics$power))
  testthat::expect_true(is.na(metrics$direction_recovery))
  testthat::expect_true(is.na(metrics$relative_mae))
  testthat::expect_true(is.na(metrics$coverage))

  fixture <- validation_four_signal_fixture()
  config <- validation_engine("wilcox")
  fixture$runs$engine <- fixture$results$engine <- "wilcox"
  fixture$results$effect_metric <- config$effect_metric
  fixture$results$truth_multiplier <- config$truth_multiplier
  fixture$results$comparable <- config$comparable
  metrics <- score_validation_results(fixture$results, fixture$truth, fixture$runs)
  testthat::expect_equal(metrics$power, 1)
  testthat::expect_equal(metrics$n_intervals, 0L)
  testthat::expect_true(is.na(metrics$bias))
  testthat::expect_true(is.na(metrics$relative_mae))
  testthat::expect_true(is.na(metrics$coverage))
})

testthat::test_that("hypothesis coverage is independent for each engine", {
  fixture <- validation_four_signal_fixture()
  second_results <- fixture$results
  second_results$engine <- "ancom"
  second_results$effect_metric <- validation_engine("ancom")$effect_metric
  second_results$truth_multiplier <- 1
  second_results$effect <- fixture$truth$true_effect
  second_run <- fixture$runs
  second_run$engine <- "ancom"
  metrics <- score_validation_results(
    rbind(fixture$results[1, ], second_results), fixture$truth,
    rbind(fixture$runs, second_run)
  )

  testthat::expect_equal(metrics$power[metrics$engine == "deseq"], 0.25)
  testthat::expect_equal(metrics$power[metrics$engine == "ancom"], 1)
  testthat::expect_equal(metrics$n_missing[metrics$engine == "deseq"], 3L)
  testthat::expect_equal(metrics$n_missing[metrics$engine == "ancom"], 0L)
})

testthat::test_that("scientific keys and effect metadata are validated before scoring", {
  fixture <- validation_four_signal_fixture()
  testthat::expect_error(
    score_validation_results(
      rbind(fixture$results, fixture$results[1, ]), fixture$truth, fixture$runs
    ), "duplicated"
  )
  invalid <- fixture$results
  invalid$taxa_id[1] <- "ghost_taxon"
  testthat::expect_error(
    score_validation_results(invalid, fixture$truth, fixture$runs), "expected universe"
  )
  invalid <- fixture$results
  invalid$effect_metric <- "unknown_scale"
  testthat::expect_error(
    score_validation_results(invalid, fixture$truth, fixture$runs), "effect metadata"
  )
  invalid <- fixture$results
  invalid$truth_multiplier <- 3
  testthat::expect_error(
    score_validation_results(invalid, fixture$truth, fixture$runs), "effect metadata"
  )
  testthat::expect_error(
    score_validation_results(
      fixture$results, rbind(fixture$truth, fixture$truth[1, ]), fixture$runs
    ), "truth.*duplicated"
  )
  testthat::expect_error(
    score_validation_results(fixture$results, fixture$truth, data.frame()), "no recorded runs"
  )
  testthat::expect_error(
    score_validation_results(fixture$results, fixture$truth, rbind(fixture$runs, fixture$runs)),
    "runs.*duplicated"
  )
  invalid <- fixture$runs
  invalid$scenario <- "missing_scenario"
  testthat::expect_error(
    score_validation_results(data.frame(), fixture$truth, invalid), "absent from the truth"
  )
  invalid <- fixture$runs
  invalid$status <- "unknown_status"
  testthat::expect_error(
    score_validation_results(data.frame(), fixture$truth, invalid), "unknown execution statuses"
  )
})

testthat::test_that("relative error gates exclude weak temporal contrasts", {
  summary <- data.frame(
    engine = "linda", scenario = "longitudinal_interaction",
    contrast_id = c(
      "condition[treated-control]@time[1]", "time[1-0]@condition[control]"
    ),
    metric = "relative_mae", estimate = 1, lower = 1, upper = 1,
    n_replicates = 30L
  )
  testthat::expect_equal(
    evaluate_validation_gates(summary)$status, c("fail", "not_applicable")
  )
})
