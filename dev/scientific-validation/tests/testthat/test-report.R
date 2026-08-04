testthat::test_that("per-engine artifacts aggregate without mixing rows", {
  fixture <- validation_metric_fixture()
  metrics <- score_validation_results(fixture$results, fixture$truth)
  summary <- summarize_validation_metrics(metrics)
  gates <- evaluate_validation_gates(summary)
  input <- tempfile("dar-validation-input-")
  first <- file.path(input, "deseq")
  second <- file.path(input, "wilcox")
  write_validation_artifacts(
    first, fixture$results, data.frame(), fixture$truth, data.frame(),
    metrics, summary, gates
  )
  wilcox_results <- fixture$results
  wilcox_results$engine <- "wilcox"
  write_validation_artifacts(
    second, wilcox_results, data.frame(), fixture$truth, data.frame(),
    metrics, summary, gates
  )
  output <- tempfile("dar-validation-aggregate-")
  aggregate <- aggregate_validation_artifacts(
    input, output, render = FALSE, strict = FALSE
  )

  testthat::expect_equal(nrow(aggregate$artifacts$results), 12L)
  testthat::expect_equal(nrow(aggregate$artifacts$truth), 6L)
  testthat::expect_setequal(aggregate$artifacts$results$engine, c("deseq", "wilcox"))
  testthat::expect_true(file.exists(file.path(output, "validation-results.rds")))
})

testthat::test_that("quick-profile baseline comparison is informational", {
  fixture <- validation_metric_fixture()
  summary <- summarize_validation_metrics(
    score_validation_results(fixture$results, fixture$truth)
  )
  baseline <- data.frame(
    summary,
    direction = "lower",
    tolerance = 0,
    source_run_id = "accepted",
    stringsAsFactors = FALSE
  )
  input <- tempfile("dar-validation-quick-")
  write_validation_artifacts(
    file.path(input, "deseq"), fixture$results, data.frame(), fixture$truth,
    data.frame(), data.frame(), summary, data.frame()
  )
  aggregate <- aggregate_validation_artifacts(
    input, tempfile("dar-validation-quick-report-"), render = FALSE,
    strict = TRUE, baseline = baseline, enforce_baseline = FALSE
  )

  testthat::expect_false(aggregate$baseline_failure)
  testthat::expect_true(all(
    aggregate$artifacts$baseline_comparison$status == "informational"
  ))
})

testthat::test_that("execution aggregation distinguishes every terminal state", {
  expected <- expand.grid(
    engine = c("deseq", "wilcox", "ancom", "linda", "aldex"),
    scenario = "cross_sectional_signal",
    stringsAsFactors = FALSE
  )
  runs <- data.frame(
    engine = expected$engine[-5],
    scenario = expected$scenario[-5],
    replicate = 1L,
    status = c("success", "expected_skip", "engine_error", "dependency_missing"),
    stringsAsFactors = FALSE
  )

  status <- summarize_validation_execution(runs, expected)
  testthat::expect_equal(status$status, c(
    "success", "skipped_incompatible", "failed_engine",
    "failed_dependency", "missing_artifact"
  ))
})

testthat::test_that("report aggregation records failures without masking rendering", {
  fixture <- validation_metric_fixture()
  input <- tempfile("dar-validation-failing-input-")
  runs <- data.frame(
    engine = "ancom", scenario = "cross_sectional_signal", replicate = 1L,
    status = "engine_error", reason = "runtime failure", warnings = "",
    duration_seconds = 1
  )
  gates <- data.frame(status = "fail")
  write_validation_artifacts(
    file.path(input, "ancom"), fixture$results, runs, fixture$truth,
    data.frame(), data.frame(), data.frame(), gates
  )

  relaxed <- aggregate_validation_artifacts(
    input, tempfile("dar-validation-relaxed-"), render = FALSE, strict = FALSE
  )
  strict <- aggregate_validation_artifacts(
    input, tempfile("dar-validation-strict-"), render = FALSE, strict = TRUE
  )

  testthat::expect_true(relaxed$engine_failure)
  testthat::expect_true(relaxed$gate_failure)
  testthat::expect_identical(relaxed$status, 0L)
  testthat::expect_identical(strict$status, 1L)
})
