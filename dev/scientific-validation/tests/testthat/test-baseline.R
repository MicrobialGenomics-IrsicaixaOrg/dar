testthat::test_that("baseline comparison tolerates Monte Carlo movement", {
  baseline <- data.frame(
    engine = "linda",
    scenario = "cross_sectional_signal",
    contrast_id = "condition[treated-control]",
    metric = "coverage",
    estimate = 0.93,
    lower = 0.90,
    upper = 0.96,
    direction = "lower",
    tolerance = 0.05,
    source_run_id = "100",
    stringsAsFactors = FALSE
  )
  current <- baseline[c(
    "engine", "scenario", "contrast_id", "metric", "estimate", "lower", "upper"
  )]
  current$estimate <- 0.88
  current$lower <- 0.84
  current$upper <- 0.92

  testthat::expect_equal(
    compare_validation_baseline(current, baseline)$status,
    "pass"
  )
  current$upper <- 0.80
  testthat::expect_equal(
    compare_validation_baseline(current, baseline)$status,
    "fail"
  )
})

testthat::test_that("missing and new baseline keys are explicit", {
  baseline <- data.frame(
    engine = "linda", scenario = "cross_sectional_signal",
    contrast_id = "condition[treated-control]", metric = "coverage",
    estimate = 0.93, lower = 0.90, upper = 0.96,
    direction = "lower", tolerance = 0.05, source_run_id = "100",
    stringsAsFactors = FALSE
  )
  current <- data.frame(
    engine = "deseq", scenario = "cross_sectional_signal",
    contrast_id = "condition[treated-control]", metric = "power",
    estimate = 0.9, lower = 0.8, upper = 1,
    stringsAsFactors = FALSE
  )
  compared <- compare_validation_baseline(current, baseline)

  testthat::expect_setequal(compared$status, c("fail", "informational"))
  testthat::expect_match(compared$reason[compared$status == "fail"], "missing")
})
