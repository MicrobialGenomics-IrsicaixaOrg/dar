testthat::test_that("CLI arguments and selections are validated", {
  parsed <- parse_validation_arguments(c(
    "--profile=full", "--engine=deseq,ancom", "--strict=false"
  ))
  testthat::expect_equal(parsed$profile, "full")
  testthat::expect_equal(parsed$engine, "deseq,ancom")
  testthat::expect_false(parsed$strict)
  testthat::expect_equal(
    validation_selection("deseq,ancom", names(validation_engine_registry())),
    c("deseq", "ancom")
  )
  testthat::expect_error(validation_selection("unknown", "deseq"), "Unknown")
})

testthat::test_that("simulation-only CLI execution writes complete artifacts", {
  output <- tempfile("dar-validation-cli-")
  outcome <- run_scientific_validation(
    profile_name = "quick",
    engines = "none",
    scenarios = "cross_sectional_null",
    output_dir = output,
    strict = TRUE
  )

  testthat::expect_equal(outcome$status, 0L)
  testthat::expect_equal(nrow(outcome$manifests), 2L)
  testthat::expect_true(direction_reversal_self_test())
  testthat::expect_true(file.exists(file.path(output, "validation-results.rds")))
})
