testthat::test_that("CLI arguments and selections are validated", {
  parsed <- parse_validation_arguments(c(
    "--profile=full", "--engine=deseq,linda", "--strict=false"
  ))
  testthat::expect_equal(parsed$profile, "full")
  testthat::expect_equal(parsed$engine, "deseq,linda")
  testthat::expect_false(parsed$strict)
  testthat::expect_equal(
    validation_selection("deseq,linda", names(validation_engine_registry())),
    c("deseq", "linda")
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

testthat::test_that("the full runner cannot approve selectively returned signals", {
  testthat::local_mocked_bindings(
    run_validation_case = function(simulation, engine, profile) {
      config <- validation_engine(engine)
      truth <- simulation$truth[1, ]
      list(
        results = data.frame(
          engine = engine, scenario = truth$scenario, replicate = truth$replicate,
          taxa_id = truth$taxa_id, contrast_id = truth$contrast_id,
          effect = truth$true_effect * config$truth_multiplier,
          padj = 0.01, std_error = 0.1, effect_metric = config$effect_metric,
          truth_multiplier = config$truth_multiplier, comparable = config$comparable
        ),
        run = validation_run_record(engine, simulation, "success")
      )
    }
  )
  outcome <- run_scientific_validation(
    profile_name = "full", engines = "deseq", scenarios = "cross_sectional_signal",
    output_dir = tempfile("dar-validation-selective-"), strict = TRUE
  )
  metrics <- outcome$metrics[outcome$metrics$effect_metric != "execution", ]

  testthat::expect_true(all(metrics$n_expected == 60L))
  testthat::expect_true(all(metrics$n_returned == 1L))
  testthat::expect_true(all(metrics$n_missing == 59L))
  testthat::expect_equal(metrics$power, rep(1 / 8, nrow(metrics)))
  testthat::expect_true(outcome$gate_failures)
  testthat::expect_identical(outcome$status, 1L)
})
