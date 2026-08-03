testthat::test_that("engine capabilities match centralized model support", {
  simple <- validation_scenario("cross_sectional_signal")
  complex <- validation_scenario("longitudinal_random_effect")

  testthat::expect_true(validation_engine_compatible(simple, "wilcox")$compatible)
  testthat::expect_false(validation_engine_compatible(complex, "wilcox")$compatible)
  testthat::expect_true(validation_engine_compatible(complex, "ancom")$compatible)
  testthat::expect_true(validation_engine_compatible(complex, "maaslin")$compatible)
  testthat::expect_false(validation_engine_compatible(complex, "deseq")$compatible)
})

testthat::test_that("ANCOM installs every runtime validation dependency", {
  testthat::expect_contains(
    validation_engine("ancom")$packages,
    "microbiome"
  )
})

testthat::test_that("normalization preserves unique scientific keys", {
  simulation <- simulate_validation_scenario("cross_sectional_signal", 1L)
  raw <- data.frame(
    taxa_id = c("taxon_01", "taxon_02"),
    contrast_id = "condition[treated-control]",
    effect_size = c(1, -1),
    effect_metric = "log2_fold_change",
    adj_p_value = c(0.01, 0.02),
    lfcSE = c(0.2, 0.3)
  )
  normalized <- normalize_validation_result(
    raw[c("taxa_id", "contrast_id", "effect_size", "effect_metric", "adj_p_value")],
    "deseq", simulation, raw_result = raw
  )

  testthat::expect_named(normalized, c(
    "engine", "scenario", "replicate", "taxa_id", "contrast_id", "effect",
    "padj", "std_error", "effect_metric", "truth_multiplier", "comparable"
  ))
  testthat::expect_equal(normalized$std_error, raw$lfcSE)
  duplicate <- rbind(raw, raw[1, ])
  testthat::expect_error(
    normalize_validation_result(
      duplicate[c("taxa_id", "contrast_id", "effect_size", "effect_metric", "adj_p_value")],
      "deseq", simulation
    ),
    "duplicated"
  )
})

testthat::test_that("incompatible engines are expected skips", {
  simulation <- simulate_validation_scenario("longitudinal_random_effect", 1L)
  outcome <- run_validation_case(simulation, "wilcox", validation_profile("quick"))

  testthat::expect_equal(outcome$run$status, "expected_skip")
  testthat::expect_match(outcome$run$reason, "confounders")
  testthat::expect_equal(nrow(outcome$results), 0L)
})
