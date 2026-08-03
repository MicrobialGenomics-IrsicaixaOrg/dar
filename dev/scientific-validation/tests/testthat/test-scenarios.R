testthat::test_that("scenario seeds and simulations are deterministic", {
  first <- simulate_validation_scenario("cross_sectional_signal", 1L)
  second <- simulate_validation_scenario("cross_sectional_signal", 1L)
  third <- simulate_validation_scenario("cross_sectional_signal", 2L)

  testthat::expect_identical(
    as(phyloseq::otu_table(first$phy), "matrix"),
    as(phyloseq::otu_table(second$phy), "matrix")
  )
  testthat::expect_identical(first$truth, second$truth)
  testthat::expect_false(identical(
    as(phyloseq::otu_table(first$phy), "matrix"),
    as(phyloseq::otu_table(third$phy), "matrix")
  ))
})

testthat::test_that("simulated counts and truth keys are valid", {
  simulation <- simulate_validation_scenario("longitudinal_interaction", 1L)
  counts <- as(phyloseq::otu_table(simulation$phy), "matrix")
  planned <- dar:::resolve_model(simulation$recipe)$contrast_plan$contrast_id

  testthat::expect_true(all(is.finite(counts)))
  testthat::expect_true(all(counts >= 0))
  testthat::expect_true(all(counts == floor(counts)))
  testthat::expect_true(all(colSums(counts) > 0))
  testthat::expect_setequal(unique(simulation$truth$contrast_id), planned)
  testthat::expect_identical(
    anyDuplicated(simulation$truth[c("taxa_id", "contrast_id")]),
    0L
  )
  testthat::expect_true(any(simulation$truth$expected_direction > 0))
  testthat::expect_true(any(simulation$truth$expected_direction < 0))
  testthat::expect_true(any(simulation$truth$is_null))
})

testthat::test_that("the global null contains no target effects", {
  simulation <- simulate_validation_scenario("cross_sectional_null", 1L)
  testthat::expect_true(all(simulation$truth$is_null))
  testthat::expect_true(all(simulation$truth$expected_direction == 0))
})

testthat::test_that("zero inflation increases sparsity without changing truth", {
  dense <- simulate_validation_scenario("cross_sectional_signal", 1L)
  sparse <- simulate_validation_scenario("sparse_zero_inflated", 1L)
  dense_counts <- as(phyloseq::otu_table(dense$phy), "matrix")
  sparse_counts <- as(phyloseq::otu_table(sparse$phy), "matrix")

  testthat::expect_gt(mean(sparse_counts == 0), mean(dense_counts == 0))
  testthat::expect_equal(
    sparse$truth[c("taxa_id", "contrast_id", "true_effect")],
    dense$truth[c("taxa_id", "contrast_id", "true_effect")]
  )
})

testthat::test_that("missing and unbalanced data retain an estimable cohort", {
  simulation <- simulate_validation_scenario("missing_unbalanced", 1L)
  resolved <- dar:::resolve_model(simulation$recipe)
  cells <- table(resolved$data$condition, resolved$data$time)

  testthat::expect_equal(length(resolved$dropped_samples), 2L)
  testthat::expect_true(all(cells > 0L))
  testthat::expect_equal(qr(resolved$matrix)$rank, ncol(resolved$matrix))
  testthat::expect_equal(
    simulation$manifest$complete_samples,
    phyloseq::nsamples(simulation$phy) - 2L
  )
})

testthat::test_that("quick and full profiles cover every registered scenario", {
  scenarios <- names(validation_scenario_registry())
  testthat::expect_setequal(names(validation_profile("quick")$replicates), scenarios)
  testthat::expect_setequal(names(validation_profile("full")$replicates), scenarios)
  testthat::expect_true(all(validation_profile("quick")$replicates == 2L))
  testthat::expect_gt(
    validation_profile("full")$replicates[["cross_sectional_null"]],
    validation_profile("quick")$replicates[["cross_sectional_null"]]
  )
})
