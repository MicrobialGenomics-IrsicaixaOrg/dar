test_that("bundled recipes use the centralized genus-level model", {
  data(test_rec)
  data(test_prep_rec)
  expect_s4_class(test_rec, "Recipe")
  expect_s4_class(test_prep_rec, "PrepRecipe")
  expect_identical(get_model(test_rec), get_model(test_prep_rec))
  expect_identical(get_model(test_prep_rec)$tax_level, "Genus")
  expect_identical(get_model(test_prep_rec)$targets, "RiskGroup2")
  expect_identical(get_model(test_prep_rec)$reference, c(RiskGroup2 = "hts"))
  expect_identical(steps_ids(test_rec), steps_ids(test_prep_rec))

  execution <- get_execution(test_prep_rec)
  expect_identical(execution$executed_steps, steps_ids(test_prep_rec, "da"))
  expect_length(execution$dropped_samples, 0L)
  expect_equal(nrow(execution$skipped_steps), 0L)
  expect_identical(execution$contrasts, dar:::resolve_model(test_prep_rec)$contrast_plan)
})

test_that("saved method results represent aggregated genera, not original features", {
  data(test_prep_rec)
  phy <- get_phy(test_prep_rec)
  expected_phy <- phyloseq::tax_glom(phy, taxrank = "Genus", NArm = FALSE)
  expected_ids <- phyloseq::taxa_names(expected_phy)
  expected_taxonomy <- as(phyloseq::tax_table(expected_phy), "matrix")
  results <- tidy_results(test_prep_rec)
  expect_lt(length(expected_ids), phyloseq::ntaxa(phy))
  expect_identical(
    anyDuplicated(results[c("taxa_id", "contrast_id", "step_id")]), 0L
  )
  expect_identical(results$taxa, unname(expected_taxonomy[results$taxa_id, "Genus"]))
  expect_false(any(results$contrast_type == "legacy_pairwise"))
  expect_false(anyNA(results$numerator))
  expect_false(anyNA(results$denominator))

  for (id in steps_ids(test_prep_rec, "da")) {
    method_results <- results[results$step_id == id, , drop = FALSE]
    expect_setequal(method_results$contrast_id, get_execution(test_prep_rec)$contrasts$contrast_id)
    for (contrast in unique(method_results$contrast_id)) {
      observed <- method_results[method_results$contrast_id == contrast, ]
      expect_setequal(observed$taxa_id, expected_ids)
      expect_equal(nrow(observed), length(expected_ids))
    }
  }

  bacteroides_id <- rownames(expected_taxonomy)[expected_taxonomy[, "Genus"] %in% "Bacteroides"]
  expect_length(bacteroides_id, 1L)
  original_taxonomy <- as(phyloseq::tax_table(phy), "matrix")
  members <- rownames(original_taxonomy)[original_taxonomy[, "Genus"] %in% "Bacteroides"]
  expect_gt(length(members), 1L)
  counts <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) counts <- t(counts)
  abundances <- tidy_abundance(test_prep_rec, taxa_ids = bacteroides_id, metadata = character())
  expect_equal(abundances$abundance, unname(colSums(counts[members, , drop = FALSE])))
  expect_identical(unique(abundances$taxa), "Bacteroides")

  for (id in steps_ids(test_prep_rec, "da")) {
    plotted <- effect_plt(
      test_prep_rec, taxa_ids = bacteroides_id, steps = id,
      contrast_id = results$contrast_id[[1]]
    )
    expect_s3_class(plotted, "ggplot")
    expect_equal(nrow(plotted$data), 1L)
    expect_identical(plotted$data$taxa, "Bacteroides")
  }
})

test_that("bundled recipe exports match the current configuration", {
  data(metaHIV_phy)
  data(test_rec)
  data(test_prep_rec)
  document <- import_steps(recipe(metaHIV_phy), system.file("extdata", "test.json", package = "dar"))
  expect_identical(
    dar:::recipe_serialization_document(document),
    dar:::recipe_serialization_document(test_rec)
  )
  with_bake <- import_steps(test_prep_rec, system.file("extdata", "test_bake.json", package = "dar"))
  expect_identical(get_model(with_bake), get_model(test_prep_rec))
  expect_length(with_bake@bakes, 1L)
  expect_identical(with_bake@bakes[[1]]$id, "bake__bundled")
})
