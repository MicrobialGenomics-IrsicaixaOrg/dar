make_analysis_universe_phy <- function(taxa_are_rows = TRUE) {
  phy <- make_multilevel_phy()
  phyloseq::taxa_names(phy) <- paste0("ASV_", seq_len(phyloseq::ntaxa(phy)))
  counts <- as(phyloseq::otu_table(phy), "matrix")
  counts[1, ] <- counts[1, ] * 2
  phyloseq::otu_table(phy) <- phyloseq::otu_table(
    if (taxa_are_rows) counts else t(counts),
    taxa_are_rows = taxa_are_rows
  )
  taxonomy <- cbind(
    Kingdom = "Bacteria",
    Genus = c("Merged", "Merged", "Genus_3", "Genus_4", "Genus_5", "Genus_6"),
    Species = c("Merged_species", "Merged_species", "Shared_species",
                "Shared_species", "Species_5", NA_character_)
  )
  rownames(taxonomy) <- rownames(counts)
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  phy
}

test_that("analysis counts and taxonomy share a stable aggregated universe", {
  withr::local_seed(130)
  expected_ids <- c("ASV_1", "ASV_3", "ASV_4", "ASV_5", "ASV_6")
  for (taxa_are_rows in c(TRUE, FALSE)) {
    phy <- make_analysis_universe_phy(taxa_are_rows)
    rec <- recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Species")
    raw <- dar:::recipe_count_matrix(rec)
    expected <- rbind(ASV_1 = raw[1, ] + raw[2, ], raw[3:6, , drop = FALSE])

    expect_equal(dar:::analysis_count_matrix(rec), expected)
    taxonomy <- dar:::analysis_tax_table(rec)
    expect_identical(taxonomy$taxa_id, expected_ids)
    expect_equal(sum(taxonomy$taxa == "Shared_species", na.rm = TRUE), 2L)
    expect_true(is.na(taxonomy$taxa[taxonomy$taxa_id == "ASV_6"]))
    expect_identical(
      phyloseq::taxa_names(dar:::model_phyloseq(rec)), expected_ids
    )
    expect_setequal(
      phyloseq::taxa_names(dar:::model_phyloseq(rec, rarefy = TRUE)),
      expected_ids
    )
    expect_identical(get_phy(rec), phy)
    expect_identical(otu_table(rec)$taxa_id, phyloseq::taxa_names(phy))
    expect_identical(
      tax_table(rec)$Species,
      unname(as(phyloseq::tax_table(phy), "matrix")[, "Species"])
    )
  }
})

test_that("modeled consensus contains only hypotheses actually aggregated", {
  phy <- make_analysis_universe_phy()
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(id = "wilcox__universe")
  prepared <- prep(rec, parallel = FALSE)
  expected_ids <- c("ASV_1", "ASV_3", "ASV_4", "ASV_5", "ASV_6")

  results <- tidy_results(prepared)
  expect_setequal(results$taxa_id, expected_ids)
  expect_equal(nrow(results), 15L)
  intersections <- intersection_df(prepared)
  expect_setequal(intersections$taxa_id, expected_ids)
  expect_equal(nrow(intersections), 15L)
  expect_false("ASV_2" %in% intersections$taxa_id)
  expect_equal(
    nrow(intersection_df(prepared, tidy = TRUE)), nrow(intersections)
  )
  expect_identical(get_phy(prepared), phy)
  expect_identical(otu_table(prepared), otu_table(rec))

  repeated <- prep(prepared, parallel = FALSE, force = TRUE)
  expect_identical(tidy_results(repeated), results)
  expect_identical(intersection_df(repeated), intersections)
  expect_identical(get_phy(repeated), phy)
})

test_that("consensus keeps homonymous analysis taxa as distinct hypotheses", {
  rec <- recipe(make_analysis_universe_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step("wilcox", id = "wilcox__one"))
  rec <- dar:::add_step(rec, dar:::step("lefse", id = "lefse__two"))
  contrast <- dar:::resolve_model(rec)$contrast_plan[1, ]
  rows <- tibble::tibble(
    taxa_id = c("ASV_3", "ASV_4"), taxa = "Shared_species",
    contrast_id = contrast$contrast_id, comparison = contrast$comparison,
    contrast_type = contrast$contrast_type, var = contrast$var,
    effect = c(1, -1), padj = 0.01, signif = TRUE
  )
  prepared <- dar:::prep_recipe(
    rec,
    list(wilcox__one = list(model = rows), lefse__two = list(model = rows)),
    list()
  )

  consensus <- cool(bake(prepared))
  expect_setequal(consensus$taxa_id, c("ASV_3", "ASV_4"))
  expect_equal(consensus$taxa, rep("Shared_species", 2))
  expect_equal(consensus$method_count, c(2, 2))
  intersections <- find_intersections(prepared)
  expect_setequal(intersections$taxa_id, c("ASV_3", "ASV_4"))
  expect_equal(intersections$sum_methods, c(2, 2))
})

test_that("model results cannot invent or mislabel analysis taxon identities", {
  rec <- recipe(make_analysis_universe_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(id = "wilcox__universe")
  prepared <- prep(rec, parallel = FALSE)

  for (invalid_id in c("ASV_2", "ghost_taxon")) {
    invalid <- prepared
    invalid@results$wilcox__universe$model$taxa_id[1] <- invalid_id
    expect_error(
      tidy_results(invalid), "outside the analysis universe",
      class = "dar_error_invalid_result_contract"
    )
  }
  invalid <- prepared
  invalid@results$wilcox__universe$model$taxa[1] <- "Wrong_species"
  expect_error(
    tidy_results(invalid), "labels inconsistent",
    class = "dar_error_invalid_result_contract"
  )
  unlabeled <- prepared
  unlabeled@results$wilcox__universe$model$taxa <- NULL
  expect_identical(tidy_results(unlabeled)$taxa, tidy_results(prepared)$taxa)
})
