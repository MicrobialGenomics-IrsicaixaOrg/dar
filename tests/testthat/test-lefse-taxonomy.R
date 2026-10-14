make_lefse_taxonomy_phy <- function(taxa_are_rows = TRUE) {
  sample_ids <- paste0("sample_", seq_len(8))
  condition <- factor(
    rep(c("control", "treated"), each = 4),
    levels = c("control", "treated")
  )
  offsets <- rep(0:3, times = 2)
  counts <- rbind(
    "ASV one" = c(rep(20, 4), rep(80, 4)) + offsets,
    "ASV|two" = c(rep(80, 4), rep(20, 4)) + offsets,
    "ASV.three" = c(rep(30, 4), rep(55, 4)) + offsets,
    "ASV`four" = c(rep(5, 4), rep(10, 4)) + offsets,
    "ÅSV-five" = 40 + offsets,
    "ASV:six" = c(rep(60, 4), rep(25, 4)) + offsets,
    "merge-a" = c(rep(12, 4), rep(24, 4)) + offsets,
    "merge-b" = c(rep(8, 4), rep(16, 4)) + offsets
  )
  colnames(counts) <- sample_ids

  taxonomy <- cbind(
    Kingdom = rep("Bacteria", nrow(counts)),
    Genus = paste0("Genus_", seq_len(nrow(counts))),
    Species = paste0("Species_", seq_len(nrow(counts)))
  )
  rownames(taxonomy) <- rownames(counts)
  taxonomy[1:2, "Species"] <- "Shared species"
  taxonomy[8, c("Genus", "Species")] <-
    taxonomy[7, c("Genus", "Species")]

  phyloseq::phyloseq(
    phyloseq::otu_table(
      if (taxa_are_rows) counts else t(counts),
      taxa_are_rows = taxa_are_rows
    ),
    phyloseq::sample_data(data.frame(
      condition = condition,
      row.names = sample_ids
    )),
    phyloseq::tax_table(taxonomy)
  )
}

test_that("modeled LEfSe preserves canonical IDs across taxonomy and orientation", {
  skip_if_not_installed("lefser")

  results <- lapply(c(TRUE, FALSE), function(taxa_are_rows) {
    rec <- recipe(make_lefse_taxonomy_phy(taxa_are_rows)) |>
      add_model(~ condition, targets = "condition", tax_level = "Species")
    analytic_phy <- dar:::model_phyloseq(rec, rarefy = FALSE)
    expected_ids <- phyloseq::taxa_names(analytic_phy)
    expected_taxa <- as(
      phyloseq::tax_table(analytic_phy), "matrix"
    )[, "Species"]

    result <- withr::with_seed(
      130,
      suppressWarnings(dar:::run_lefse_model(
        rec, 1, 1, 0, 1L, TRUE, FALSE
      ))
    ) |>
      dar:::flatten_model_output()

    expect_length(expected_ids, 7L)
    expect_setequal(result$taxa_id, expected_ids)
    expect_identical(
      anyDuplicated(dplyr::select(result, "taxa_id", "contrast_id")),
      0L
    )
    expect_equal(result$taxa, unname(expected_taxa[result$taxa_id]))
    expect_equal(sum(result$taxa == "Shared species"), 2L)
    special_ids <- c(
      "ASV one", "ASV|two", "ASV.three", "ASV`four", "ÅSV-five", "ASV:six"
    )
    expect_true(all(special_ids %in% result$taxa_id))

    dplyr::arrange(result, .data$taxa_id)
  })

  expect_equal(results[[2]], results[[1]], tolerance = 1e-12)
})
