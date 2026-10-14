make_single_rank_analysis_phy <- function(two_ranks = FALSE) {
  taxa_ids <- c(
    "ASV_high", "ASV_low", "ASV_missing", "ASV_literal_na",
    "ASV_blank_empty_1", "ASV_blank_empty_2", "ASV_blank_space",
    "ASV_blank_tab", "ASV_single"
  )
  sample_ids <- paste0("sample_", seq_len(4))
  counts <- rbind(
    ASV_high = c(20, 18, 22, 24),
    ASV_low = c(5, 4, 6, 5),
    ASV_missing = c(15, 14, 16, 15),
    ASV_literal_na = c(2, 3, 2, 3),
    ASV_blank_empty_1 = c(8, 7, 8, 7),
    ASV_blank_empty_2 = c(3, 2, 3, 2),
    ASV_blank_space = c(4, 4, 4, 4),
    ASV_blank_tab = c(6, 6, 6, 6),
    ASV_single = c(9, 10, 9, 10)
  )
  colnames(counts) <- sample_ids
  labels <- c("same", "same", NA, "NA", "", "", " ", "\t", "solo")
  taxonomy <- matrix(
    labels, ncol = 1L,
    dimnames = list(taxa_ids, "Species")
  )
  if (two_ranks) {
    taxonomy <- cbind(taxonomy, Reference_rank = "constant")
  }
  tree <- ape::read.tree(
    text = paste0("(", paste(taxa_ids, collapse = ","), ");")
  )
  sequences <- Biostrings::DNAStringSet(stats::setNames(
    c("ACGT", "ACGA", "ACGC", "ACGG", "ATGT", "ATGA", "ATGC", "ATGG", "AAGT"),
    taxa_ids
  ))

  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(data.frame(
      condition = factor(c("control", "control", "treated", "treated")),
      row.names = sample_ids
    )),
    phyloseq::phy_tree(tree),
    sequences
  )
}

test_that("single-rank aggregation matches native NArm false semantics", {
  single_rank <- make_single_rank_analysis_phy()
  reference <- make_single_rank_analysis_phy(two_ranks = TRUE) |>
    phyloseq::tax_glom(taxrank = "Species", NArm = FALSE)
  rec <- recipe(single_rank) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")

  observed <- dar:::analysis_phyloseq(rec)
  expected_ids <- c(
    "ASV_high", "ASV_missing", "ASV_blank_empty_1",
    "ASV_blank_empty_2", "ASV_blank_space", "ASV_blank_tab", "ASV_single"
  )

  expect_identical(phyloseq::taxa_names(observed), expected_ids)
  expect_identical(
    phyloseq::taxa_names(observed), phyloseq::taxa_names(reference)
  )
  expect_equal(
    as(phyloseq::otu_table(observed), "matrix"),
    as(phyloseq::otu_table(reference), "matrix")
  )
  expect_identical(
    as(phyloseq::tax_table(observed), "matrix")[, "Species"],
    as(phyloseq::tax_table(reference), "matrix")[, "Species"]
  )
  expect_true(is.na(
    as(phyloseq::tax_table(observed), "matrix")["ASV_missing", "Species"]
  ))
  expect_identical(
    phyloseq::phy_tree(observed), phyloseq::phy_tree(reference)
  )
  expect_identical(
    phyloseq::refseq(observed), phyloseq::refseq(reference)
  )
  expect_identical(
    phyloseq::sample_data(observed), phyloseq::sample_data(reference)
  )
  expect_identical(get_phy(rec), single_rank)
})

test_that("single-rank aggregation is a no-op for one taxon", {
  phy <- suppressWarnings(phyloseq::prune_taxa(
    "ASV_single", make_single_rank_analysis_phy()
  ))
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")

  expect_identical(dar:::analysis_phyloseq(rec), phy)
})
