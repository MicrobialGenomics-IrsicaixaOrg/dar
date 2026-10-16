make_legacy_maaslin_phy <- function(taxa_are_rows = TRUE) {
  counts <- rbind(
    ASV_1 = c(40, 50, 60, 70, 80, 90),
    ASV_2 = c(10, 12, 14, 16, 18, 20),
    ASV_3 = c(25, 30, 35, 40, 45, 50),
    ASV_4 = c(15, 20, 25, 30, 35, 40)
  )
  colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
  taxonomy <- cbind(
    Kingdom = "Bacteria",
    Family = c("Family_A", "Family_A", "Family_B", "Family_C"),
    Genus = c("Shared_genus", "Shared_genus", "Shared_genus", "Other_genus"),
    Species = paste0("Species_", seq_len(nrow(counts)))
  )
  rownames(taxonomy) <- rownames(counts)
  metadata <- data.frame(
    condition = rep(c("control", "rescue", "treated"), each = 2L),
    row.names = colnames(counts)
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(
      if (taxa_are_rows) counts else t(counts),
      taxa_are_rows = taxa_are_rows
    ),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(metadata)
  )
}

capture_legacy_maaslin_inputs <- function() {
  calls <- list()
  testthat::local_mocked_bindings(
    maaslin3_quietly = function(input_data, input_metadata, ...) {
      calls[[length(calls) + 1L]] <<- list(
        counts = input_data, metadata = input_metadata
      )
      list(result = list(fit_data_abundance = list(results = data.frame(
        feature = rownames(input_data),
        coef = seq_len(nrow(input_data)),
        stderr = 0.1,
        pval_individual = 0.01,
        qval_individual = 0.02
      ))))
    },
    .package = "dar",
    .env = parent.frame()
  )
  function() calls
}

run_legacy_maaslin_counts <- function(rec, rarefy = FALSE) {
  dar:::run_maaslin(
    rec = rec,
    min_abundance = 0,
    min_prevalence = 0.1,
    min_variance = 0,
    normalization = "TSS",
    transform = "LOG",
    max_significance = 0.1,
    random_effects = NULL,
    correction = "BH",
    standardize = TRUE,
    reference = NULL,
    median_comparison_abundance = TRUE,
    rarefy = rarefy,
    id = "maaslin__legacy_aggregation"
  )
}

test_that("legacy MaAsLin fits aggregated counts in either OTU orientation", {
  calls <- capture_legacy_maaslin_inputs()
  for (taxa_are_rows in c(TRUE, FALSE)) {
    phy <- make_legacy_maaslin_phy(taxa_are_rows)
    rec <- suppressWarnings(recipe(phy, "condition", "Genus"))
    original <- get_phy(rec)
    counts <- as(phyloseq::otu_table(phy), "matrix")
    if (!taxa_are_rows) counts <- t(counts)
    expected <- rbind(
      ASV_1 = counts["ASV_1", ] + counts["ASV_2", ],
      counts[c("ASV_3", "ASV_4"), , drop = FALSE]
    )

    previous_calls <- length(calls())
    result <- run_legacy_maaslin_counts(rec)$condition
    observed <- tail(calls(), length(calls()) - previous_calls)

    expect_length(observed, 3L)
    for (call in observed) {
      expect_identical(colnames(call$counts), rownames(call$metadata))
      expect_equal(
        call$counts,
        expected[, rownames(call$metadata), drop = FALSE]
      )
      expect_equal(ncol(call$counts), 4L)
      expect_length(unique(call$metadata$condition), 2L)
    }
    expect_equal(nrow(result), 9L)
    expect_setequal(result$taxa_id, rownames(expected))
    expect_false("ASV_2" %in% result$taxa_id)
    expect_equal(sum(result$taxa == "Shared_genus"), 6L)
    expect_equal(
      result$taxa,
      unname(c(ASV_1 = "Shared_genus", ASV_3 = "Shared_genus",
               ASV_4 = "Other_genus")[result$taxa_id])
    )
    expect_identical(get_phy(rec), original)
    expect_identical(get_phy(rec), phy)
  }
})

test_that("legacy MaAsLin fits rarefied rather than stored counts", {
  withr::local_seed(130)
  calls <- capture_legacy_maaslin_inputs()
  for (taxa_are_rows in c(TRUE, FALSE)) {
    phy <- make_legacy_maaslin_phy(taxa_are_rows)
    rec <- suppressWarnings(recipe(phy, "condition", "Genus"))
    expected_phy <- phyloseq::rarefy_even_depth(
      phy, rngseed = 1234, verbose = FALSE
    ) |>
      phyloseq::tax_glom(taxrank = "Genus", NArm = FALSE)
    expected <- as(phyloseq::otu_table(expected_phy), "matrix")
    if (!taxa_are_rows) expected <- t(expected)
    target_depth <- min(phyloseq::sample_sums(phy))

    previous_calls <- length(calls())
    result <- run_legacy_maaslin_counts(rec, rarefy = TRUE)$condition
    observed <- tail(calls(), length(calls()) - previous_calls)

    expect_length(observed, 3L)
    for (call in observed) {
      expect_equal(
        call$counts,
        expected[, rownames(call$metadata), drop = FALSE]
      )
      expect_equal(unname(colSums(call$counts)), rep(target_depth, 4L))
      expect_identical(colnames(call$counts), rownames(call$metadata))
    }
    expect_setequal(result$taxa_id, rownames(expected))
    expect_false(anyNA(result$taxa))
    expect_identical(get_phy(rec), phy)
  }
})
