make_analysis_plot_phy <- function(taxa_are_rows = TRUE) {
  counts <- rbind(
    asv_a = c(40, 30, 50, 20, 35, 45),
    asv_b = c(10, 15, 5, 20, 12, 8),
    asv_c = c(50, 40, 35, 60, 25, 55)
  )
  colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
  taxonomy <- cbind(
    Kingdom = "Bacteria",
    Phylum = c("Phylum_A", "Phylum_A", "Phylum_C"),
    Genus = c("Genus_A", "Genus_A", "Genus_C"),
    Species = c("Species_a", "Species_b", "Species_c")
  )
  rownames(taxonomy) <- rownames(counts)
  metadata <- data.frame(
    condition = factor(rep(c("control", "treated"), each = 3)),
    time = c(0, 1, 2, 0, 1, NA),
    batch = rep(c("A", "B"), 3),
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

analysis_plot_expected_counts <- function() {
  counts <- rbind(
    "asv_a|Genus_A" = c(50, 45, 55, 40, 47, 53),
    "asv_c|Genus_C" = c(50, 40, 35, 60, 25, 55)
  )
  colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
  counts
}

analysis_plot_data_matrix <- function(data) {
  data |>
    dplyr::select("taxa", "sample_id", "value") |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "value") |>
    data.frame(row.names = 1, check.names = FALSE) |>
    as.matrix()
}

test_that("annotated plot counts sum all ASVs at the modeled rank", {
  for (taxa_are_rows in c(TRUE, FALSE)) {
    rec <- recipe(make_analysis_plot_phy(taxa_are_rows)) |>
      add_model(~ condition, targets = "condition", tax_level = "Genus")
    original_phy <- get_phy(rec)
    original_counts <- otu_table(rec)
    original_taxonomy <- tax_table(rec)

    annotated <- dar:::.annotated_counts(rec) |>
      tidyr::unite("taxa", c("taxa_id", "taxa"), sep = "|")
    expect_equal(
      analysis_plot_data_matrix(annotated),
      analysis_plot_expected_counts()
    )
    expect_identical(get_phy(rec), original_phy)
    expect_identical(otu_table(rec), original_counts)
    expect_identical(tax_table(rec), original_taxonomy)
    expect_equal(nrow(original_counts), 3L)
  }
})

test_that("boxplots aggregate counts before abundance transformations", {
  skip_if_not_installed("microbiome")
  expected_counts <- analysis_plot_expected_counts()
  expected_compositional <- sweep(
    expected_counts, 2, colSums(expected_counts), "/"
  )
  expected_clr <- sweep(
    log(expected_counts), 2, colMeans(log(expected_counts)), "-"
  )
  expected <- list(
    raw = expected_counts,
    compositional = expected_compositional,
    clr = expected_clr
  )

  for (taxa_are_rows in c(TRUE, FALSE)) {
    rec <- recipe(make_analysis_plot_phy(taxa_are_rows)) |>
      add_model(~ condition, targets = "condition", tax_level = "Genus")
    prepared <- dar:::prep_recipe(rec, results = list(), bakes = list())
    original_phy <- get_phy(prepared)
    original_counts <- otu_table(prepared)
    for (transformation in names(expected)) {
      plot <- dar:::.abundance_boxplot(
        prepared, taxa_ids = c("asv_a", "asv_c"),
        transform = if (transformation == "raw") NULL else transformation,
        scale = 1, target = "condition", contrast_id = NULL
      )
      expect_equal(
        analysis_plot_data_matrix(plot$data),
        expected[[transformation]],
        info = paste("Transformation:", transformation)
      )
    }
    expect_identical(get_phy(prepared), original_phy)
    expect_identical(otu_table(prepared), original_counts)
  }
})

test_that("heatmaps use the same aggregated transformed abundance matrix", {
  skip_if_not_installed("microbiome")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("gplots")
  withr::defer(ComplexHeatmap::ht_opt(RESET = TRUE))
  rec <- recipe(make_analysis_plot_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Genus")
  prepared <- dar:::prep_recipe(rec, results = list(), bakes = list())
  original_phy <- get_phy(prepared)
  expected_counts <- analysis_plot_expected_counts()
  expected <- list(
    raw = expected_counts,
    compositional = sweep(expected_counts, 2, colSums(expected_counts), "/"),
    clr = sweep(log(expected_counts), 2, colMeans(log(expected_counts)), "-")
  )

  for (transformation in names(expected)) {
    heatmap <- dar:::.abundance_heatmap(
      prepared, taxa_ids = c("asv_a", "asv_c"),
      transform = if (transformation == "raw") NULL else transformation,
      scale = 1, top_n = 20, target = "condition", contrast_id = NULL
    )
    expect_equal(
      heatmap@matrix, expected[[transformation]],
      info = paste("Transformation:", transformation)
    )
  }
  expect_identical(get_phy(prepared), original_phy)
})

test_that("plot annotations follow the prepared model cohort", {
  rec <- recipe(make_analysis_plot_phy()) |>
    add_model(
      ~ condition + time, targets = "condition", tax_level = "Genus",
      time = "time", na_action = "drop"
    )
  original_phy <- get_phy(rec)
  modeled <- dar:::apply_model_to_recipe(rec)$rec
  prepared <- dar:::prep_recipe(modeled, results = list(), bakes = list())
  original_counts <- otu_table(prepared)
  original_metadata <- sample_data(prepared)
  annotated <- dar:::.annotated_counts(prepared)

  expect_setequal(annotated$sample_id, paste0("sample_", seq_len(5)))
  expect_false(anyNA(annotated$condition))
  expect_identical(
    annotated$condition,
    original_metadata$condition[
      match(annotated$sample_id, original_metadata$sample_id)
    ]
  )
  expect_equal(nrow(annotated), 2L * 5L)
  expect_named(original_metadata, c("sample_id", "condition", "time", "batch"))
  expect_identical(otu_table(prepared), original_counts)
  expect_identical(sample_data(prepared), original_metadata)
  expect_identical(get_phy(rec), original_phy)
  expect_equal(nrow(otu_table(rec)), 3L)
  expect_equal(ncol(otu_table(rec)) - 1L, 6L)
})
