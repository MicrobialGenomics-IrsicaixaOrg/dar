make_analysis_data_phy <- function(taxa_are_rows = TRUE) {
  counts <- rbind(
    ASV_high = c(12, 18, 24, 30),
    ASV_low = c(3, 2, 1, 4),
    ASV_other_lineage = c(7, 8, 9, 10),
    ASV_missing = c(2, 1, 3, 2)
  )
  colnames(counts) <- c("sample_3", "sample_1", "sample_4", "sample_2")
  taxonomy <- cbind(
    Kingdom = "Bacteria",
    Genus = c("First", "First", "Second", "Third"),
    Species = c("Shared", "Shared", "Shared", NA_character_)
  )
  rownames(taxonomy) <- rownames(counts)
  metadata <- data.frame(
    condition = factor(c("treated", "control", "treated", "control"),
                       levels = c("control", "treated")),
    visit = c(1, 0, NA, 1),
    subject = c("p2", "p1", "p2", "p1"),
    Species = c("metadata one", "metadata two", "metadata three", "metadata four"),
    abundance = seq_len(4),
    taxa = letters[1:4],
    taxa_id = letters[5:8],
    row.names = colnames(counts),
    check.names = FALSE
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(if (taxa_are_rows) counts else t(counts),
                       taxa_are_rows = taxa_are_rows),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(metadata)
  )
}

test_that("tidy abundances use shared aggregation with full lineage and metadata", {
  for (taxa_are_rows in c(TRUE, FALSE)) {
    phy <- make_analysis_data_phy(taxa_are_rows)
    rec <- recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Species")
    original <- get_phy(rec)
    observed <- tidy_abundance(rec)
    expected_ids <- c("ASV_high", "ASV_other_lineage", "ASV_missing")
    expected_counts <- rbind(
      ASV_high = c(15, 20, 25, 34),
      ASV_other_lineage = c(7, 8, 9, 10),
      ASV_missing = c(2, 1, 3, 2)
    )
    expect_named(observed, c(
      "taxa_id", "taxa", "sample_id", "abundance",
      "taxonomy_Kingdom", "taxonomy_Genus", "taxonomy_Species",
      "metadata_condition", "metadata_visit", "metadata_subject",
      "metadata_Species", "metadata_abundance", "metadata_taxa",
      "metadata_taxa_id"
    ))
    expect_identical(observed$taxa_id, rep(expected_ids, each = 4))
    expect_identical(observed$sample_id, rep(phyloseq::sample_names(phy), 3))
    expect_equal(observed$abundance, as.numeric(t(expected_counts)))
    expect_equal(unique(observed$taxa_id[observed$taxa %in% "Shared"]),
                 c("ASV_high", "ASV_other_lineage"))
    expect_equal(unique(observed$taxonomy_Genus), c("First", "Second", "Third"))
    expect_true(all(is.na(observed$taxa[observed$taxa_id == "ASV_missing"])))
    expect_s3_class(observed$metadata_condition, "factor")
    expect_identical(levels(observed$metadata_condition), c("control", "treated"))
    expect_identical(observed$metadata_visit, rep(c(1, 0, NA, 1), 3))
    expect_identical(get_phy(rec), original)
  }
})

test_that("tidy abundances select taxa and metadata without changing denominators", {
  rec <- recipe(make_analysis_data_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  all_fractions <- tidy_abundance(rec, transform = "compositional")
  selected <- tidy_abundance(
    rec, taxa_ids = c("ASV_missing", "ASV_high"),
    transform = "compositional", metadata = c("subject", "visit")
  )
  expect_named(selected, c(
    "taxa_id", "taxa", "sample_id", "abundance", "taxonomy_Kingdom",
    "taxonomy_Genus", "taxonomy_Species", "metadata_subject", "metadata_visit"
  ))
  expect_identical(unique(selected$taxa_id), c("ASV_high", "ASV_missing"))
  expect_equal(selected$abundance, all_fractions$abundance[
    all_fractions$taxa_id %in% c("ASV_high", "ASV_missing")
  ])
  expect_equal(as.vector(tapply(all_fractions$abundance, all_fractions$sample_id, sum)),
               rep(1, 4))
  expect_equal(selected$abundance[1], 15 / 24)

  without_metadata <- tidy_abundance(rec, metadata = character(0))
  expect_false(any(startsWith(names(without_metadata), "metadata_")))
  empty <- tidy_abundance(rec, taxa_ids = character(0), metadata = "visit")
  expect_equal(nrow(empty), 0L)
  expect_type(empty$taxa_id, "character")
  expect_type(empty$abundance, "double")
  expect_named(empty, c(
    "taxa_id", "taxa", "sample_id", "abundance", "taxonomy_Kingdom",
    "taxonomy_Genus", "taxonomy_Species", "metadata_visit"
  ))
  expect_error(tidy_abundance(rec, taxa_ids = "ASV_low"), "shared analysis view",
               class = "dar_error_invalid_taxa_id")
  expect_error(tidy_abundance(rec, metadata = "missing_column"), "metadata column",
               class = "dar_error_invalid_input")
  expect_error(tidy_abundance(rec, taxa_ids = c("ASV_high", "ASV_high")))
  expect_error(tidy_abundance(rec, transform = "unknown"))
  expect_error(tidy_abundance(NULL), class = "dar_error_invalid_input")
})

test_that("tidy abundances work without a model or optional annotations", {
  phy <- make_analysis_data_phy()
  rec <- recipe(phy)
  observed <- tidy_abundance(rec)
  expect_identical(unique(observed$taxa_id), phyloseq::taxa_names(phy))
  expect_identical(observed$taxa, observed$taxa_id)

  metadata <- as(phyloseq::sample_data(phy), "data.frame")
  metadata$sample_id <- paste0("metadata_", seq_len(nrow(metadata)))
  phyloseq::sample_data(phy) <- phyloseq::sample_data(metadata)
  observed <- tidy_abundance(recipe(phy), metadata = "sample_id")
  expect_identical(observed$sample_id[seq_len(4)], phyloseq::sample_names(phy))
  expect_identical(observed$metadata_sample_id[seq_len(4)], metadata$sample_id)

  minimal_phy <- methods::new("phyloseq", otu_table = phyloseq::otu_table(phy))
  minimal <- tidy_abundance(recipe(minimal_phy))
  expect_named(minimal, c("taxa_id", "taxa", "sample_id", "abundance"))
  expect_equal(nrow(minimal), 16L)
  expect_error(tidy_abundance(recipe(minimal_phy), metadata = "subject"),
               class = "dar_error_invalid_input")

  legacy <- suppressWarnings(recipe(phy, var_info = "condition", tax_info = "Species"))
  expect_identical(unique(tidy_abundance(legacy)$taxa_id), phyloseq::taxa_names(phy))
  expect_equal(tidy_abundance(legacy)$taxa[1], "Shared")
})

test_that("compositional extraction handles zero libraries explicitly", {
  phy <- make_analysis_data_phy()
  counts <- as(phyloseq::otu_table(phy), "matrix")
  counts[, 2] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  rec <- recipe(phy)
  expect_warning(
    observed <- tidy_abundance(rec, transform = "compositional"),
    "sample_1", class = "dar_warning_zero_library"
  )
  expect_true(all(is.na(observed$abundance[observed$sample_id == "sample_1"])))
  expect_true(all(is.finite(observed$abundance[observed$sample_id != "sample_1"])))
  expect_equal(tidy_abundance(rec)$abundance[2], 0)

  counts[1, 1] <- -1
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  expect_error(tidy_abundance(recipe(phy), transform = "compositional"),
               "non-negative", class = "dar_error_invalid_input")
})

test_that("prepared abundance extraction uses stored preprocessed cohort only", {
  phy <- make_analysis_data_phy()
  prepared_phy <- phyloseq::prune_samples(c("sample_3", "sample_2"), phy)
  rec <- recipe(prepared_phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  prepared <- dar:::prep_recipe(rec, results = list(), bakes = list())
  expect_identical(tidy_abundance(prepared), tidy_abundance(rec))
  expect_setequal(unique(tidy_abundance(prepared)$sample_id),
                  c("sample_3", "sample_2"))
})

test_that("execution accessor returns recorded information without inference", {
  rec <- recipe(make_analysis_data_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  execution <- list(
    contrasts = tibble::tibble(contrast_id = "condition[treated-control]"),
    dropped_samples = "sample_4",
    executed_steps = "wilcox__one",
    skipped_steps = tibble::tibble(
      step_id = "lefse__two", compatible = FALSE, reason = "unsupported design"
    )
  )
  prepared <- dar:::prep_recipe(rec, results = list(), bakes = list(),
                               execution = execution)
  expect_identical(get_execution(prepared), execution)
  extracted <- get_execution(prepared)
  extracted$dropped_samples <- "changed"
  expect_identical(get_execution(prepared)$dropped_samples, "sample_4")

  legacy <- dar:::prep_recipe(rec, results = list(), bakes = list())
  expect_identical(get_execution(legacy), list(
    contrasts = NULL, dropped_samples = NULL, executed_steps = NULL,
    skipped_steps = NULL
  ))
  legacy@execution <- list(executed_steps = character())
  expect_identical(get_execution(legacy)$executed_steps, character())
  expect_null(get_execution(legacy)$dropped_samples)
  expect_error(get_execution(rec), class = "dar_error_invalid_input")
  expect_error(get_execution(list()), class = "dar_error_invalid_input")
})

test_that("public execution and abundance accessors follow actual prep cohorts", {
  rec <- recipe(make_analysis_data_phy()) |>
    add_model(
      ~ condition + visit, targets = "condition", tax_level = "Species",
      na_action = "drop"
    ) |>
    step_filter_by_prevalence(threshold = 0)
  original <- rec
  prepared <- prep(rec, parallel = FALSE)
  execution <- get_execution(prepared)
  expect_identical(execution$dropped_samples, "sample_4")
  expect_identical(execution$contrasts, dar:::resolve_model(prepared)$contrast_plan)
  expect_length(execution$executed_steps, 0L)
  expect_equal(nrow(execution$skipped_steps), 0L)
  observed <- tidy_abundance(prepared, metadata = c("subject", "visit"))
  expect_setequal(observed$sample_id, c("sample_3", "sample_1", "sample_2"))
  expect_false(anyNA(observed$metadata_visit))
  expect_identical(rec, original)
  expect_true("sample_4" %in% tidy_abundance(rec)$sample_id)
})
