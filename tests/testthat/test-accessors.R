make_accessor_phy <- function(taxa_are_rows = TRUE) {
  counts <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("taxon_c", "taxon_a", "taxon_b"),
      c("sample_3", "sample_1", "sample_4", "sample_2")
    )
  )
  metadata <- data.frame(
    condition = factor(
      c("treated", "control", "treated", "control"),
      levels = c("control", "treated")
    ),
    time = c(2, 0, NA, 1),
    batch = c("b", "a", "b", "a"),
    subject = c("s2", "s1", "s2", "s1"),
    row.names = colnames(counts),
    check.names = FALSE
  )
  taxonomy <- cbind(
    Phylum = c("Firmicutes", "Bacteroidota", NA),
    Species = c("Species c", "Species a", "Species b")
  )
  rownames(taxonomy) <- rownames(counts)

  phyloseq::phyloseq(
    phyloseq::otu_table(
      if (taxa_are_rows) counts else t(counts),
      taxa_are_rows = taxa_are_rows
    ),
    phyloseq::sample_data(metadata),
    phyloseq::tax_table(taxonomy)
  )
}

test_that("public accessors return complete model-independent components", {
  phy <- make_accessor_phy()
  model_free <- recipe(phy)
  modeled <- add_model(
    model_free,
    ~ condition + batch,
    targets = c("condition", "batch"),
    tax_level = "Species"
  )
  prepared <- dar:::prep_recipe(modeled, results = list(), bakes = list())

  for (rec in list(model_free, modeled, prepared)) {
    metadata <- sample_data(rec)
    taxonomy <- tax_table(rec)
    counts <- otu_table(rec)

    expect_named(
      metadata,
      c("sample_id", "condition", "time", "batch", "subject")
    )
    expect_identical(metadata$sample_id, phyloseq::sample_names(phy))
    expect_s3_class(metadata$condition, "factor")
    expect_identical(levels(metadata$condition), c("control", "treated"))
    expect_type(metadata$time, "double")
    expect_true(anyNA(metadata$time))

    expect_named(taxonomy, c("taxa_id", "Phylum", "Species"))
    expect_identical(taxonomy$taxa_id, phyloseq::taxa_names(phy))
    expect_true(anyNA(taxonomy$Phylum))

    expect_named(counts, c("taxa_id", phyloseq::sample_names(phy)))
    expect_identical(counts$taxa_id, phyloseq::taxa_names(phy))
  }
})

test_that("otu_table has one canonical orientation", {
  taxa_rows <- otu_table(recipe(make_accessor_phy(TRUE)))
  taxa_columns <- otu_table(recipe(make_accessor_phy(FALSE)))

  expect_identical(taxa_columns, taxa_rows)
  expect_equal(unname(as.matrix(taxa_rows[-1])), matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    nrow = 3,
    byrow = TRUE
  ))
})

test_that("missing metadata and taxonomy produce ID-only tibbles", {
  phy <- make_accessor_phy()
  counts <- phyloseq::otu_table(phy)

  without_metadata <- recipe(phyloseq::phyloseq(
    counts,
    phyloseq::tax_table(phy)
  ))
  without_taxonomy <- recipe(phyloseq::phyloseq(
    counts,
    phyloseq::sample_data(phy)
  ))
  expect_identical(
    sample_data(without_metadata),
    tibble::tibble(sample_id = phyloseq::sample_names(phy))
  )
  expect_identical(
    tax_table(without_taxonomy),
    tibble::tibble(taxa_id = phyloseq::taxa_names(phy))
  )
})

test_that("internal analysis views retain their reduced contracts", {
  rec <- recipe(make_accessor_phy()) |>
    add_model(
      ~ condition + batch,
      targets = c("condition", "batch"),
      tax_level = "Species"
    )

  metadata <- dar:::analysis_sample_data(rec)
  taxonomy <- dar:::analysis_tax_table(rec)

  expect_named(metadata, c("sample_id", "condition", "batch"))
  expect_named(taxonomy, c("taxa_id", "taxa"))
  expect_identical(taxonomy$taxa, tax_table(rec)$Species)
})

test_that("model cohort preparation preserves complete metadata types", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition",
      tax_level = "Species",
      time = "time"
    )

  modeled <- dar:::apply_model_to_recipe(rec)$rec
  metadata <- sample_data(modeled)

  expect_named(
    metadata,
    c("sample_id", "condition", "time", "subject", "batch", "age")
  )
  expect_true(is.numeric(metadata$time))
  expect_true(is.numeric(metadata$age))
})
