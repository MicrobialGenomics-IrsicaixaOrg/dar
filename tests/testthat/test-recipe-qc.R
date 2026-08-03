make_qc_phy <- function(taxa_are_rows = TRUE) {
  counts <- rbind(
    taxon_1 = c(0, 1, 2, 0),
    taxon_2 = c(0, 0, 2, 2),
    taxon_3 = c(1, 1, 1, 1)
  )
  colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))

  metadata <- data.frame(
    cohort = factor(c("A", "A", "B", "B"), levels = c("B", "A")),
    time = c(0, 0, 1, 1),
    target_a = factor(c("x", "x", "y", "y")),
    target_b = factor(c("u", "v", "u", "v")),
    literal = factor(
      c("all", "control", NA, "control"),
      levels = c("control", "all")
    ),
    row.names = colnames(counts)
  )
  taxonomy <- matrix(
    paste0("Species_", seq_len(nrow(counts))),
    ncol = 1,
    dimnames = list(rownames(counts), "Species")
  )
  otu <- if (taxa_are_rows) counts else t(counts)

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = taxa_are_rows),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(metadata)
  )
}

qc_columns <- c(
  "group_by", "group", "n", "n_zero", "pct_zero", "pct_all_zero",
  "pct_singletons", "pct_doubletons", "n_samples", "lib_size_min",
  "lib_size_max", "count_mean", "count_max"
)

test_that("recipe_qc returns exact overall metrics without a model", {
  rec <- recipe(make_qc_phy())
  result <- recipe_qc(rec)

  expect_named(result, qc_columns)
  expect_equal(nrow(result), 1L)
  expect_identical(result$group_by, NA_character_)
  expect_identical(result$group, NA_character_)
  expect_equal(
    result[setdiff(qc_columns, c("group_by", "group"))],
    tibble::tibble(
      n = 12, n_zero = 4, pct_zero = 100 / 3,
      pct_all_zero = 0, pct_singletons = 0, pct_doubletons = 200 / 3,
      n_samples = 4, lib_size_min = 1, lib_size_max = 5,
      count_mean = 11 / 12, count_max = 1.5
    )
  )
  expect_true(all(vapply(result[-c(1, 2)], is.numeric, logical(1))))
})

test_that("overall recipe_qc does not require sample metadata", {
  phy <- make_qc_phy()
  without_metadata <- phyloseq::phyloseq(
    phyloseq::otu_table(phy),
    phyloseq::tax_table(phy)
  )
  rec <- recipe(without_metadata)

  expect_equal(recipe_qc(rec)$n, 12)
  expect_error(
    recipe_qc(rec, group_by = "cohort"),
    class = "dar_error_invalid_qc_group"
  )
})

test_that("recipe_qc adds exact grouped metrics in deterministic order", {
  result <- recipe_qc(recipe(make_qc_phy()), group_by = "cohort")

  expect_named(result, qc_columns)
  expect_identical(result$group_by, c(NA_character_, "cohort", "cohort"))
  expect_identical(result$group, c(NA_character_, "B", "A"))
  expect_equal(result$n, c(12, 6, 6))
  expect_equal(result$n_zero, c(4, 1, 3))
  expect_equal(result$pct_zero, c(100 / 3, 100 / 6, 50))
  expect_equal(result$pct_all_zero, c(0, 0, 100 / 3))
  expect_equal(result$pct_singletons, c(0, 100 / 3, 100 / 3))
  expect_equal(result$pct_doubletons, c(200 / 3, 200 / 3, 100 / 3))
  expect_equal(result$n_samples, c(4, 2, 2))
  expect_equal(result$lib_size_min, c(1, 3, 1))
  expect_equal(result$lib_size_max, c(5, 5, 2))
  expect_equal(result$count_mean, c(11 / 12, 4 / 3, 0.5))
  expect_equal(result$count_max, c(1.5, 2, 1))
})

test_that("recipe_qc handles numeric, missing, and literal all groups", {
  rec <- recipe(make_qc_phy())

  numeric <- recipe_qc(rec, group_by = "time")
  expect_identical(numeric$group, c(NA_character_, "0", "1"))

  literal <- recipe_qc(rec, group_by = "literal")
  expect_identical(
    literal$group_by,
    c(NA_character_, "literal", "literal", "literal")
  )
  expect_identical(
    literal$group,
    c(NA_character_, "control", "all", NA_character_)
  )
  expect_equal(literal$n_samples, c(4, 2, 1, 1))
})

test_that("recipe_qc is independent of centralized model targets", {
  phy <- make_qc_phy()
  plain <- recipe(phy)
  modeled <- recipe(phy) |>
    add_model(
      ~ target_a + target_b,
      targets = c("target_a", "target_b"),
      tax_level = "Species"
    )

  expect_equal(recipe_qc(modeled), recipe_qc(plain))
  expect_equal(
    recipe_qc(modeled, group_by = "cohort"),
    recipe_qc(plain, group_by = "cohort")
  )
})

test_that("recipe_qc accepts prepared recipes after preprocessing", {
  prepared <- recipe(make_qc_phy()) |>
    step_filter_taxa(.f = function(x) sum(x) >= 4) |>
    prep(parallel = FALSE)
  result <- recipe_qc(prepared, group_by = "cohort")

  expect_s4_class(prepared, "PrepRecipe")
  expect_equal(result$n, c(8, 4, 4))
  expect_equal(result$n_samples, c(4, 2, 2))
})

test_that("recipe_qc is invariant to the phyloseq count orientation", {
  taxa_rows <- recipe_qc(recipe(make_qc_phy(TRUE)), group_by = "cohort")
  taxa_columns <- recipe_qc(recipe(make_qc_phy(FALSE)), group_by = "cohort")

  expect_equal(taxa_columns, taxa_rows)
})

test_that("recipe_qc extracts counts and grouping metadata once per call", {
  native_qc_count_data <- dar:::qc_count_data
  native_qc_sample_metadata <- dar:::qc_sample_metadata
  count_calls <- 0L
  metadata_calls <- 0L
  testthat::local_mocked_bindings(
    qc_count_data = function(rec) {
      count_calls <<- count_calls + 1L
      native_qc_count_data(rec)
    },
    qc_sample_metadata = function(rec) {
      metadata_calls <<- metadata_calls + 1L
      native_qc_sample_metadata(rec)
    },
    .package = "dar"
  )

  recipe_qc(recipe(make_qc_phy()), group_by = "cohort")
  expect_identical(count_calls, 1L)
  expect_identical(metadata_calls, 1L)
})

test_that("recipe_qc validates recipes and grouping specifications", {
  rec <- recipe(make_qc_phy())

  expect_error(recipe_qc(1), class = "dar_error_invalid_input")
  expect_error(
    recipe_qc(rec, group_by = c("cohort", "time")),
    class = "dar_error_invalid_qc_group"
  )
  expect_error(
    recipe_qc(rec, group_by = NA_character_),
    class = "dar_error_invalid_qc_group"
  )
  expect_error(
    recipe_qc(rec, group_by = "missing"),
    class = "dar_error_invalid_qc_group"
  )

  phy <- make_qc_phy()
  metadata <- as(phyloseq::sample_data(phy), "data.frame")
  metadata$nested <- I(lapply(seq_len(nrow(metadata)), list))
  phyloseq::sample_data(phy) <- phyloseq::sample_data(metadata)
  expect_error(
    recipe_qc(recipe(phy), group_by = "nested"),
    class = "dar_error_invalid_qc_group"
  )

  corrupt <- rec
  corrupt@var_info <- tibble::tibble(wrong = "cohort")
  expect_error(recipe_qc(corrupt), class = "dar_error_invalid_recipe")
})

test_that("phy_qc is defunct and points to recipe_qc", {
  expect_error(
    phy_qc(recipe(make_qc_phy())),
    regexp = "recipe_qc",
    class = "defunctError"
  )
})
