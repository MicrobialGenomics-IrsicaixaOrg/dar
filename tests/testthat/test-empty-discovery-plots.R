make_empty_discovery_recipe <- function(discoveries = FALSE) {
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step("deseq", id = "deseq__one"))
  rec <- dar:::add_step(rec, dar:::step("maaslin", id = "maaslin__two"))
  plan <- dar:::resolve_model(rec)$contrast_plan
  rows <- tidyr::crossing(
    dar:::analysis_tax_table(rec),
    dplyr::select(plan, "contrast_id", "comparison", "contrast_type", "var")
  ) |>
    dplyr::mutate(
      effect = 1,
      signif = discoveries & .data$contrast_id == plan$contrast_id[[1]],
      padj = ifelse(.data$signif, 0.01, 0.99)
    )
  dar:::prep_recipe(
    rec,
    list(deseq__one = list(model = rows), maaslin__two = list(model = rows)),
    list()
  )
}

test_that("valid selectors do not depend on significant discoveries", {
  for (discoveries in c(FALSE, TRUE)) {
    rec <- make_empty_discovery_recipe(discoveries)
    contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[2]]
    expect_no_error(
      selected <- dar:::selected_significant(
        rec, target = "condition", contrast_id = contrast
      )
    )
    expect_equal(nrow(selected), 0L)
  }
})

test_that("empty discovery plots report no data instead of invalid selectors", {
  rec <- make_empty_discovery_recipe()
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[2]]
  for (type in c("boxplot", "heatmap")) {
    expect_message(
      plot <- abundance_plt(
        rec, target = "condition", contrast_id = contrast, type = type
      ),
      "No significant taxa"
    )
    expect_null(plot)
  }
  expect_message(
    plot <- mutual_plt(
      rec, target = "condition", contrast_id = contrast, count_cutoff = 2
    ),
    "No significant taxa"
  )
  expect_null(plot)
  expect_message(
    plot <- intersection_plt(rec, target = "condition", contrast_id = contrast),
    "No significant taxa"
  )
  expect_null(plot)
})

test_that("a contrast without discoveries remains valid when others have hits", {
  rec <- make_empty_discovery_recipe(discoveries = TRUE)
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[2]]
  expect_message(
    plot <- mutual_plt(rec, count_cutoff = 2, comparisons = contrast),
    "No significant taxa"
  )
  expect_null(plot)
  expect_message(
    plot <- abundance_plt(rec, target = "condition", contrast_id = contrast),
    "No significant taxa"
  )
  expect_null(plot)
  selected <- intersection_df(rec, target = "condition", contrast_id = contrast)
  expect_equal(nrow(selected), phyloseq::ntaxa(get_phy(rec)))
  expect_true(all(selected$deseq__one == 0L))
  expect_true(all(selected$maaslin__two == 0L))
  plot <- exclusion_plt(rec, target = "condition", contrast_id = contrast)
  expect_equal(nrow(plot$data), 0L)
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

test_that("invalid selectors are still rejected without discoveries", {
  rec <- make_empty_discovery_recipe()
  for (value in list("missing", NA_character_, character(), c("a", "b"))) {
    expect_error(
      dar:::selected_significant(rec, target = value),
      class = "dar_error_invalid_contrast_selector"
    )
    expect_error(
      dar:::selected_significant(rec, contrast_id = value),
      class = "dar_error_invalid_contrast_selector"
    )
  }
  expect_error(
    abundance_plt(rec, target = "condition", contrast_id = "missing"),
    class = "dar_error_invalid_contrast_selector"
  )
  expect_error(
    mutual_plt(rec, count_cutoff = 2, comparisons = "missing"),
    class = "dar_error_invalid_comparison"
  )
  expect_error(
    intersection_plt(rec, target = "missing"),
    class = "dar_error_invalid_contrast_selector"
  )
})

test_that("targets and contrasts must identify the same planned hypotheses", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + batch, targets = c("condition", "batch"),
      tax_level = "Species"
    )
  rec <- dar:::add_step(rec, dar:::step("deseq", id = "deseq__one"))
  plan <- dar:::resolve_model(rec)$contrast_plan
  rows <- tidyr::crossing(
    dar:::analysis_tax_table(rec),
    dplyr::select(plan, "contrast_id", "comparison", "contrast_type", "var")
  ) |>
    dplyr::mutate(effect = 1, padj = 0.99, signif = FALSE)
  rec <- dar:::prep_recipe(rec, list(deseq__one = list(model = rows)), list())
  contrast <- plan$contrast_id[plan$var == "batch"][[1]]
  expect_error(
    dar:::selected_significant(rec, target = "condition", contrast_id = contrast),
    class = "dar_error_invalid_contrast_selector"
  )
  expect_no_error(
    selected <- dar:::selected_significant(rec, target = "batch", contrast_id = contrast)
  )
  expect_equal(nrow(selected), 0L)
})

test_that("nonempty discoveries and explicit abundance selections still plot", {
  skip_if_not_installed("microbiome")
  plot_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  withr::defer(grDevices::dev.off())
  rec <- make_empty_discovery_recipe(discoveries = TRUE)
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  plot <- mutual_plt(rec, count_cutoff = 2, contrast_id = contrast)
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  expect_identical(unique(plot$data$contrast_id), contrast)
  expect_true(all(plot$data$method_count == 2L))

  empty <- make_empty_discovery_recipe()
  original <- get_phy(empty)
  plot <- abundance_plt(empty, taxa_ids = "taxon_1", transform = NULL)
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), phyloseq::nsamples(original))
  expect_identical(get_phy(empty), original)
})
