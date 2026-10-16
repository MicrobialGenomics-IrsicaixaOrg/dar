make_abundance_selection_recipe <- function() {
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  for (id in c("deseq_a", "deseq_b", "wilcox_c")) {
    engine <- if (id == "wilcox_c") "wilcox" else "deseq"
    rec <- dar:::add_step(rec, dar:::step(
      engine, id = id, padj_cutoff = if (id == "deseq_b") 0.1 else 0.05
    ))
  }
  plan <- dar:::resolve_model(rec)$contrast_plan
  rows <- tidyr::crossing(
    dar:::analysis_tax_table(rec),
    dplyr::select(plan, "contrast_id", "comparison", "contrast_type", "var")
  ) |>
    dplyr::mutate(effect = 1, padj = 0.9, signif = FALSE)
  first <- rows
  second <- rows
  third <- rows
  selected <- rows$contrast_id == plan$contrast_id[[1]]
  first$signif <- selected & rows$taxa_id %in% c(
    "taxon_1", "taxon_2", "taxon_4", "taxon_5"
  )
  second$signif <- selected & rows$taxa_id %in% c(
    "taxon_1", "taxon_3", "taxon_4", "taxon_5"
  )
  third$signif <- selected & rows$taxa_id == "taxon_6"
  first$effect[first$taxa_id == "taxon_4"] <- -1
  first$effect[first$taxa_id == "taxon_5"] <- NA_real_
  second$effect[second$taxa_id == "taxon_5"] <- NA_real_
  # A different contrast must not contribute to the selected one.
  first$signif[rows$contrast_id == plan$contrast_id[[2]]] <- TRUE
  second$signif[rows$contrast_id == plan$contrast_id[[2]]] <- TRUE
  # Preserve native decisions, rather than recomputing them from padj.
  first$padj[first$signif] <- 0.8
  second$padj[second$signif] <- 0.8
  dar:::prep_recipe(
    rec,
    list(
      deseq_a = list(model = first), deseq_b = list(model = second),
      wilcox_c = list(model = third)
    ),
    list()
  )
}

test_that("abundance selection supports individual configured steps", {
  rec <- make_abundance_selection_recipe()
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  select <- function(steps, selection = "shared", top_n = 20) {
    dar:::select_abundance_taxa(
      rec, steps, selection, "condition", contrast, top_n
    )
  }
  expect_identical(
    select("deseq_a"), c("taxon_1", "taxon_2", "taxon_4", "taxon_5")
  )
  expect_identical(
    select("deseq_b"), c("taxon_1", "taxon_3", "taxon_4", "taxon_5")
  )
  expect_identical(select("wilcox_c"), "taxon_6")
  expect_identical(select(c("deseq_a", "deseq_b")), "taxon_1")
  expect_identical(select(steps_ids(rec, "da")), character())
  expect_setequal(
    select(c("deseq_a", "deseq_b"), "union"), paste0("taxon_", 1:5)
  )
  expect_length(select(c("deseq_a", "deseq_b"), "union", 2), 2L)
})

test_that("abundance plots do not silently replace empty shared selections", {
  rec <- make_abundance_selection_recipe()
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  for (type in c("boxplot", "heatmap")) {
    expect_message(
      plot <- abundance_plt(rec, contrast_id = contrast, type = type),
      "No significant taxa are shared"
    )
    expect_null(plot)
  }
})

test_that("boxplots and heatmaps use the same selected-method discoveries", {
  skip_if_not_installed("microbiome")
  grDevices::pdf(withr::local_tempfile(fileext = ".pdf"))
  withr::defer(grDevices::dev.off())
  rec <- make_abundance_selection_recipe()
  original <- rec
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  args <- list(
    rec = rec, steps = "deseq_a", contrast_id = contrast, transform = NULL
  )
  box <- do.call(abundance_plt, args)
  ids <- c("taxon_1", "taxon_2", "taxon_4", "taxon_5")
  expect_setequal(box$data$taxa_id, ids)
  union <- abundance_plt(
    rec, steps = c("deseq_a", "deseq_b"), selection = "union",
    contrast_id = contrast, transform = NULL
  )
  expect_setequal(union$data$taxa_id, paste0("taxon_", 1:5))
  explicit <- abundance_plt(
    rec, taxa_ids = c("taxon_3", "taxon_6"), steps = "deseq_a",
    top_n = 1, transform = NULL
  )
  expect_setequal(explicit$data$taxa_id, c("taxon_3", "taxon_6"))
  expect_identical(rec, original)
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("gplots")
  withr::defer(ComplexHeatmap::ht_opt(RESET = TRUE))
  heat <- do.call(abundance_plt, c(args, list(type = "heatmap")))
  expect_s4_class(heat, "HeatmapList")
  expected <- box$data |>
    dplyr::select("taxa", "sample_id", "value") |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "value") |>
    data.frame(row.names = 1, check.names = FALSE) |>
    as.matrix()
  expect_identical(heat@ht_list[[1]]@matrix, expected)
})

test_that("abundance selection rejects invalid inputs before plotting", {
  rec <- make_abundance_selection_recipe()
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  for (steps in list(character(), "missing", c("deseq_a", "deseq_a"), NA_character_)) {
    expect_error(
      abundance_plt(rec, steps = steps, contrast_id = contrast),
      class = "dar_error_invalid_steps"
    )
  }
  for (top_n in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(
      abundance_plt(rec, top_n = top_n, contrast_id = contrast),
      class = "dar_error_invalid_top_n"
    )
  }
  for (taxa_ids in list("missing", NA_character_, c("taxon_1", "taxon_1"))) {
    expect_error(
      abundance_plt(rec, taxa_ids = taxa_ids),
      class = "dar_error_invalid_taxa_ids"
    )
  }
  expect_error(abundance_plt(rec), class = "dar_error_ambiguous_plot_contrast")
  expect_error(abundance_plt(rec, contrast_id = contrast, selection = "typo"))
  expect_message(
    result <- abundance_plt(rec, taxa_ids = character()), "No taxa were selected"
  )
  expect_null(result)
})
