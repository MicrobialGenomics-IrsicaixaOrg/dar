make_effect_plot_recipe <- function(multilevel = FALSE, legacy = FALSE,
                                     duplicated_labels = FALSE) {
  phy <- if (multilevel) make_multilevel_phy() else make_longitudinal_phy()
  if (duplicated_labels) {
    taxonomy <- as(phyloseq::tax_table(phy), "matrix")
    taxonomy[, "Species"] <- c("Repeated", "Repeated", paste0("Species_", 3:6))
    taxonomy <- cbind(
      Genus = paste0("Genus_", seq_len(nrow(taxonomy))), taxonomy
    )
    phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  }
  rec <- if (legacy) {
    suppressWarnings(recipe(phy, var_info = "condition", tax_info = "Species"))
  } else {
    recipe(phy)
  }
  if (!legacy) {
    rec <- add_model(
      rec, ~ condition, targets = "condition", tax_level = "Species"
    )
  }
  rec <- dar:::add_step(rec, dar:::step("deseq", id = "deseq__effects"))
  rec <- dar:::add_step(rec, dar:::step("wilcox", id = "wilcox__effects"))
  plan <- if (legacy) {
    tibble::tibble(
      comparison = "treated-control", var = "condition"
    )
  } else {
    dar:::resolve_model(rec)$contrast_plan |>
      dplyr::select("contrast_id", "comparison", "contrast_type", "var")
  }
  rows <- tidyr::crossing(dar:::analysis_tax_table(rec), plan) |>
    dplyr::mutate(
      effect = rep(c(-4, 3, 2, 1, 0, NA_real_), each = nrow(plan)),
      signif = rep(c(TRUE, FALSE, NA, FALSE, FALSE, TRUE), each = nrow(plan)),
      padj = rep(c(0.9, 0.001, NA_real_, 0.5, 0.9, 0.01), each = nrow(plan))
    )
  wilcox_rows <- dplyr::mutate(
    rows, effect = rep(c(0.01, -0.02, 0.03, -0.04, 0, NA_real_), each = nrow(plan))
  )
  dar:::prep_recipe(
    rec,
    list(
      deseq__effects = list(model = rows),
      wilcox__effects = list(model = wilcox_rows)
    ),
    list()
  )
}

test_that("effect plots preserve native estimates, decisions and separate scales", {
  rec <- make_effect_plot_recipe()
  before <- rec
  expect_no_message(plot <- effect_plt(rec))
  expect_s3_class(plot, "ggplot")
  expect_identical(rec, before)
  expect_equal(nrow(plot$data), 10L)
  expect_setequal(
    unique(plot$data$effect_metric),
    c("log2_fold_change", "median_compositional_difference")
  )
  expect_identical(plot$facet$params$free, list(x = TRUE, y = TRUE))
  expect_length(unique(plot$data$panel), 2L)
  expect_false(any(grepl("ymin|ymax|xmin|xmax", names(plot$mapping))))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  significant <- subset(plot$data, step_id == "deseq__effects" & taxa_id == "taxon_1")
  nonsignificant <- subset(plot$data, step_id == "deseq__effects" & taxa_id == "taxon_2")
  unavailable <- subset(plot$data, step_id == "deseq__effects" & taxa_id == "taxon_3")
  expect_true(significant$significant)
  expect_identical(as.character(significant$decision), "Significant")
  expect_gt(significant$adj_p_value, 0.05)
  expect_false(nonsignificant$significant)
  expect_identical(as.character(nonsignificant$decision), "Not significant")
  expect_lt(nonsignificant$adj_p_value, 0.05)
  expect_true(is.na(unavailable$significant))
  expect_identical(as.character(unavailable$decision), "Not available")
})

test_that("automatic effect plot selection ranks separately within each panel", {
  rec <- make_effect_plot_recipe()
  plot <- suppressMessages(effect_plt(rec, top_n = 2))
  expect_equal(nrow(plot$data), 4L)
  expect_identical(
    plot$data$taxa_id[plot$data$step_id == "deseq__effects"],
    c("taxon_1", "taxon_2")
  )
  expect_identical(
    plot$data$taxa_id[plot$data$step_id == "wilcox__effects"],
    c("taxon_4", "taxon_3")
  )
  expect_identical(as.integer(table(plot$data$panel)), c(2L, 2L))
  reversed <- suppressMessages(effect_plt(
    rec, steps = c("wilcox__effects", "deseq__effects"), top_n = 1
  ))
  expect_identical(unique(reversed$data$step_id), c("wilcox__effects", "deseq__effects"))
  tied <- rec
  tied@results$deseq__effects$model$effect[1:2] <- 4
  plot <- suppressMessages(effect_plt(tied, steps = "deseq__effects", top_n = 1))
  expect_identical(plot$data$taxa_id, "taxon_1")
})

test_that("explicit taxa and step selection do not require shared discoveries", {
  rec <- make_effect_plot_recipe()
  plot <- effect_plt(
    rec, taxa_ids = c("taxon_2", "taxon_4", "taxon_5"),
    steps = "wilcox__effects", top_n = 1
  )
  expect_equal(nrow(plot$data), 3L)
  expect_identical(unique(plot$data$step_id), "wilcox__effects")
  expect_true(all(!plot$data$significant))
  expect_setequal(plot$data$taxa_id, c("taxon_2", "taxon_4", "taxon_5"))
  expect_true(any(plot$data$effect_size == 0))
  absent <- rec
  absent@results$wilcox__effects$model <- subset(
    absent@results$wilcox__effects$model, taxa_id != "taxon_2"
  )
  plot <- effect_plt(
    absent, taxa_ids = c("taxon_2", "taxon_4"), steps = "wilcox__effects"
  )
  expect_identical(plot$data$taxa_id, "taxon_4")
  expect_message(
    empty <- effect_plt(absent, taxa_ids = "taxon_2", steps = "wilcox__effects"),
    "No finite effects"
  )
  expect_null(empty)
  expect_message(empty <- effect_plt(rec, taxa_ids = character()), "No finite effects")
  expect_null(empty)
})

test_that("effect plot omission never converts nonfinite estimates to zero", {
  rec <- make_effect_plot_recipe()
  rec@results$deseq__effects$model$effect <- c(NA_real_, Inf, -Inf, NA, Inf, NA)
  expect_message(
    result <- effect_plt(rec, steps = "deseq__effects"),
    "No finite effects"
  )
  expect_null(result)
  expect_message(
    effect_plt(rec, taxa_ids = "taxon_6", steps = "wilcox__effects"),
    "No finite effects"
  )
  empty <- dar:::prep_recipe(
    recipe(make_longitudinal_phy()) |>
      add_model(~ condition, targets = "condition", tax_level = "Species"),
    list(), list()
  )
  expect_message(result <- effect_plt(empty), "No finite effects")
  expect_null(result)
})

test_that("effect plots require explicit selection of ambiguous contrasts", {
  rec <- make_effect_plot_recipe(multilevel = TRUE)
  plan <- dar:::resolve_model(rec)$contrast_plan
  expect_error(effect_plt(rec), class = "dar_error_ambiguous_plot_contrast")
  expect_error(
    effect_plt(rec, target = "condition"),
    class = "dar_error_ambiguous_plot_contrast"
  )
  selected <- plan$contrast_id[[2]]
  plot <- suppressMessages(effect_plt(rec, contrast_id = selected))
  expect_identical(unique(plot$data$contrast_id), selected)
  for (selector in list("missing", NA_character_, character(), c("one", "two"))) {
    expect_error(
      effect_plt(rec, contrast_id = selector),
      class = "dar_error_invalid_contrast_selector"
    )
    expect_error(
      effect_plt(rec, target = selector),
      class = "dar_error_invalid_contrast_selector"
    )
  }
})

test_that("effect plots select contrast IDs across multiple targets", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + batch, targets = c("condition", "batch"),
      tax_level = "Species"
    )
  rec <- dar:::add_step(rec, dar:::step("deseq", id = "deseq__effects"))
  plan <- dar:::resolve_model(rec)$contrast_plan
  rows <- tidyr::crossing(
    dar:::analysis_tax_table(rec),
    dplyr::select(plan, "contrast_id", "comparison", "contrast_type", "var")
  ) |>
    dplyr::mutate(effect = 1, signif = FALSE, padj = 0.9)
  rec <- dar:::prep_recipe(rec, list(deseq__effects = list(model = rows)), list())
  contrast <- plan$contrast_id[plan$var == "batch"][[1]]
  plot <- effect_plt(rec, contrast_id = contrast)
  expect_identical(unique(plot$data$var), "batch")
  plot <- effect_plt(rec, target = "batch")
  expect_identical(unique(plot$data$contrast_id), contrast)
  expect_error(
    effect_plt(rec, target = "condition", contrast_id = contrast),
    class = "dar_error_invalid_contrast_selector"
  )
})

test_that("legacy effect plots use actual result contrasts", {
  rec <- make_effect_plot_recipe(legacy = TRUE)
  plot <- suppressMessages(effect_plt(rec, target = "condition"))
  expect_identical(unique(plot$data$contrast_id), "condition[treated-control]")
  expect_error(
    effect_plt(rec, contrast_id = "missing"),
    class = "dar_error_invalid_contrast_selector"
  )
  expect_error(
    effect_plt(rec, target = "missing"),
    class = "dar_error_invalid_contrast_selector"
  )
})

test_that("taxon display labels do not merge distinct lineages", {
  rec <- make_effect_plot_recipe(duplicated_labels = TRUE)
  plot <- effect_plt(
    rec, taxa_ids = c("taxon_1", "taxon_2"), steps = "deseq__effects"
  )
  expect_equal(nrow(plot$data), 2L)
  expect_equal(length(unique(plot$data$.plot_taxon)), 2L)
  labels <- plot$scales$get_scales("y")$labels
  expect_setequal(unname(labels), c("Repeated [taxon_1]", "Repeated [taxon_2]"))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

test_that("effect plot arguments reject invalid selections", {
  rec <- make_effect_plot_recipe()
  for (value in list(0, -1, NA_real_, Inf, 1.5, c(1, 2), "two")) {
    expect_error(effect_plt(rec, top_n = value))
  }
  for (value in list(NA_character_, "", c("taxon_1", "taxon_1"))) {
    expect_error(effect_plt(rec, taxa_ids = value))
  }
  expect_error(
    effect_plt(rec, taxa_ids = "missing"),
    class = "dar_error_invalid_taxa_ids"
  )
  expect_error(
    effect_plt(rec, steps = "missing"),
    class = "dar_error_invalid_step_id"
  )
  expect_error(effect_plt(recipe(make_longitudinal_phy())))
})
