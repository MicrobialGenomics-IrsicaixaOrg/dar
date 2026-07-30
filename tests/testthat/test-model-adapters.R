model_smoke_recipe <- function() {
  recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
}

expect_complete_model_result <- function(result, rec) {
  table <- dar:::flatten_model_output(result)
  expect_setequal(
    unique(table$contrast_id),
    dar:::resolve_model(rec)$contrast_plan$contrast_id
  )
  expect_true(all(c(
    "taxa_id", "effect", "padj", "signif", "contrast_id", "comparison",
    "contrast_type", "var"
  ) %in% names(table)))
}

test_that("DESeq2 compiles and executes a centralized contrast", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  rec <- model_smoke_recipe()
  result <- dar:::run_deseq_model(
    rec, "Wald", "mean", FALSE, "normal", 0.05, 0, FALSE
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("ALDEx2 compiles and executes a centralized contrast", {
  skip_if_not_installed("ALDEx2")
  withr::local_seed(130)
  rec <- model_smoke_recipe()
  result <- suppressWarnings(
    dar:::run_aldex_model(rec, 0.05, 8, "all", FALSE)
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("ANCOM-BC2 compiles and executes a centralized contrast", {
  skip_if_not_installed("ANCOMBC")
  rec <- model_smoke_recipe()
  result <- suppressWarnings(
    dar:::run_ancom_model(
      rec, "holm", 0.1, 0, 0.05, FALSE, FALSE, 0.05, 1, FALSE, FALSE
    )
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("corncob compiles and executes a centralized contrast", {
  skip_if_not_installed("corncob")
  rec <- model_smoke_recipe()
  result <- dar:::run_corncob_model(
    rec, ~ 1, "logit", "logit", TRUE, 0.05, "fdr", 0, FALSE
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("MaAsLin3 compiles and executes a centralized contrast", {
  skip_if_not_installed("maaslin3")
  rec <- model_smoke_recipe()
  result <- dar:::run_maaslin_model(
    rec, 0, 0.1, 0, "TSS", "LOG", 0.1, "BH", TRUE, TRUE, FALSE
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("Wilcoxon executes only its supported simple model", {
  skip_if_not_installed("microbiome")
  rec <- model_smoke_recipe()
  result <- dar:::run_wilcox_model(
    rec, "compositional", 0.05, "BH", FALSE
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("LEfSe executes only its supported simple model", {
  skip_if_not_installed("lefser")
  withr::local_seed(130)
  rec <- model_smoke_recipe()
  result <- suppressWarnings(
    dar:::run_lefse_model(rec, 0.05, 0.05, 2, 1L, FALSE, FALSE)
  )
  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("LEfSe executes every pairwise contrast for multilevel targets", {
  skip_if_not_installed("lefser")
  withr::local_seed(130)
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  result <- suppressWarnings(
    dar:::run_lefse_model(rec, 0.05, 0.05, 2, 1L, FALSE, FALSE)
  )

  expect_complete_model_result(result, rec)
  table <- dar:::flatten_model_output(result)
  expect_length(unique(table$contrast_id), 3L)
  expect_identical(
    anyDuplicated(dplyr::select(table, "taxa_id", "contrast_id")),
    0L
  )
  expect_true(all(is.finite(table$effect)))
})
