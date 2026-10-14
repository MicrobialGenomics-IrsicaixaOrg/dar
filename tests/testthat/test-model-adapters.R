model_smoke_recipe <- function() {
  recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
}

record_engine_stages <- function() {
  calls <- character()
  native_stage_function <- dar:::engine_stage_function
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      native <- native_stage_function(engine, stage)
      function(...) {
        calls <<- c(calls, paste(engine, stage, sep = ":"))
        native(...)
      }
    },
    .package = "dar",
    .env = parent.frame()
  )
  function() calls
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
  table <- dar:::flatten_model_output(result)
  expect_false(anyNA(table$taxa_id))
  expect_setequal(table$taxa_id, phyloseq::taxa_names(get_phy(rec)))
})

test_that("ANCOM-BC2 preserves aggregated taxon IDs across planned contrasts", {
  skip_if_not_installed("ANCOMBC")
  phy <- make_multilevel_phy()
  phyloseq::taxa_names(phy) <- paste0("ASV_", seq_len(phyloseq::ntaxa(phy)))
  taxonomy <- cbind(
    Genus = paste0("Genus_", seq_len(phyloseq::ntaxa(phy))),
    as(phyloseq::tax_table(phy), "matrix")
  )
  taxonomy[1:2, "Species"] <- "Shared_species"
  taxonomy[4, c("Genus", "Species")] <- taxonomy[3, c("Genus", "Species")]
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  aggregated <- phyloseq::tax_glom(phy, taxrank = "Species", NArm = FALSE)
  expected_ids <- phyloseq::taxa_names(aggregated)
  expected_taxa <- as(phyloseq::tax_table(aggregated), "matrix")[, "Species"]
  rec <- suppressWarnings(
    recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Species") |>
      step_ancom(prv_cut = 0, id = "ancom__taxonomy")
  )

  prepared <- suppressWarnings(prep(rec, parallel = FALSE))
  table <- tidy_results(prepared)

  expect_false(anyNA(table$taxa_id))
  expect_setequal(table$taxa_id, expected_ids)
  expect_length(expected_ids, 5L)
  expect_equal(nrow(table), length(expected_ids) * 3L)
  expect_equal(sum(table$taxa == "Shared_species"), 6L)
  expect_identical(
    anyDuplicated(dplyr::select(table, "taxa_id", "contrast_id")),
    0L
  )
  expect_setequal(
    table$contrast_id,
    dar:::resolve_model(rec)$contrast_plan$contrast_id
  )
  expect_equal(table$taxa, unname(expected_taxa[table$taxa_id]))
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

test_that("DESeq2 routes advanced arguments through fit and shrink", {
  skip_if_not_installed("DESeq2")
  stages <- record_engine_stages()
  rec <- model_smoke_recipe()
  result <- dar:::run_deseq_model(
    rec, "Wald", "mean", FALSE, "normal", 0.05, 0, FALSE,
    engine_args = dar:::normalize_engine_args("deseq", list(
      fit = list(minReplicatesForReplace = Inf),
      shrink = list(lfcThreshold = 0)
    ))
  )
  expect_complete_model_result(result, rec)
  expect_setequal(
    unique(stages()), c("deseq:size_factors", "deseq:fit", "deseq:shrink")
  )
})

test_that("DESeq2 selects results for interaction contrasts and rejects shrink", {
  skip_if_not_installed("DESeq2")
  stages <- record_engine_stages()
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time, targets = "condition", tax_level = "Species",
      time = "time"
    )
  result <- dar:::run_deseq_model(
    rec, "Wald", "mean", FALSE, "normal", 0.05, 0, FALSE,
    engine_args = dar:::normalize_engine_args(
      "deseq", list(results = list(pAdjustMethod = "BH"))
    )
  )
  expect_complete_model_result(result, rec)
  expect_true("deseq:results" %in% stages())
  expect_false("deseq:shrink" %in% stages())

  expect_error(
    dar:::run_deseq_model(
      rec, "Wald", "mean", FALSE, "normal", 0.05, 0, FALSE,
      engine_args = dar:::normalize_engine_args(
        "deseq", list(shrink = list(lfcThreshold = 1))
      )
    ),
    class = "dar_error_invalid_engine_args"
  )
})

test_that("ALDEx2 routes advanced arguments through all three stages", {
  skip_if_not_installed("ALDEx2")
  stages <- record_engine_stages()
  withr::local_seed(130)
  rec <- model_smoke_recipe()
  result <- suppressWarnings(dar:::run_aldex_model(
    rec, 0.05, 8, "all", FALSE,
    engine_args = dar:::normalize_engine_args("aldex", list(
      clr = list(useMC = FALSE), fit = list(verbose = FALSE),
      effect = list(include.sample.summary = TRUE)
    ))
  ))
  expect_complete_model_result(result, rec)
  expect_setequal(
    unique(stages()), c("aldex:clr", "aldex:fit", "aldex:effect")
  )
})

test_that("ANCOM-BC2 routes advanced arguments through fit", {
  skip_if_not_installed("ANCOMBC")
  stages <- record_engine_stages()
  rec <- model_smoke_recipe()
  result <- suppressWarnings(dar:::run_ancom_model(
    rec, "holm", 0.1, 0, 0.05, FALSE, FALSE, 0.05, 1, FALSE, FALSE,
    engine_args = dar:::normalize_engine_args("ancom", list(
      fit = list(iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE))
    ))
  ))
  expect_complete_model_result(result, rec)
  expect_identical(unique(stages()), "ancom:fit")
})

test_that("corncob routes advanced arguments through fit", {
  skip_if_not_installed("corncob")
  stages <- record_engine_stages()
  rec <- model_smoke_recipe()
  result <- dar:::run_corncob_model(
    rec, ~ 1, "logit", "logit", TRUE, 0.05, "fdr", 0, FALSE,
    engine_args = dar:::normalize_engine_args(
      "corncob", list(fit = list(robust = FALSE))
    )
  )
  expect_complete_model_result(result, rec)
  expect_identical(unique(stages()), "corncob:fit")
})

test_that("MaAsLin3 routes advanced arguments through fit and contrast", {
  skip_if_not_installed("maaslin3")
  stages <- record_engine_stages()
  rec <- model_smoke_recipe()
  result <- dar:::run_maaslin_model(
    rec, 0, 0.1, 0, "TSS", "LOG", 0.1, "BH", TRUE, TRUE, FALSE,
    engine_args = dar:::normalize_engine_args("maaslin", list(
      fit = list(cores = 1), contrast = list(subtract_median = FALSE)
    ))
  )
  expect_complete_model_result(result, rec)
  expect_setequal(unique(stages()), c("maaslin:fit", "maaslin:contrast"))
})

test_that("Wilcoxon routes advanced arguments through test", {
  skip_if_not_installed("microbiome")
  stages <- record_engine_stages()
  rec <- model_smoke_recipe()
  result <- dar:::run_wilcox_model(
    rec, "compositional", 0.05, "BH", FALSE,
    engine_args = dar:::normalize_engine_args(
      "wilcox", list(test = list(paired = FALSE))
    )
  )
  expect_complete_model_result(result, rec)
  expect_identical(unique(stages()), "wilcox:test")
})

test_that("LEfSe routes advanced arguments through transform and fit", {
  skip_if_not_installed("lefser")
  stages <- record_engine_stages()
  withr::local_seed(130)
  rec <- model_smoke_recipe()
  result <- suppressWarnings(dar:::run_lefse_model(
    rec, 0.05, 0.05, 2, 1L, FALSE, FALSE,
    engine_args = dar:::normalize_engine_args(
      "lefse", list(fit = list(checkAbundances = FALSE))
    )
  ))
  expect_complete_model_result(result, rec)
  expect_setequal(unique(stages()), c("lefse:transform", "lefse:fit"))
})
