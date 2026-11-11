capture_deseq_stages <- function() {
  calls <- list()
  native_stage_function <- dar:::engine_stage_function
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      native <- native_stage_function(engine, stage)
      function(...) {
        arguments <- list(...)
        result <- native(...)
        calls[[length(calls) + 1L]] <<- list(
          stage = stage, arguments = arguments, result = result
        )
        result
      }
    },
    .package = "dar",
    .env = parent.frame()
  )
  function(stage) Filter(function(call) call$stage == stage, calls)
}

make_sparse_deseq_phy <- function() {
  phy <- make_multilevel_phy()
  counts <- as(phyloseq::otu_table(phy), "matrix")
  counts[cbind(seq_len(nrow(counts)), seq_len(nrow(counts)))] <- 0
  counts[nrow(counts), ] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  phy
}

test_that("DESeq2 uses native poscounts and passes aligned results to ashr", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("ashr")
  rec <- recipe(make_sparse_deseq_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_deseq(fitType = "mean", id = "deseq__poscounts")
  calls <- capture_deseq_stages()

  prepared <- prep(rec, parallel = FALSE)
  table <- tidy_results(prepared)
  contrasts <- get_execution(prepared)$contrasts
  expect_length(calls("results"), nrow(contrasts))
  expect_true(all(rowSums(as(phyloseq::otu_table(get_phy(rec)), "matrix") == 0) > 0))

  for (index in seq_len(nrow(contrasts))) {
    size_factors <- calls("size_factors")[[index]]
    extracted <- calls("results")[[index]]
    shrunk <- calls("shrink")[[index]]
    expect_identical(size_factors$arguments$type, "poscounts")
    expect_false("geoMeans" %in% names(size_factors$arguments))
    native_dds <- DESeq2::estimateSizeFactors(
      size_factors$arguments$object, type = "poscounts"
    )
    expect_equal(
      DESeq2::sizeFactors(size_factors$result), DESeq2::sizeFactors(native_dds)
    )
    expect_identical(extracted$arguments$alpha, 0.05)
    native <- DESeq2::results(
      extracted$arguments$object, name = extracted$arguments$name, alpha = 0.05
    )
    expect_equal(extracted$result, native)
    expect_identical(shrunk$arguments$res, native)
    expect_equal(shrunk$result$pvalue, native$pvalue)
    expect_equal(shrunk$result$padj, native$padj)
    expect_true(anyNA(native$padj))

    rows <- table[table$contrast_id == contrasts$contrast_id[[index]], ]
    rows <- rows[match(rownames(native), rows$taxa_id), ]
    expect_equal(rows$effect_size, shrunk$result$log2FoldChange)
    expect_equal(rows$p_value, native$pvalue)
    expect_equal(rows$adj_p_value, native$padj)
    expect_identical(rows$significant, !is.na(native$padj) & native$padj < 0.05)
  }
})

test_that("DESeq2 aligns nondefault alpha and results controls before shrinkage", {
  skip_if_not_installed("DESeq2")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_deseq(
      type = "normal", fitType = "mean", max_significance = 0.2,
      engine_args = list(results = list(
        independentFiltering = FALSE, cooksCutoff = FALSE, pAdjustMethod = "BY"
      ))
    )
  calls <- capture_deseq_stages()
  prepared <- prep(rec, parallel = FALSE)
  table <- tidy_results(prepared)
  expect_length(calls("results"), nrow(get_execution(prepared)$contrasts))

  for (index in seq_along(calls("results"))) {
    extracted <- calls("results")[[index]]
    shrunk <- calls("shrink")[[index]]
    native <- DESeq2::results(
      extracted$arguments$object, name = extracted$arguments$name, alpha = 0.2,
      independentFiltering = FALSE, cooksCutoff = FALSE, pAdjustMethod = "BY"
    )
    expect_identical(extracted$arguments$alpha, 0.2)
    expect_equal(extracted$result, native)
    expect_identical(shrunk$arguments$res, native)
    expect_equal(shrunk$result$pvalue, native$pvalue)
    expect_equal(shrunk$result$padj, native$padj)
  }
  expect_identical(
    table$significant, !is.na(table$adj_p_value) & table$adj_p_value < 0.2
  )
})

test_that("DESeq2 checks only the package for the configured shrinkage", {
  checked <- list()
  local_mocked_bindings(
    recipes_pkg_check = function(pkg, ...) {
      checked[[length(checked) + 1L]] <<- pkg
    },
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy())
  expected <- list(
    ashr = c("bioc::DESeq2", "ashr"),
    apeglm = c("bioc::DESeq2", "bioc::apeglm"),
    normal = "bioc::DESeq2"
  )
  for (type in names(expected)) {
    dar:::required_deps(step_deseq(rec, type = type))
  }
  expect_identical(checked, rep(unname(expected), each = 2L))
})

test_that("DESeq2 reserves normalization and alpha for the public contract", {
  for (arguments in list(
    list(size_factors = list(type = "ratio")),
    list(size_factors = list(geoMeans = 1)),
    list(results = list(alpha = 0.1))
  )) {
    expect_error(
      dar:::normalize_engine_args("deseq", arguments),
      class = "dar_error_invalid_engine_args"
    )
  }
})
