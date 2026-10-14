corncob_threshold_fixture <- function(effects, padj = rep(0.01, length(effects))) {
  list(
    significant_models = lapply(effects, function(effect) {
      list(
        np.mu = 2L,
        coefficients = rbind(
          intercept = c(Estimate = 0, Std.Error = 0.1),
          conditiontreated = c(Estimate = effect, Std.Error = 0.1)
        )
      )
    }),
    p_fdr = padj
  )
}

test_that("legacy corncob applies the requested absolute effect threshold", {
  effects <- c(-2.1, -2, -0.2, 0, 0.2, 2, 2.1)
  ids <- paste0("taxon_", seq_along(effects))
  result <- dar:::corncob_stats_tbl(
    corncob_threshold_fixture(effects), paste0("Species_", seq_along(effects), " (", ids, ")"),
    c("control", "treated"), "condition", log2FC = 2, fdr_cutoff = 0.05
  )
  expect_identical(result$signif, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(result$effect, effects)
  expect_equal(result$log2FC, effects)
  expect_identical(result$taxa_id, ids)
  expect_true(all(result$comparison == "control_treated"))
})

test_that("legacy corncob still requires the adjusted p-value cutoff", {
  effects <- c(-3, -3, -3, 3, 3, 3)
  padj <- rep(c(0.01, 0.05, 0.2), 2)
  result <- dar:::corncob_stats_tbl(
    corncob_threshold_fixture(effects, padj),
    paste0("Species_", seq_along(effects), " (taxon_", seq_along(effects), ")"),
    c("control", "treated"), "condition", log2FC = 2, fdr_cutoff = 0.05
  )
  expect_identical(result$signif, c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(result$padj, padj)
})

test_that("legacy corncob threshold changes only significance flags", {
  effects <- c(-0.2, 0, 0.2)
  native <- corncob_threshold_fixture(effects)
  taxa <- paste0("Species_", seq_along(effects), " (taxon_", seq_along(effects), ")")
  default <- dar:::corncob_stats_tbl(
    native, taxa, c("control", "treated"), "condition", 0, 0.05
  )
  strict <- dar:::corncob_stats_tbl(
    native, taxa, c("control", "treated"), "condition", 2, 0.05
  )
  expect_true(all(default$signif))
  expect_false(any(strict$signif))
  expect_identical(dplyr::select(default, -signif), dplyr::select(strict, -signif))
})

test_that("legacy prep honors the effect threshold with native corncob fits", {
  skip_if_not_installed("corncob")
  withr::local_seed(130)
  phy <- make_longitudinal_phy()
  taxonomy <- cbind(
    Genus = paste0("Genus_", seq_len(phyloseq::ntaxa(phy))),
    as(phyloseq::tax_table(phy), "matrix")
  )
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  rec <- suppressWarnings(recipe(phy, "condition", "Species")) |>
    step_corncob(
      log2FC = 1e6, fdr_cutoff = 1, filter_discriminant = FALSE,
      id = "corncob__legacy_threshold"
    )
  prepared <- suppressWarnings(prep(rec, parallel = FALSE))
  result <- tidy_results(prepared)
  expect_null(get_model(prepared))
  expect_gt(nrow(result), 0L)
  expect_true(all(is.finite(result$effect_size)))
  expect_true(any(result$adj_p_value < 1))
  expect_false(any(result$significant))
  expect_true(all(result$significant == (
    result$adj_p_value < 1 & abs(result$effect_size) >= 1e6
  )))
  expect_identical(get_phy(prepared), phy)
})
