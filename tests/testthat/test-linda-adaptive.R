linda_adaptive_fixture <- function(depth_associated = TRUE, zeros = TRUE) {
  counts <- outer(seq_len(8), seq_len(20), function(taxon, sample) {
    120 + 30 * taxon + ((taxon * sample * 11) %% 37) + 8 * (sample %% 4)
  })
  if (zeros) {
    counts[1, c(2, 4, 12, 16)] <- 0
    counts[2, c(3, 8, 14, 20)] <- 0
    counts[3, c(5, 17)] <- 0
  }
  if (depth_associated) {
    counts <- t(t(counts) * rep(c(1, 8), each = 10))
  }
  dimnames(counts) <- list(
    paste0("taxon_", seq_len(nrow(counts))),
    paste0("sample_", seq_len(ncol(counts)))
  )
  metadata <- data.frame(
    condition = factor(rep(c("control", "treated"), each = 10)),
    row.names = colnames(counts)
  )
  list(counts = counts, metadata = metadata)
}

fit_linda_fixture <- function(fun, fixture, ...) {
  fun(
    feature.dat = fixture$counts,
    meta.dat = fixture$metadata,
    formula = "~ condition",
    verbose = FALSE,
    ...
  )
}

linda_fixture_depth_pvalue <- function(fixture) {
  stats::coef(summary(stats::lm(
    log(colSums(fixture$counts)) ~ condition, fixture$metadata
  )))["conditiontreated", "Pr(>|t|)"]
}

test_that("LinDA repairs only a verified local copy of the native function", {
  skip_if_not_installed("MicrobiomeStat")
  skip_if_not(identical(as.character(utils::packageVersion("MicrobiomeStat")), "1.4"))
  native <- getExportedValue("MicrobiomeStat", "linda")
  original_body <- body(native)
  corrected <- dar:::linda_engine_function()

  expect_identical(formals(corrected), formals(native))
  expect_identical(environment(corrected), environment(native))
  expect_identical(body(getExportedValue("MicrobiomeStat", "linda")), original_body)
  expect_false(identical(body(corrected), original_body))
  expect_identical(dar:::engine_stage_function("linda", "fit"), corrected)
  native_lines <- deparse(original_body, width.cutoff = 500L)
  expect_equal(sum(grepl(
    'zero.handling <- "Imputation"', native_lines, fixed = TRUE
  )), 1L)
  expected_body <- parse(text = sub(
    'zero.handling <- "Imputation"', 'zero.handling <- "imputation"',
    native_lines, fixed = TRUE
  ))[[1L]]
  expect_identical(body(corrected), expected_body)
})

test_that("LinDA's local correction preserves NULL and missing call arguments", {
  native <- function(x) {
    zero.handling <- "Imputation"
    covariance <- NULL
    passed <- zero.handling == "imputation"
    list(covariance = covariance, passed = passed, data = x[, 1, drop = FALSE])
  }
  corrected <- dar:::linda_engine_function(native, version = "1.4")
  expected_body <- quote({
    zero.handling <- "imputation"
    covariance <- NULL
    passed <- zero.handling == "imputation"
    list(covariance = covariance, passed = passed, data = x[, 1, drop = FALSE])
  })
  data <- matrix(seq_len(4), nrow = 2)

  expect_identical(body(corrected), expected_body)
  expect_identical(corrected(data), list(
    covariance = NULL, passed = TRUE, data = data[, 1, drop = FALSE]
  ))
})

test_that("LinDA refuses unverified implementations with the capitalization bug", {
  native <- function(zero.handling = NULL) {
    zero.handling <- "Imputation"
    zero.handling == "imputation"
  }
  corrected <- dar:::linda_engine_function(native, version = "1.4")

  expect_identical(corrected(), TRUE)
  expect_identical(native(), FALSE)
  expect_identical(formals(corrected), formals(native))
  expect_identical(environment(corrected), environment(native))
  expect_error(
    dar:::linda_engine_function(native, version = "1.5"),
    class = "dar_error_model_engine"
  )
  expect_error(
    dar:::linda_engine_function(function() {
      zero.handling <- "Imputation"
      zero.handling <- "Imputation"
      zero.handling == "imputation"
    }, version = "1.4"),
    class = "dar_error_model_engine"
  )
  expect_error(
    dar:::linda_engine_function(function() {
      zero.handling <- "Imputation"
      zero.handling
    }, version = "1.4"),
    class = "dar_error_model_engine"
  )
  expect_error(
    dar:::linda_engine_function(function() {
      zero.handling <- "Imputation"
      zero.handling == "imputation"
      zero.handling == "imputation"
    }, version = "1.4"),
    class = "dar_error_model_engine"
  )
})

test_that("LinDA leaves fixed or unrelated native implementations untouched", {
  fixed <- function() {
    zero.handling <- "imputation"
    zero.handling == "imputation"
  }
  unrelated <- function(...) "native result"

  expect_identical(dar:::linda_engine_function(fixed, version = "1.5"), fixed)
  expect_identical(
    dar:::linda_engine_function(unrelated, version = "1.5"), unrelated
  )
})

test_that("adaptive LinDA imputes zeros when depth is associated with design", {
  skip_if_not_installed("MicrobiomeStat")
  skip_if_not(identical(as.character(utils::packageVersion("MicrobiomeStat")), "1.4"))
  fixture <- linda_adaptive_fixture()
  expect_lt(linda_fixture_depth_pvalue(fixture), 0.1)
  corrected <- fit_linda_fixture(
    dar:::linda_engine_function(), fixture, adaptive = TRUE, is.winsor = FALSE
  )
  imputed <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture, adaptive = FALSE,
    zero.handling = "imputation", is.winsor = FALSE
  )
  buggy <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture, adaptive = TRUE, is.winsor = FALSE
  )
  pseudocount <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture, adaptive = FALSE,
    zero.handling = "pseudo-count", is.winsor = FALSE
  )

  expect_equal(corrected, imputed)
  expect_equal(buggy, pseudocount)
  expect_false(isTRUE(all.equal(corrected$feature.dat.use, buggy$feature.dat.use)))
  expect_false(isTRUE(all.equal(corrected$output, buggy$output)))
})

test_that("adaptive LinDA retains the native pseudocount branch", {
  skip_if_not_installed("MicrobiomeStat")
  fixture <- linda_adaptive_fixture(depth_associated = FALSE)
  expect_gt(linda_fixture_depth_pvalue(fixture), 0.1)
  corrected <- fit_linda_fixture(
    dar:::linda_engine_function(), fixture, adaptive = TRUE, is.winsor = FALSE
  )
  native <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture, adaptive = TRUE, is.winsor = FALSE
  )
  pseudocount <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture, adaptive = FALSE,
    zero.handling = "pseudo-count", is.winsor = FALSE
  )

  expect_equal(corrected, native)
  expect_equal(corrected, pseudocount)
})

test_that("LinDA without zeros keeps the native analysis unchanged", {
  skip_if_not_installed("MicrobiomeStat")
  fixture <- linda_adaptive_fixture(zeros = FALSE)
  corrected <- fit_linda_fixture(dar:::linda_engine_function(), fixture)
  native <- fit_linda_fixture(MicrobiomeStat::linda, fixture)

  expect_equal(corrected, native)
})

test_that("adaptive LinDA still chooses zero handling after native preprocessing", {
  skip_if_not_installed("MicrobiomeStat")
  fixture <- linda_adaptive_fixture(depth_associated = FALSE)
  rare <- numeric(ncol(fixture$counts))
  rare[16:20] <- 1e6
  fixture$counts <- rbind(fixture$counts, rare = rare)
  expect_lt(linda_fixture_depth_pvalue(fixture), 0.1)
  corrected <- fit_linda_fixture(
    dar:::linda_engine_function(), fixture,
    prev.filter = 0.3, adaptive = TRUE,
    is.winsor = TRUE, outlier.pct = 0.1
  )
  pseudocount <- fit_linda_fixture(
    MicrobiomeStat::linda, fixture,
    prev.filter = 0.3, adaptive = FALSE, zero.handling = "pseudo-count",
    is.winsor = TRUE, outlier.pct = 0.1
  )

  expect_false("rare" %in% rownames(corrected$feature.dat.use))
  expect_equal(corrected, pseudocount)
  filtered <- fixture$counts[rownames(fixture$counts) != "rare", , drop = FALSE]
  winsorized <- getFromNamespace("winsor.fun", "MicrobiomeStat")(
    filtered, 0.9, "count"
  )
  expect_equal(corrected$feature.dat.use, winsorized + 0.5)
})

test_that("centralized LinDA routes adaptive fitting and native reject decisions", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  compiled <- dar:::compile_model(rec, "linda")
  coefficient <- dar:::reparameterize_model_contrast(
    compiled, compiled$contrasts[1, , drop = FALSE]
  )$coefficient
  native <- data.frame(
    log2FoldChange = c(1, -1), pvalue = c(0.01, 0.02),
    padj = c(0.01, 0.02), reject = c(FALSE, TRUE),
    row.names = c("taxon_1", "taxon_2")
  )
  observed <- list()
  testthat::local_mocked_bindings(
    linda_engine_function = function(...) {
      function(...) {
        observed <<- list(...)
        list(output = stats::setNames(list(native), coefficient))
      }
    },
    .package = "dar"
  )
  prepared <- prep(
    step_linda(rec, id = "linda__adaptive", engine_args = list(
      fit = list(verbose = FALSE)
    )),
    parallel = FALSE
  )

  expect_identical(observed$adaptive, TRUE)
  expect_identical(observed$feature.dat.type, "count")
  expect_identical(observed$zero.handling, "pseudo-count")
  expect_identical(observed$verbose, FALSE)
  expect_identical(tidy_results(prepared)$significant, c(FALSE, TRUE))
  expect_identical(
    tidy_results(prepared, significant_only = TRUE)$taxa_id, "taxon_2"
  )
})
