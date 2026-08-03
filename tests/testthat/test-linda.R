test_that("step_linda requires a centralized model", {
  skip_if_not_installed("MicrobiomeStat")
  legacy <- recipe(make_longitudinal_phy())

  expect_error(
    step_linda(legacy),
    class = "dar_error_model_required"
  )
})

test_that("step_linda stores validated public and advanced arguments", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  configured <- step_linda(
    rec,
    prev_filter = 0.1,
    winsorize = FALSE,
    zero_handling = "imputation",
    n_cpus = 2L,
    id = "linda__configured",
    engine_args = list(fit = list(verbose = TRUE))
  )
  step <- configured@steps[[1L]]

  expect_s3_class(step, "step_linda")
  expect_identical(step$prev_filter, 0.1)
  expect_false(step$winsorize)
  expect_identical(step$zero_handling, "imputation")
  expect_identical(step$n_cpus, 2L)
  expect_true(step$engine_args$fit$verbose)
  expect_identical(tail(names(formals(step_linda)), 1L), "engine_args")
  expect_error(
    step_linda(dar:::prep_recipe(rec, list(), list())),
    class = "dar_error_invalid_input"
  )
  expect_error(step_linda(rec, prev_filter = 2), class = "simpleError")
  expect_error(step_linda(rec, pseudo_count = 0))
  expect_identical(dar:::required_pkgs_linda(), "MicrobiomeStat")
})

test_that("LinDA configuration survives step export and import", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_linda(
      winsorize = FALSE,
      engine_args = list(fit = list(verbose = TRUE)),
      id = "linda__roundtrip"
    )
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- import_steps(recipe(make_longitudinal_phy()), path)

  expect_identical(imported@steps[[1]], rec@steps[[1]])
})

mock_linda_fit <- function(call_log = NULL, omit = FALSE) {
  force(call_log)
  function(feature.dat, meta.dat, formula, ...) {
    if (!is.null(call_log)) {
      call_log$items[[length(call_log$items) + 1L]] <- list(
        feature.dat = feature.dat, meta.dat = meta.dat,
        formula = formula, arguments = list(...)
      )
    }
    design <- stats::model.matrix(
      reformulas::nobars(stats::as.formula(formula)),
      data = meta.dat
    )
    variables <- colnames(design)[-1L]
    if (omit) {
      variables <- setdiff(variables, variables[[1L]])
    }
    output <- stats::setNames(lapply(seq_along(variables), function(index) {
      data.frame(
        baseMean = rowMeans(feature.dat),
        log2FoldChange = rep(index, nrow(feature.dat)),
        lfcSE = 0.2,
        stat = 1,
        pvalue = 0.01,
        padj = 0.02,
        reject = TRUE,
        df = nrow(meta.dat) - ncol(design),
        row.names = rownames(feature.dat)
      )
    }), variables)
    list(variables = variables, output = output)
  }
}

test_that("LinDA receives the full model and executes every temporal contrast", {
  skip_if_not_installed("MicrobiomeStat")
  calls <- new.env(parent = emptyenv())
  calls$items <- list()
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      expect_identical(engine, "linda")
      expect_identical(stage, "fit")
      mock_linda_fit(calls)
    },
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time",
      contrasts = "all"
    )
  result <- dar:::run_linda_model(
    rec, 0, 0, 0, FALSE, 0.03, FALSE, "pseudo-count", 0.5, 0.1,
    "BH", 0.05, 1L, FALSE,
    engine_args = list(fit = list(verbose = TRUE))
  )
  table <- dar:::flatten_model_output(result)
  plan <- dar:::resolve_model(rec)$contrast_plan

  expect_setequal(unique(table$contrast_id), plan$contrast_id)
  expect_identical(anyDuplicated(table[c("taxa_id", "contrast_id")]), 0L)
  expect_equal(nrow(table), phyloseq::ntaxa(rec@phyloseq) * nrow(plan))
  expect_length(calls$items, nrow(plan))
  purrr::walk(calls$items, function(call) {
    expect_match(call$formula, "condition \\* time \\+ batch")
    expect_identical(colnames(call$feature.dat), rownames(call$meta.dat))
    expect_true(call$arguments$verbose)
  })
})

test_that("LinDA rejects missing planned coefficients", {
  skip_if_not_installed("MicrobiomeStat")
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) mock_linda_fit(omit = TRUE),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")

  expect_error(
    dar:::run_linda_model(
      rec, 0, 0, 0, FALSE, 0.03, FALSE, "pseudo-count", 0.5, 0.1,
      "BH", 0.05, 1L, FALSE
    ),
    class = "dar_error_invalid_result_contract"
  )
})

test_that("LinDA recovers known simple-model effect directions", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  result <- dar:::run_linda_model(
    rec, 0, 0, 0, FALSE, 0.03, FALSE, "pseudo-count", 0.5, 0.1,
    "BH", 0.05, 1L, FALSE
  )

  expect_complete_model_result(result, rec)
  expect_model_directions(result)
})

test_that("LinDA integrates with prep and the tidy result contract", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_linda(
      winsorize = FALSE, adaptive = FALSE, id = "linda__integration"
    )
  prepared <- prep(rec, parallel = FALSE)
  result <- tidy_results(prepared)

  expect_s4_class(prepared, "PrepRecipe")
  expect_true(all(result$method == "linda"))
  expect_true(all(
    result$effect_metric == "bias_corrected_log2_fold_change"
  ))
  expect_identical(unique(result$step_id), "linda__integration")
})

test_that("LinDA fits an estimable random-intercept model", {
  skip_if_not_installed("MicrobiomeStat")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + time + (1 | subject),
      targets = "condition", tax_level = "Species", time = "time"
    )
  result <- suppressWarnings(dar:::run_linda_model(
    rec, 0, 0, 0, FALSE, 0.03, FALSE, "pseudo-count", 0.5, 0.1,
    "BH", 0.05, 1L, FALSE
  ))

  expect_complete_model_result(result, rec)
})
