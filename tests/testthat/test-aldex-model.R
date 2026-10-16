capture_aldex_fits <- function() {
  fits <- list()
  stages <- character()
  native_stage_function <- dar:::engine_stage_function
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      stopifnot(identical(engine, "aldex2"))
      stages <<- c(stages, stage)
      if (!stage %in% c("clr", "fit")) {
        stop("The centralized ALDEx2 path must not call post-hoc effects.")
      }
      native <- native_stage_function(engine, stage)
      if (identical(stage, "fit")) {
        return(function(clr, ...) {
          fit <- native(clr = clr, ...)
          fits[[length(fits) + 1L]] <<- list(clr = clr, fit = fit)
          fit
        })
      }
      native
    },
    .package = "dar",
    .env = parent.frame()
  )
  list(fits = function() fits, stages = function() stages)
}

expect_aldex_native_contrasts <- function(rec, result, fits, alpha = 0.05) {
  compiled <- dar:::compile_model(rec, "aldex2")
  table <- dar:::flatten_model_output(result)
  expect_length(fits, nrow(compiled$contrasts))

  for (index in seq_len(nrow(compiled$contrasts))) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    parameterized <- dar:::reparameterize_model_contrast(compiled, contrast)
    captured <- fits[[index]]
    native <- captured$fit
    coefficient <- parameterized$coefficient
    observed <- table[table$contrast_id == contrast$contrast_id, , drop = FALSE]
    observed <- observed[match(rownames(native), observed$taxa_id), , drop = FALSE]

    expect_equal(
      observed$effect,
      native[[paste0(coefficient, ":Est")]] * parameterized$sign
    )
    expect_equal(observed$padj, native[[paste0(coefficient, ":pval.padj")]])
    expect_identical(observed$signif, !is.na(observed$padj) & observed$padj < alpha)
    expect_identical(
      observed$effect_metric, rep("aldex_glm_coefficient", nrow(observed))
    )

    # With the default Gaussian identity model, averaging the Monte Carlo
    # coefficients equals fitting their mean log-ratio response. Use the same
    # draws and the original contrast weights as an independent estimand check.
    instances <- ALDEx2::getMonteCarloInstances(captured$clr)
    response <- vapply(instances, rowMeans, numeric(nrow(native)))
    coefficients <- qr.coef(qr(compiled$matrix), t(response))
    expected <- as.numeric(crossprod(contrast$weights[[1]], coefficients))
    expect_equal(observed$effect, expected, tolerance = 1e-10)
  }
}

test_that("multigroup ALDEx2 reports each planned GLM coefficient", {
  skip_if_not_installed("ALDEx2")
  withr::local_seed(130)
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_aldex2(mc.samples = 8, id = "aldex2__coefficient")
  captured <- capture_aldex_fits()

  result <- suppressWarnings(dar:::run_aldex2_model(rec, 0.05, 8, "all", FALSE))

  expect_aldex_native_contrasts(rec, result, captured$fits())
  expect_setequal(captured$stages(), c("clr", "fit"))
  prepared <- dar:::prep_recipe(
    rec, list(aldex2__coefficient = result), list()
  )
  tidy <- tidy_results(prepared)
  compiled <- dar:::compile_model(rec, "aldex2")
  for (index in seq_len(nrow(compiled$contrasts))) {
    contrast <- compiled$contrasts[index, , drop = FALSE]
    coefficient <- dar:::reparameterize_model_contrast(compiled, contrast)$coefficient
    native <- captured$fits()[[index]]$fit
    rows <- tidy[tidy$contrast_id == contrast$contrast_id, , drop = FALSE]
    rows <- rows[match(rownames(native), rows$taxa_id), , drop = FALSE]
    expect_equal(rows$p_value, native[[paste0(coefficient, ":pval")]])
    expect_equal(rows$adj_p_value, native[[paste0(coefficient, ":pval.padj")]])
    expect_identical(rows$significant, rows$adj_p_value < 0.05)
    expect_identical(rows$effect_metric, rep("aldex_glm_coefficient", nrow(rows)))
  }
  path <- tempfile(fileext = ".rds")
  saveRDS(prepared, path)
  expect_identical(tidy_results(readRDS(path)), tidy)
})

test_that("adjusted ALDEx2 effects come from the covariate-adjusted GLM", {
  skip_if_not_installed("ALDEx2")
  withr::local_seed(131)
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species")
  captured <- capture_aldex_fits()

  result <- suppressWarnings(dar:::run_aldex2_model(rec, 0.05, 8, "all", FALSE))

  expect_aldex_native_contrasts(rec, result, captured$fits())
  expect_true("age" %in% colnames(dar:::compile_model(rec, "aldex2")$matrix))
})

test_that("ALDEx2 temporal effects include the planned interaction contrasts", {
  skip_if_not_installed("ALDEx2")
  withr::local_seed(132)
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time, targets = "condition", tax_level = "Species",
      time = "time", contrasts = "all"
    )
  captured <- capture_aldex_fits()

  result <- suppressWarnings(dar:::run_aldex2_model(rec, 0.05, 8, "all", FALSE))

  expect_aldex_native_contrasts(rec, result, captured$fits())
  expect_equal(
    sum(dar:::compile_model(rec, "aldex2")$contrasts$contrast_type ==
      "difference_in_differences"),
    3L
  )
})
test_that("reversing the ALDEx2 reference reverses effects, not significance", {
  skip_if_not_installed("ALDEx2")
  forward <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species")
  reverse <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + age, targets = "condition", tax_level = "Species",
      reference = c(condition = "treated")
    )
  forward_result <- withr::with_seed(133, suppressWarnings(
    dar:::run_aldex2_model(forward, 0.05, 8, "all", FALSE)
  ))$model
  reverse_result <- withr::with_seed(133, suppressWarnings(
    dar:::run_aldex2_model(reverse, 0.05, 8, "all", FALSE)
  ))$model

  expect_identical(forward_result$taxa_id, reverse_result$taxa_id)
  expect_identical(
    unique(as.character(reverse_result$contrast_id)), "condition[control-treated]"
  )
  expect_equal(reverse_result$effect, -forward_result$effect, tolerance = 1e-10)
  expect_equal(reverse_result$padj, forward_result$padj, tolerance = 1e-10)
  expect_identical(reverse_result$signif, forward_result$signif)
})

test_that("centralized ALDEx2 rejects the obsolete post-hoc effect stage", {
  expect_error(
    dar:::normalize_engine_args("aldex2", list(effect = list())),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args(
      "aldex2", list(effect = list(include.sample.summary = TRUE))
    ),
    class = "dar_error_invalid_engine_args"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  expect_error(
    dar:::run_aldex2_model(
      rec, 0.05, 8, "all", FALSE, engine_args = list(effect = list(CI = TRUE))
    ),
    class = "dar_error_invalid_engine_args"
  )
})
