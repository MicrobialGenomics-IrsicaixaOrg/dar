test_that("corncob thresholds absolute effects without changing estimates", {
  skip_if_not_installed("corncob")
  effects <- c(-3, -2, -0.2, 0.2, 2, 3)
  index <- 0L
  local_mocked_bindings(
    hessian = function(fit) diag(c(100, 100, 100)),
    .package = "corncob"
  )
  local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      function(...) {
        index <<- index + 1L
        list(
          np.mu = 2L, np.phi = 1L, b.mu = c(0, effects[[index]]),
          robust = FALSE, sep_da = FALSE
        )
      }
    },
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  run <- function(threshold, cutoff = 0.05) {
    index <<- 0L
    dar:::run_corncob_model(
      rec, ~ 1, "logit", "logit", FALSE, cutoff, "BH", threshold, FALSE
    )$model
  }
  default <- run(0)
  strict <- run(2)
  rejected <- run(2, cutoff = 0)
  expect_equal(strict$effect, effects)
  expect_identical(strict$signif, c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(
    dplyr::select(default, -signif), dplyr::select(strict, -signif)
  )
  expect_false(any(rejected$signif))
  expect_identical(strict$taxa_id, paste0("taxon_", seq_along(effects)))
})

test_that("corncob incompatible tests are reported rather than executed", {
  skip_if_not_installed("corncob")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  for (settings in list(list(boot = TRUE, B = 7L), list(test = "LRT"))) {
    configured <- suppressWarnings(do.call(step_corncob, c(list(rec), settings)))
    status <- dar:::model_steps_status(configured)
    expect_false(status$compatible)
    expect_match(status$reason, "non-bootstrap Wald")
  }
})
