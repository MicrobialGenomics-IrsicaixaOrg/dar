test_that("additive time is an adjustment whether or not it is declared", {
  for (time in list(NULL, "time")) {
    rec <- recipe(make_longitudinal_phy()) |>
      add_model(
        ~ condition + time, targets = "condition", tax_level = "Species",
        time = time
      )
    resolved <- dar:::resolve_model(rec)
    expect_identical(resolved$confounders, "time")

    for (engine in names(dar:::model_engine_capabilities)) {
      configured_step <- if (engine == "corncob") {
        dar:::step(engine, id = engine, test = "Wald", boot = FALSE)
      } else {
        dar:::step(engine, id = engine)
      }
      status <- dar:::model_step_status(rec, configured_step, resolved)
      supports_adjustment <- !engine %in% c("wilcox", "lefse")
      expect_identical(status$compatible, supports_adjustment, info = engine)
      if (!supports_adjustment) {
        expect_match(status$reason, "adjusted confounders", info = engine)
      }
    }
  }
})

test_that("simple inference on condition or time remains compatible", {
  phy <- make_longitudinal_phy()
  metadata <- as(phyloseq::sample_data(phy), "data.frame")
  metadata$time <- factor(metadata$time)
  phyloseq::sample_data(phy) <- phyloseq::sample_data(metadata)
  recipes <- list(
    recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Species"),
    recipe(phy) |>
      add_model(~ time, targets = "time", tax_level = "Species", time = "time")
  )

  for (rec in recipes) {
    expect_length(dar:::resolve_model(rec)$confounders, 0L)
    for (engine in c("wilcox", "lefse")) {
      expect_silent(
        configured <- dar:::add_step(rec, dar:::step(engine, id = engine))
      )
      status <- dar:::model_steps_status(configured)
      expect_true(status$compatible)
      expect_true(is.na(status$reason))
    }
  }
})

test_that("time-adjustment warnings are independent of construction order", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  constructors <- list(wilcox = step_wilcox, lefse = step_lefse)
  base <- recipe(make_longitudinal_phy())
  modeled <- add_model(
    base, ~ condition + time, targets = "condition", tax_level = "Species",
    time = "time"
  )

  for (engine in names(constructors)) {
    constructor <- constructors[[engine]]
    expect_warning(
      model_first <- constructor(modeled, id = engine),
      "adjusted confounders", class = "dar_warning_model_step_skipped"
    )
    step_first <- constructor(base, id = engine)
    expect_warning(
      step_first <- add_model(
        step_first, ~ condition + time, targets = "condition",
        tax_level = "Species", time = "time"
      ),
      "adjusted confounders", class = "dar_warning_model_step_skipped"
    )
    expect_identical(
      dar:::model_steps_status(model_first),
      dar:::model_steps_status(step_first)
    )
    expect_false(dar:::model_steps_status(model_first)$compatible)
  }
})

test_that("prep never executes group comparisons with time adjustment", {
  executed <- character()
  local_mocked_bindings(
    run_wilcox = function(...) {
      executed <<- c(executed, "wilcox")
      stop("An incompatible Wilcoxon step must not execute.")
    },
    run_lefse = function(...) {
      executed <<- c(executed, "lefse")
      stop("An incompatible LEfSe step must not execute.")
    },
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + time, targets = "condition", tax_level = "Species",
      time = "time"
    )
  for (engine in c("wilcox", "lefse")) {
    rec <- suppressWarnings(
      dar:::add_step(rec, dar:::step(engine, id = engine))
    )
  }

  expect_length(dar:::model_dependency_steps(rec), 0L)
  for (parallel in c(FALSE, TRUE)) {
    expect_warning(
      expect_error(
        prep(rec, parallel = parallel, workers = 1),
        class = "dar_error_no_compatible_steps"
      ),
      "adjusted confounders", class = "dar_warning_model_step_skipped"
    )
  }
  expect_length(executed, 0L)
})

test_that("prep keeps adjusted engines and records skipped group comparisons", {
  skip_if_not_installed("MicrobiomeStat")
  withr::local_seed(130)
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + time, targets = "condition", tax_level = "Species",
      time = "time"
    ) |>
    step_linda(prev_filter = 0, winsorize = FALSE, adaptive = FALSE,
               id = "linda__adjusted")
  expect_warning(
    rec <- step_wilcox(rec, id = "wilcox__skipped"),
    class = "dar_warning_model_step_skipped"
  )
  expect_warning(
    rec <- step_lefse(rec, id = "lefse__skipped"),
    class = "dar_warning_model_step_skipped"
  )
  expect_warning(
    prepared <- prep(rec, parallel = FALSE),
    "adjusted confounders", class = "dar_warning_model_step_skipped"
  )

  expect_identical(prepared@execution$executed_steps, "linda__adjusted")
  expect_identical(names(prepared@results), "linda__adjusted")
  expect_identical(steps_ids(prepared, "da"), "linda__adjusted")
  expect_identical(
    prepared@execution$skipped_steps$step_id,
    c("wilcox__skipped", "lefse__skipped")
  )
  expect_complete_model_result(prepared@results$linda__adjusted, rec)
  expect_identical(get_phy(prepared), get_phy(rec))
})
