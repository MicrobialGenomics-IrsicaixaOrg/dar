test_that("engine_args validates its nested stage contract", {
  expect_identical(dar:::normalize_engine_args("wilcox", list()), list())
  expect_error(
    dar:::normalize_engine_args("wilcox", 1),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args("wilcox", list(list(exact = FALSE))),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args("wilcox", list(unknown = list())),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args("wilcox", list(test = FALSE)),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args("wilcox", list(test = list(x = 1))),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args(
      "wilcox", list(test = list(max_significance = 0.1))
    ),
    class = "dar_error_invalid_engine_args"
  )
  expect_error(
    dar:::normalize_engine_args(
      "deseq", list(results = list(not_a_native_argument = TRUE))
    ),
    class = "dar_error_invalid_engine_args"
  )
})

test_that("engine_args is the final public argument of every DA step", {
  constructors <- list(
    step_deseq, step_aldex, step_ancom, step_corncob, step_maaslin,
    step_wilcox, step_lefse
  )
  purrr::walk(constructors, function(constructor) {
    expect_identical(tail(names(formals(constructor)), 1L), "engine_args")
  })
})

test_that("the stage registry covers all seven engines", {
  registry <- dar:::engine_args_registry()
  expect_identical(
    names(registry),
    c("deseq", "aldex", "ancom", "corncob", "maaslin", "wilcox", "lefse")
  )
  expect_identical(
    lapply(registry, function(engine) names(engine$stages)),
    list(
      deseq = c("size_factors", "fit", "results", "shrink"),
      aldex = c("clr", "fit", "effect"), ancom = "fit",
      corncob = "fit", maaslin = c("fit", "contrast"),
      wilcox = "test", lefse = c("transform", "fit")
    )
  )
  purrr::iwalk(registry, function(spec, engine) {
    purrr::walk(names(spec$stages), function(stage) {
      expect_true(is.function(dar:::engine_stage_function(engine, stage)))
    })
  })
})

test_that("native stages combine fixed, default and advanced arguments", {
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      function(x, option = 1, advanced = FALSE) {
        list(x = x, option = option, advanced = advanced)
      }
    },
    .package = "dar"
  )
  result <- dar:::exec_engine_stage(
    "wilcox", "test",
    engine_args = list(test = list(option = 3, advanced = TRUE)),
    fixed = list(x = 2), defaults = list(option = 1)
  )
  expect_identical(result, list(x = 2, option = 3, advanced = TRUE))
})

test_that("configured but unexecuted stages are rejected", {
  expect_error(
    dar:::check_unused_engine_args(
      "deseq", list(results = list(alpha = 0.1)), "shrink"
    ),
    class = "dar_error_invalid_engine_args"
  )
  expect_silent(
    dar:::check_unused_engine_args(
      "deseq", list(shrink = list(lfcThreshold = 1)), "shrink"
    )
  )
})

test_that("engine_args distinguishes repeated instances of one engine", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(
      id = "wilcox__one",
      engine_args = list(test = list(conf.int = FALSE))
    ) |>
    step_wilcox(
      id = "wilcox__two",
      engine_args = list(test = list(conf.int = TRUE))
    )

  expect_length(rec@steps, 2L)
  expect_false(identical(rec@steps[[1]]$engine_args, rec@steps[[2]]$engine_args))
})

test_that("legacy recipes reject non-empty engine_args during prep", {
  rec <- suppressWarnings(
    recipe(make_longitudinal_phy(), "condition", "Species")
  )
  rec <- dar:::add_step(rec, dar:::step(
    "wilcox", norm_method = "identity", max_significance = 0.05,
    p_adj_method = "BH", rarefy = FALSE, id = "wilcox__legacy",
    engine_args = list(test = list(conf.int = TRUE))
  ))

  expect_error(
    prep(rec, parallel = FALSE),
    class = "dar_error_engine_args_requires_model"
  )
  rec@steps[[1]]$engine_args <- list()
  expect_silent(dar:::validate_engine_args_model(rec))
})

test_that("nested engine_args survive export and import", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(
      id = "wilcox__roundtrip",
      engine_args = list(
        test = list(conf.int = TRUE, conf.level = 0.9, null.value = NULL)
      )
    )
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- import_steps(recipe(make_longitudinal_phy()), path)

  expect_identical(
    imported@steps[[1]]$engine_args,
    rec@steps[[1]]$engine_args
  )
})

test_that("step export rejects values that cannot be reconstructed", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step(
    "wilcox", id = "wilcox__environment",
    engine_args = list(test = list(callback_state = new.env()))
  ))

  expect_error(
    export_steps(rec, tempfile(fileext = ".json")),
    class = "dar_error_unserializable_step"
  )
})

test_that("recursive step serialization supports compatible values", {
  value <- list(
    fit = list(control = list(tolerance = 1e-4), transform = ~ x + y),
    effect = list(callback = function(x) x)
  )
  expression <- dar:::step_value_expr(value)
  restored <- eval(parse(text = expression))

  expect_equal(restored$fit$control$tolerance, 1e-4)
  expect_equal(deparse(restored$fit$transform), deparse(~ x + y))
  expect_true(is.function(restored$effect$callback))

  threshold <- 2
  expect_error(
    dar:::step_value_expr(list(fit = list(callback = function(x) x > threshold))),
    class = "dar_error_unserializable_step"
  )
})
