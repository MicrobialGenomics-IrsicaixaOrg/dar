test_that("every DA engine requires add_model before any execution", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    required_deps = function(...) stop("Dependency checks must not run"),
    run_filter_taxa = function(...) stop("Preprocessing must not run"),
    engine_stage_function = function(...) stop("Native engines must not run"),
    .package = "dar"
  )
  withr::local_seed(130)
  phy <- make_longitudinal_phy()
  engines <- c(
    "aldex2", "aldex3", "ancom", "corncob", "deseq", "lefse",
    "linda", "maaslin", "wilcox"
  )
  for (selectors in c(FALSE, TRUE)) {
    base <- if (selectors) {
      suppressWarnings(recipe(phy, "condition", "Species"))
    } else {
      recipe(phy)
    }
    base <- step_filter_taxa(base, function(x) sum(x) > 0)
    for (engine in engines) {
      constructor <- getExportedValue("dar", paste0("step_", engine))
      rec <- constructor(base, id = paste0(engine, "__required"))
      seed <- .Random.seed
      expect_error(
        prep(rec, parallel = FALSE),
        "add_model",
        class = "dar_error_model_required",
        info = engine
      )
      runner <- get(paste0("run_", engine), asNamespace("dar"))
      # Other arguments must stay unevaluated until the model is checked.
      expect_error(runner(rec), class = "dar_error_model_required", info = engine)
      expect_error(
        dar:::compile_model(rec, engine),
        class = "dar_error_model_required",
        info = engine
      )
      expect_identical(.Random.seed, seed)
      expect_identical(get_phy(rec), phy)
      expect_null(get_model(rec))
      configured <- add_model(
        rec, ~ condition, targets = "condition", tax_level = "Species"
      )
      expect_identical(get_model(configured)$targets, "condition")
      expect_length(configured@steps, 2L)
    }
  }
})

test_that("removing a configured model cannot restore model-free DA", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox() |>
    add_model(NULL)
  expect_null(get_model(rec))
  expect_error(prep(rec, parallel = FALSE), class = "dar_error_model_required")
})

test_that("preprocessing alone remains executable without a model", {
  rec <- recipe(make_longitudinal_phy()) |>
    step_filter_taxa(function(x) sum(x) > 0)
  prepared <- prep(rec, parallel = FALSE)
  expect_s4_class(prepared, "PrepRecipe")
  expect_null(get_model(prepared))
  expect_identical(get_phy(prepared), get_phy(rec))
})

test_that("import cannot execute DA without an explicit model", {
  skip_if_not_installed("microbiome")
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- step_wilcox(recipe(make_longitudinal_phy()))
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- import_steps(recipe(make_longitudinal_phy()), path)
  expect_null(get_model(imported))
  expect_error(prep(imported, parallel = FALSE), class = "dar_error_model_required")
  configured <- add_model(
    imported, ~ condition, targets = "condition", tax_level = "Species"
  )
  expect_s4_class(prep(configured, parallel = FALSE), "PrepRecipe")
})
