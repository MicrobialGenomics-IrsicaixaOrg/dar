test_that("default engine designs do not produce model override warnings", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  constructors <- list(
    step_aldex, step_ancom, step_corncob, step_deseq,
    step_lefse, step_linda, step_maaslin, step_wilcox
  )

  for (constructor in constructors) {
    expect_silent(constructor(rec))
  }
  expect_silent(step_ancom(rec, rarefy = TRUE))
  expect_silent(step_ancom(rec, fix_formula = NULL))
  expect_silent(step_corncob(rec, formula_null = NULL, phi.formula_null = NULL))
})

test_that("rarefaction is silent during configuration and reported by prep", {
  skip_if_not_installed("ANCOMBC")
  withr::local_seed(130)
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")

  expect_silent(configured <- step_ancom(rec, rarefy = TRUE))
  expect_true(configured@steps[[1]]$rarefy)
  expect_silent(step_ancom(rec, rarefy = "no_seed"))
  expect_silent(step_ancom(rec, rarefy = FALSE))
  expect_message(
    prepared <- suppressWarnings(prep(configured, parallel = FALSE)),
    "not without controversy"
  )
  expect_s4_class(prepared, "PrepRecipe")
})

test_that("adding a model after default steps does not warn about overrides", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    step_ancom() |>
    step_corncob()

  expect_silent(
    add_model(rec, ~ condition, targets = "condition", tax_level = "Species")
  )
})

test_that("explicit engine designs still produce classified override warnings", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  designs <- list(
    function() step_ancom(rec, fix_formula = "condition"),
    function() step_ancom(rec, rand_formula = "(1 | subject)"),
    function() step_ancom(rec, group = "condition"),
    function() step_ancom(rec, global = TRUE),
    function() step_corncob(rec, formula_null = ~ 1),
    function() step_corncob(rec, phi.formula_null = ~ 1),
    function() step_maaslin(rec, random_effects = "subject"),
    function() step_maaslin(rec, reference = c("condition", "control")),
    function() step_lefse(rec, subclassCol = "batch")
  )

  for (design in designs) {
    expect_warning(design(), class = "dar_warning_model_override")
  }

  configured <- step_ancom(
    recipe(make_longitudinal_phy()), fix_formula = "condition"
  )
  expect_warning(
    add_model(configured, ~ condition, targets = "condition", tax_level = "Species"),
    class = "dar_warning_model_override"
  )
})

test_that("recipe files preserve default and explicit engine designs", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  phy <- make_longitudinal_phy()
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_ancom() |>
    step_corncob()
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)

  expect_silent(imported <- import_steps(recipe(phy), path))
  expect_identical(
    dar:::recipe_serialization_document(imported),
    dar:::recipe_serialization_document(rec)
  )

  configurations <- list(
    step_ancom(recipe(phy), fix_formula = "condition"),
    step_corncob(recipe(phy), formula_null = ~ 1),
    step_corncob(recipe(phy), phi.formula_null = ~ 1)
  )
  for (configured in configurations) {
    explicit <- suppressWarnings(
      add_model(configured, ~ condition, targets = "condition", tax_level = "Species")
    )
    export_steps(explicit, path)
    expect_warning(
      imported <- import_steps(recipe(phy), path),
      class = "dar_warning_model_override"
    )
    expect_identical(
      dar:::recipe_serialization_document(imported),
      dar:::recipe_serialization_document(explicit)
    )
  }
})
