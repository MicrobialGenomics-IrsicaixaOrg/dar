write_recipe_json <- function(document, path) {
  jsonlite::write_json(
    document, path, auto_unbox = TRUE, pretty = TRUE,
    null = "null", na = "null", digits = NA
  )
  path
}

read_recipe_json <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

test_that("schema v1 round trips model-free preprocessing recipes", {
  data(metaHIV_phy)
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(metaHIV_phy) |>
    step_subset_taxa(
      tax_level = "Kingdom", taxa = c("Bacteria", "Archaea"),
      id = "subset__safe"
    ) |>
    step_filter_taxa(
      function(x) sum(x > 0, na.rm = TRUE) >= 0.3 * length(x),
      id = "filter__safe"
    ) |>
    step_filter_by_prevalence(threshold = 0.2, id = "prevalence__safe")

  path <- tempfile(fileext = ".json")
  expect_identical(export_steps(rec, path), invisible(path))
  raw <- read_recipe_json(path)
  expect_identical(raw$format, "dar-recipe")
  expect_equal(raw$version, 1)
  expect_null(raw$model)

  imported <- import_steps(recipe(metaHIV_phy), path)
  expect_identical(
    dar:::recipe_serialization_document(imported),
    dar:::recipe_serialization_document(rec)
  )
})

test_that("schema v1 preserves centralized models and nested data arguments", {
  rec <- suppressWarnings(recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch + (1 | subject),
      targets = "condition", tax_level = "Species", time = "time",
      reference = c(condition = "control"), contrasts = "all",
      na_action = "drop", max_time_levels = 8
    ) |>
    step_wilcox(
      id = "wilcox__nested",
      engine_args = list(
        test = list(conf.int = TRUE, conf.level = 0.9, null.value = NULL)
      )
    ))
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- suppressWarnings(import_steps(recipe(make_longitudinal_phy()), path))

  expect_equal(deparse(get_model(imported)$formula), deparse(get_model(rec)$formula))
  expect_identical(get_model(imported)[-1], get_model(rec)[-1])
  expect_identical(imported@steps[[1]]$engine_args, rec@steps[[1]]$engine_args)
})

test_that("every registered constructor round trips through the allowlist", {
  data(metaHIV_phy)
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(metaHIV_phy) |>
    add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
  constructors <- list(
    function(x) step_filter_by_abundance(x, id = "abundance__x"),
    function(x) step_filter_by_prevalence(x, id = "prevalence__x"),
    function(x) step_filter_by_rarity(x, id = "rarity__x"),
    function(x) step_filter_by_variance(x, id = "variance__x"),
    function(x) step_filter_taxa(x, function(x) x > -1, id = "filter__x"),
    function(x) step_subset_taxa(x, "Kingdom", "Bacteria", id = "subset__x"),
    function(x) step_rarefaction(x, id = "rarefaction__x"),
    function(x) step_aldex(x, id = "aldex__x"),
    function(x) step_ancom(x, id = "ancom__x"),
    function(x) step_corncob(x, id = "corncob__x"),
    function(x) step_deseq(x, id = "deseq__x"),
    function(x) suppressWarnings(step_lefse(x, id = "lefse__x")),
    function(x) step_linda(x, id = "linda__x"),
    function(x) step_maaslin(x, id = "maaslin__x"),
    function(x) step_wilcox(x, id = "wilcox__x")
  )
  rec <- suppressWarnings(Reduce(function(x, constructor) constructor(x), constructors, rec))
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- suppressWarnings(import_steps(recipe(metaHIV_phy), path))

  expect_identical(
    dar:::recipe_serialization_document(imported),
    dar:::recipe_serialization_document(rec)
  )
})

test_that("bake imports require explicit preparation", {
  data(metaHIV_phy)
  rec <- recipe(metaHIV_phy) |>
    add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
  prepared <- dar:::prep_recipe(rec, results = list(), bakes = list()) |>
    bake(id = "bake__safe")
  path <- tempfile(fileext = ".json")
  export_steps(prepared, path)

  expect_error(
    import_steps(recipe(metaHIV_phy), path),
    class = "dar_error_import_requires_prepare"
  )

  prep_calls <- 0L
  local_mocked_bindings(
    prep = function(rec, ...) {
      prep_calls <<- prep_calls + 1L
      dar:::prep_recipe(rec, results = list(), bakes = list())
    },
    .package = "dar"
  )
  imported <- import_steps(recipe(metaHIV_phy), path, prepare = TRUE)
  expect_s4_class(imported, "PrepRecipe")
  expect_equal(prep_calls, 1L)
  expect_identical(imported@bakes[[1]]$id, "bake__safe")

  prep_calls <- 0L
  imported_prepared <- import_steps(
    dar:::prep_recipe(rec, results = list(), bakes = list()), path
  )
  expect_equal(prep_calls, 0L)
  expect_identical(imported_prepared@bakes[[1]]$id, "bake__safe")
})

test_that("legacy and malformed files are rejected without execution", {
  data(metaHIV_phy)
  legacy <- tempfile(fileext = ".json")
  writeLines(c(
    "{",
    '  \"value\": [c(\"options(dar.import.executed = TRUE)\")],',
    '  \"id\": [c(\"filter_taxa__legacy\")],',
    "}"
  ), legacy)
  old_option <- getOption("dar.import.executed")
  on.exit(options(dar.import.executed = old_option), add = TRUE)
  options(dar.import.executed = FALSE)

  expect_error(
    import_steps(recipe(metaHIV_phy), legacy),
    class = "dar_error_legacy_recipe_schema"
  )
  expect_false(getOption("dar.import.executed"))

  malformed <- tempfile(fileext = ".json")
  writeLines("{not-json", malformed)
  expect_error(
    import_steps(recipe(metaHIV_phy), malformed),
    class = "dar_error_invalid_recipe_file"
  )
})

test_that("formula and predicate payloads cannot execute arbitrary calls", {
  data(metaHIV_phy)
  rec <- recipe(metaHIV_phy) |>
    add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
    step_filter_taxa(function(x) sum(x) > 0, id = "filter__safe")
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  document <- read_recipe_json(path)

  marker <- tempfile()
  document$model$formula$text <- paste0("~ system(", deparse(paste("touch", marker)), ")")
  write_recipe_json(document, path)
  expect_error(
    import_steps(recipe(metaHIV_phy), path),
    class = "dar_error_unsafe_recipe_value"
  )
  expect_false(file.exists(marker))

  export_steps(rec, path)
  document <- read_recipe_json(path)
  document$steps[[1]]$arguments$.f$body <- list(
    node = "call", operator = "system",
    arguments = list(list(name = NULL, value = list(node = "symbol", name = "x")))
  )
  write_recipe_json(document, path)
  expect_error(
    import_steps(recipe(metaHIV_phy), path),
    class = "dar_error_unsafe_recipe_value"
  )
})

test_that("unknown schema content fails transactionally with context", {
  data(metaHIV_phy)
  original <- recipe(metaHIV_phy)
  configured <- original |>
    step_filter_by_prevalence(id = "prevalence__safe")
  path <- tempfile(fileext = ".json")
  export_steps(configured, path)
  document <- read_recipe_json(path)
  document$steps[[1]]$constructor <- "step_system"
  write_recipe_json(document, path)

  error <- expect_error(
    import_steps(original, path),
    class = "dar_error_invalid_recipe_record"
  )
  expect_identical(error$record, "steps[[1]]")
  expect_identical(error$field, "constructor")
  expect_length(original@steps, 0L)

  document <- read_recipe_json(export_steps(configured, path))
  document$unexpected <- TRUE
  write_recipe_json(document, path)
  expect_error(
    import_steps(original, path),
    class = "dar_error_invalid_recipe_record"
  )

  document <- read_recipe_json(export_steps(configured, path))
  document$steps[[1]]$arguments$unexpected <- list(
    type = "logical", values = list(TRUE), missing = list(FALSE), names = NULL
  )
  write_recipe_json(document, path)
  expect_error(
    import_steps(original, path),
    class = "dar_error_invalid_recipe_record"
  )

  document <- read_recipe_json(export_steps(configured, path))
  document$steps[[1]]$id <- ""
  write_recipe_json(document, path)
  expect_error(
    import_steps(original, path),
    class = "dar_error_invalid_recipe_record"
  )

  document <- read_recipe_json(export_steps(configured, path))
  document$version <- 999
  write_recipe_json(document, path)
  expect_error(
    import_steps(original, path),
    class = "dar_error_unsupported_recipe_schema"
  )

  duplicate <- tempfile(fileext = ".json")
  writeLines(
    '{"format":"dar-recipe","format":"dar-recipe","version":1,"model":null,"steps":[],"bakes":[]}',
    duplicate
  )
  expect_error(
    import_steps(original, duplicate),
    class = "dar_error_invalid_recipe_record"
  )
})

test_that("function values in engine_args cannot be exported", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step(
    "wilcox", id = "wilcox__callback",
    engine_args = list(test = list(callback = function(x) x))
  ))
  expect_error(
    export_steps(rec, tempfile(fileext = ".json")),
    class = "dar_error_unserializable_step"
  )
})

test_that("filter predicates enforce the safe single-expression subset", {
  data(metaHIV_phy)
  expect_s4_class(
    recipe(metaHIV_phy) |>
      step_filter_taxa("function(x) sum(x > 0) >= 0.3 * length(x)"),
    "Recipe"
  )
  threshold <- 2
  expect_error(
    recipe(metaHIV_phy) |>
      step_filter_taxa(function(x) sum(x) > threshold),
    class = "dar_error_invalid_filter_predicate"
  )
  expect_error(
    recipe(metaHIV_phy) |>
      step_filter_taxa(function(x) system("true")),
    class = "dar_error_invalid_filter_predicate"
  )
  expect_error(
    recipe(metaHIV_phy) |>
      step_filter_taxa(function(x) { sum(x) > 0 }),
    class = "dar_error_invalid_filter_predicate"
  )
})
