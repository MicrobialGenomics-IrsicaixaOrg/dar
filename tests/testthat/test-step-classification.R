test_that("every supported constructor has a class-based execution role", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )

  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")

  preprocessing <- list(
    filter_by_abundance = function(x, id) step_filter_by_abundance(x, id = id),
    filter_by_prevalence = function(x, id) step_filter_by_prevalence(x, id = id),
    filter_by_rarity = function(x, id) step_filter_by_rarity(x, id = id),
    filter_by_variance = function(x, id) step_filter_by_variance(x, id = id),
    filter_taxa = function(x, id) step_filter_taxa(x, function(value) TRUE, id),
    subset_taxa = function(x, id) {
      step_subset_taxa(x, tax_level = "Kingdom", taxa = "Bacteria", id = id)
    },
    rarefaction = function(x, id) step_rarefaction(x, id = id)
  )
  da <- list(
    aldex = function(x, id) step_aldex(x, id = id),
    ancom = function(x, id) step_ancom(x, id = id),
    corncob = function(x, id) step_corncob(x, id = id),
    deseq = function(x, id) step_deseq(x, id = id),
    lefse = function(x, id) suppressWarnings(step_lefse(x, id = id)),
    linda = function(x, id) step_linda(x, id = id),
    maaslin = function(x, id) step_maaslin(x, id = id),
    wilcox = function(x, id) step_wilcox(x, id = id)
  )

  purrr::iwalk(preprocessing, function(constructor, method) {
    configured <- constructor(rec, paste0("arbitrary-preprocessing-", method))
    configured_step <- configured@steps[[1L]]
    expect_identical(dar:::step_method(configured_step), method)
    expect_identical(dar:::step_role(configured_step), "preprocessing")
    expect_true(is.function(dar:::step_runner(configured_step)))
  })
  purrr::iwalk(da, function(constructor, method) {
    configured <- suppressWarnings(
      constructor(rec, paste0("filter_named_da_", method))
    )
    configured_step <- configured@steps[[1L]]
    expect_identical(dar:::step_method(configured_step), method)
    expect_identical(dar:::step_role(configured_step), "da")
    expect_identical(dar:::step_engine(configured_step), method)
    expect_true(is.function(dar:::step_runner(configured_step)))
  })
  registered_da <- purrr::keep(
    dar:::recipe_step_registry,
    ~ identical(.x$role, "da")
  )
  expect_setequal(
    purrr::map_chr(registered_da, "engine"),
    names(dar:::model_engine_capabilities)
  )

  prepared <- dar:::prep_recipe(rec, list(), list(), execution = list())
  baked <- bake(prepared, id = "looks_like_filter")
  expect_identical(dar:::step_role(baked@bakes[[1L]]), "bake")
  expect_identical(dar:::step_method(baked@bakes[[1L]]), "bake")
})

test_that("misleading IDs never change classification", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_filter_by_prevalence(id = "custom_preprocessing") |>
    step_deseq(id = "filter_named_da") |>
    step_rarefaction(id = "neutral_operation")

  expect_identical(
    steps_ids(rec, "prepro"),
    c("custom_preprocessing", "neutral_operation")
  )
  expect_identical(steps_ids(rec, "da"), "filter_named_da")
  expect_true(contains_rarefaction(rec))
  expect_identical(
    dar:::model_steps_status(rec)$step_id,
    "filter_named_da"
  )

  printed <- capture.output(show(rec))
  prepro_line <- grep("custom_preprocessing", printed, value = TRUE)
  da_line <- grep("filter_named_da", printed, value = TRUE)
  expect_match(prepro_line, "step_filter_by_prevalence")
  expect_match(da_line, "step_deseq")
  expect_true(any(grepl("Preprocessing steps", printed, fixed = TRUE)))
})

test_that("preprocessing steps execute in configured order", {
  observed <- character()
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    required_deps = function(...) invisible(NULL),
    run_rarefaction = function(rec, id) {
      observed <<- c(observed, id)
      rec
    },
    run_filter_by_prevalence = function(rec, threshold, id) {
      observed <<- c(observed, id)
      rec
    },
    run_subset_taxa = function(rec, tax_level, taxa) {
      observed <<- c(observed, "third_subset")
      rec
    },
    .package = "dar"
  )

  rec <- recipe(make_longitudinal_phy()) |>
    step_rarefaction(id = "first_neutral") |>
    step_filter_by_prevalence(id = "second_neutral") |>
    step_subset_taxa(
      tax_level = "Kingdom", taxa = "Bacteria", id = "third_subset"
    )
  prepared <- prep(rec, parallel = FALSE)

  expect_identical(
    observed,
    c("first_neutral", "second_neutral", "third_subset")
  )
  expect_identical(steps_ids(prepared, "prepro"), observed)
  expect_length(steps_ids(prepared, "da"), 0L)
})

test_that("dependency checks use the same compatible DA classification", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + (1 | subject),
      targets = "condition", tax_level = "Species"
    )
  rec <- suppressWarnings(
    dar:::add_step(
      rec,
      dar:::step(
        "filter_by_prevalence", threshold = 0.1, id = "neutral_preprocessing"
      )
    ) |>
      dar:::add_step(dar:::step("deseq", id = "filter_named_da")) |>
      dar:::add_step(dar:::step("ancom", id = "neutral_da"))
  )

  status <- dar:::model_steps_status(rec)
  dependency_steps <- dar:::model_dependency_steps(rec, status)
  dependency_ids <- purrr::map_chr(dependency_steps, "id")

  expect_false(status$compatible[status$step_id == "filter_named_da"])
  expect_true(status$compatible[status$step_id == "neutral_da"])
  expect_identical(
    dependency_ids,
    c("neutral_preprocessing", "neutral_da")
  )
})

test_that("step validation rejects malformed and misplaced objects", {
  rec <- recipe(make_longitudinal_phy())
  unknown <- structure(list(id = "unknown"), class = c("step_unknown", "step"))
  conflicting <- structure(
    list(id = "conflicting"),
    class = c("step_deseq", "step_filter_taxa", "step")
  )

  expect_error(dar:::add_step(rec, list(id = "plain")),
               class = "dar_error_invalid_step")
  expect_error(dar:::add_step(rec, unknown), class = "dar_error_invalid_step")
  expect_error(
    dar:::add_step(rec, dar:::step("deseq")),
    class = "dar_error_invalid_step"
  )
  expect_error(dar:::add_step(rec, conflicting),
               class = "dar_error_invalid_step")
  expect_error(
    dar:::add_step(rec, dar:::step("bake", id = "bake_in_recipe")),
    class = "dar_error_invalid_step"
  )

  rec <- dar:::add_step(
    rec,
    dar:::step("filter_by_prevalence", threshold = 0.1, id = "duplicate")
  )
  expect_error(
    dar:::add_step(
      rec,
      dar:::step("filter_by_prevalence", threshold = 0.2, id = "duplicate")
    ),
    class = "dar_error_invalid_step"
  )

  prepared <- dar:::prep_recipe(rec, list(), list(), execution = list())
  expect_error(
    dar:::add_step(prepared, dar:::step("deseq", id = "da_in_bakes")),
    class = "dar_error_invalid_step"
  )
})

test_that("direct step-slot corruption is detected by Recipe validity", {
  rec <- recipe(make_longitudinal_phy())
  corrupt <- rec
  corrupt@steps <- list(
    structure(list(id = "unknown"), class = c("step_unknown", "step"))
  )
  expect_match(methods::validObject(corrupt, test = TRUE), "supported")
  expect_error(steps_ids(corrupt), class = "dar_error_invalid_recipe")

  duplicated <- rec
  duplicated@steps <- list(
    dar:::step("deseq", id = "same"),
    dar:::step("wilcox", id = "same")
  )
  expect_match(methods::validObject(duplicated, test = TRUE), "unique")

  prepared <- dar:::prep_recipe(rec, list(), list(), execution = list())
  prepared@bakes <- list(dar:::step("deseq", id = "wrong_role"))
  expect_match(methods::validObject(prepared, test = TRUE), "expected one of: bake")
})
