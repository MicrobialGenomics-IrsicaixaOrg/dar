reserved_da_ids <- c(
  "taxa_id", "contrast_id", "comparison", "contrast_type", "var", "effect",
  "sum", ".tax_order", ".contrast_order", ".effect_order"
)

test_that("DA constructors reject reserved IDs before modifying recipes", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  constructors <- list(
    step_aldex2, step_aldex3, step_ancom, step_corncob, step_deseq, step_lefse,
    step_linda, step_maaslin, step_rademu, step_wilcox
  )
  for (constructor in constructors) {
    for (id in reserved_da_ids) {
      expect_error(constructor(rec, id = id), "reserved", class = "dar_error_invalid_step")
    }
  }
  expect_length(rec@steps, 0L)
  expect_error(
    step_wilcox(recipe(make_multilevel_phy()), id = "effect"),
    "da_effect", class = "dar_error_invalid_step"
  )
})

test_that("reserved DA IDs are detected in supplied and corrupted objects", {
  phy <- make_multilevel_phy()
  rec <- recipe(phy)
  prepared <- dar:::prep_recipe(rec, list(), list())
  for (id in reserved_da_ids) {
    configured_step <- dar:::step("wilcox", id = id)
    expect_error(
      recipe(phy, steps = list(configured_step)),
      "reserved"
    )
    corrupt <- rec
    corrupt@steps <- list(configured_step)
    expect_match(methods::validObject(corrupt, test = TRUE), "reserved")
    expect_error(prep(corrupt, parallel = FALSE), class = "dar_error_invalid_recipe")
    prepared@steps <- list(configured_step)
    expect_match(methods::validObject(prepared, test = TRUE), "reserved")
    expect_error(tidy_results(prepared), class = "dar_error_invalid_recipe")
  }
})

test_that("reserved names remain valid for preprocessing and bake IDs", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy())
  prepared <- dar:::prep_recipe(rec, list(), list())
  for (id in reserved_da_ids) {
    filtered <- step_filter_by_prevalence(rec, id = id)
    expect_identical(steps_ids(filtered, "prepro"), id)
    expect_true(methods::validObject(filtered, test = TRUE))
    baked <- bake(prepared, id = id)
    expect_identical(baked@bakes[[1]]$id, id)
    expect_true(methods::validObject(baked, test = TRUE))
  }
})

test_that("JSON imports reject reserved DA IDs transactionally with context", {
  local_mocked_bindings(
    recipes_pkg_check = function(...) invisible(NULL),
    .package = "dar"
  )
  original <- recipe(make_multilevel_phy())
  configured <- original |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_filter_by_prevalence(id = "filter__safe") |>
    step_wilcox(id = "wilcox__safe")
  path <- withr::local_tempfile(fileext = ".json")
  export_steps(configured, path)
  document <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  applied <- 0L
  local_mocked_bindings(
    apply_serialized_step = function(...) {
      applied <<- applied + 1L
      stop("Constructors must not execute before all records are validated.")
    },
    .package = "dar"
  )
  for (id in reserved_da_ids) {
    document$steps[[2]]$id <- id
    jsonlite::write_json(document, path, auto_unbox = TRUE, null = "null")
    error <- expect_error(
      import_steps(original, path), "reserved",
      class = "dar_error_invalid_recipe_record"
    )
    expect_identical(error$record, "steps[[2]]")
    expect_identical(error$field, "id")
  }
  expect_identical(applied, 0L)
  expect_length(original@steps, 0L)
  expect_null(get_model(original))
})

test_that("safe IDs preserve consensus columns and overlap calculations", {
  safe_ids <- c(
    ".custom", "wilcox_effect", "EFFECT", "name", "value", "method",
    "taxa", "method_count", "sum_methods"
  )
  for (id in safe_ids) {
    rec <- recipe(make_multilevel_phy()) |>
      add_model(~ condition, targets = "condition", tax_level = "Species")
    rec <- dar:::add_step(rec, dar:::step("wilcox", id = id))
    rec <- dar:::add_step(rec, dar:::step("lefse", id = "reference"))
    plan <- dar:::resolve_model(rec)$contrast_plan[1, ]
    rows <- tibble::tibble(
      taxa_id = c("taxon_1", "taxon_2"), taxa = c("Species_1", "Species_2"),
      contrast_id = plan$contrast_id, comparison = plan$comparison,
      contrast_type = plan$contrast_type, var = plan$var,
      effect = 1, padj = 0.01, signif = TRUE
    )
    prepared <- dar:::prep_recipe(
      rec,
      stats::setNames(
        list(list(model = rows[1, ]), list(model = rows[2, ])),
        c(id, "reference")
      ),
      list()
    )
    wide <- intersection_df(prepared, contrast_id = plan$contrast_id)
    expect_true(id %in% names(wide))
    expect_identical(wide[[id]], c(1L, rep(0L, 5)))
    expect_false(any(c(".tax_order", ".contrast_order", ".effect_order") %in% names(wide)))
    tidy <- intersection_df(prepared, tidy = TRUE, contrast_id = plan$contrast_id)
    expect_setequal(tidy$name, c(id, "reference"))
    overlap <- overlap_df(prepared, type = "da", contrast_id = plan$contrast_id)
    expect_equal(overlap[id, "reference"], 0)
    expect_equal(overlap[id, id], 1)
    expect_setequal(find_intersections(prepared)$step_ids, c(id, "reference"))
  }
})

test_that("dotted custom IDs survive native execution and serialization", {
  skip_if_not_installed("microbiome")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(id = ".custom")
  path <- withr::local_tempfile(fileext = ".json")
  export_steps(rec, path)
  restored <- import_steps(recipe(make_multilevel_phy()), path)
  expect_identical(restored@steps[[1]]$id, ".custom")
  prepared <- prep(restored, parallel = FALSE)
  expect_identical(steps_ids(prepared, "da"), ".custom")
  expect_identical(unique(tidy_results(prepared)$step_id), ".custom")
  wide <- intersection_df(prepared)
  expect_identical(sum(wide$.custom), sum(tidy_results(prepared)$significant))
  expect_true(is.numeric(overlap_df(prepared)[".custom", ".custom"]))
  expect_identical(get_phy(prepared), get_phy(rec))
})
