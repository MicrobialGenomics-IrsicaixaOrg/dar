test_that("add_model owns targets and taxonomic resolution", {
  phy <- make_longitudinal_phy()
  rec <- recipe(phy)
  modeled <- add_model(
    rec,
    ~ condition * time + batch + (1 | subject),
    targets = "condition",
    tax_level = "species",
    time = "time",
    reference = c(condition = "control", time = "0")
  )

  expect_equal(get_model(modeled)$targets, "condition")
  expect_equal(get_model(modeled)$tax_level, "Species")
  expect_equal(modeled@var_info$vars, "condition")
  expect_equal(modeled@tax_info$tax_lev, "Species")
  expect_equal(get_model(modeled)$time, "time")
  expect_equal(get_model(modeled)$contrasts, "simple")
  expect_true(methods::validObject(modeled, test = TRUE))
  expect_null(get_model(add_model(modeled, NULL)))
})

test_that("add_model validates formula variables, targets and references", {
  phy <- make_longitudinal_phy()
  rec <- recipe(phy)

  add <- function(formula, targets = "condition", tax_level = "Species", ...) {
    add_model(rec, formula, targets = targets, tax_level = tax_level, ...)
  }

  expect_error(add(response ~ condition), class = "dar_error_invalid_model")
  expect_error(add(~ missing), class = "dar_error_missing_vars")
  expect_error(add(~ batch), class = "dar_error_invalid_model")
  expect_error(add(~ condition, targets = NULL), class = "dar_error_invalid_model")
  expect_error(add(~ condition, targets = c("condition", "condition")),
               class = "dar_error_invalid_model")
  expect_error(add(~ condition, targets = "missing"), class = "dar_error_missing_vars")
  expect_error(add(~ age, targets = "age"), class = "dar_error_invalid_model")
  expect_error(add(~ condition, tax_level = NULL), class = "dar_error_invalid_model")
  expect_error(add(~ condition, tax_level = "Missing_rank"),
               class = "dar_error_missing_tax")
  expect_error(
    add(~ condition, reference = c(condition = "missing")),
    class = "dar_error_invalid_model"
  )
  expect_error(
    add(~ condition + age, reference = c(age = "35")),
    class = "dar_error_invalid_model"
  )
  expect_error(
    add(~ condition:batch),
    class = "dar_error_invalid_model"
  )
  expect_error(
    add(~ scale(age) + condition),
    class = "dar_error_invalid_model"
  )
})

test_that("legacy selectors remain compatible with classified warnings", {
  phy <- make_longitudinal_phy()
  expect_warning(
    legacy <- recipe(phy, var_info = "condition", tax_info = "Species"),
    class = "dar_warning_deprecated_recipe_argument"
  )
  expect_warning(
    legacy <- add_var(legacy, "batch"),
    class = "dar_warning_deprecated_selector"
  )
  expect_warning(
    legacy <- add_tax(legacy, "species"),
    class = "dar_warning_deprecated_selector"
  )
  expect_warning(get_var(legacy), class = "dar_warning_deprecated_selector")
  expect_warning(get_tax(legacy), class = "dar_warning_deprecated_selector")
  expect_null(get_model(legacy))
  expect_true(methods::validObject(legacy, test = TRUE))
})

test_that("legacy centralized models are normalized on access", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec@model$targets <- NULL
  rec@model$tax_level <- NULL

  expect_true(methods::validObject(rec, test = TRUE))
  expect_equal(get_model(rec)$targets, "condition")
  expect_equal(get_model(rec)$tax_level, "Species")
})

test_that("deprecated setters cannot diverge from a centralized model", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + batch,
      targets = "condition",
      tax_level = "Species"
    )

  expect_warning(
    rec <- add_var(rec, "batch"),
    class = "dar_warning_deprecated_selector"
  )
  expect_equal(get_model(rec)$targets, "batch")
  expect_equal(rec@var_info$vars, "batch")

  expect_warning(
    rec <- add_tax(rec, "species"),
    class = "dar_warning_deprecated_selector"
  )
  expect_equal(get_model(rec)$tax_level, "Species")
  expect_equal(rec@tax_info$tax_lev, "Species")

  corrupt <- rec
  corrupt@var_info$vars <- "condition"
  expect_match(methods::validObject(corrupt, test = TRUE), "must agree")
})

test_that("removing a model preserves legacy selectors for migration", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    add_model(NULL)

  expect_null(get_model(rec))
  expect_equal(rec@var_info$vars, "condition")
  expect_equal(rec@tax_info$tax_lev, "Species")
})

test_that("model validity detects direct slot corruption and is inherited", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  corrupt <- rec
  corrupt@model$contrasts <- "unsupported"

  expect_match(methods::validObject(corrupt, test = TRUE), "contrasts")
  prepped <- dar:::prep_recipe(rec, list(), list(), execution = list())
  expect_true(methods::validObject(prepped, test = TRUE))
})

test_that("temporal contrast plan contains both simple-effect directions", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time"
    )
  resolved <- dar:::resolve_model(rec)
  plan <- resolved$contrast_plan

  expect_equal(sum(plan$contrast_type == "condition_within_time"), 3)
  expect_equal(sum(plan$contrast_type == "time_within_condition"), 6)
  expect_false(any(plan$contrast_type == "difference_in_differences"))
  expect_equal(levels(resolved$data$time), c("0", "1", "2"))
  expect_true(all(lengths(plan$weights) == ncol(resolved$matrix)))

  all_rec <- add_model(
    recipe(make_longitudinal_phy()),
    ~ condition * time + batch,
    targets = "condition",
    tax_level = "Species",
    time = "time",
    contrasts = "all"
  )
  all_plan <- dar:::resolve_model(all_rec)$contrast_plan
  expect_equal(sum(all_plan$contrast_type == "difference_in_differences"), 3)
})

test_that("temporal contrast weights reproduce an independent algebraic truth", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition",
      tax_level = "Species",
      time = "time",
      contrasts = "all"
    )
  compiled <- dar:::compile_model(rec, "deseq")
  plan <- compiled$contrasts

  coefficients <- stats::setNames(
    rep(0, ncol(compiled$matrix)), colnames(compiled$matrix)
  )
  coefficients[c(
    "conditiontreated", "time1", "time2", "batchB",
    "conditiontreated:time1", "conditiontreated:time2"
  )] <- c(2, 1, 3, 17, 1, 3)

  expected <- c(
    "condition[treated-control]@time[0]" = 2,
    "condition[treated-control]@time[1]" = 3,
    "condition[treated-control]@time[2]" = 5,
    "time[1-0]@condition[control]" = 1,
    "time[2-0]@condition[control]" = 3,
    "time[2-1]@condition[control]" = 2,
    "time[1-0]@condition[treated]" = 2,
    "time[2-0]@condition[treated]" = 6,
    "time[2-1]@condition[treated]" = 4,
    "condition[treated-control]:time[1-0]" = 1,
    "condition[treated-control]:time[2-0]" = 3,
    "condition[treated-control]:time[2-1]" = 2
  )
  observed <- vapply(
    plan$weights,
    function(weights) sum(weights * coefficients),
    numeric(1)
  )
  names(observed) <- as.character(plan$contrast_id)

  expect_setequal(names(observed), names(expected))
  expect_equal(observed[names(expected)], expected, tolerance = 1e-12)

  parameterized <- purrr::map(
    seq_len(nrow(plan)),
    ~ dar:::reparameterize_model_contrast(
      compiled, plan[.x, , drop = FALSE]
    )
  )
  expected_coefficients <- c(
    rep("conditiontreated", 3),
    "time1", "time2", "time2",
    "time1", "time2", "time2",
    "conditiontreated:time1",
    "conditiontreated:time2",
    "conditiontreated:time2"
  )
  expect_equal(
    vapply(parameterized, `[[`, character(1), "coefficient"),
    expected_coefficients
  )
  expect_equal(
    vapply(parameterized, `[[`, numeric(1), "sign"),
    rep(1, nrow(plan))
  )
})

test_that("time factorization is bounded and does not mutate the input", {
  phy <- make_longitudinal_phy()
  original <- as(phyloseq::sample_data(phy), "data.frame")$time
  rec <- recipe(phy) |>
    add_model(
      ~ condition * time, targets = "condition", tax_level = "Species",
      time = "time", max_time_levels = 2
    )

  expect_error(dar:::resolve_model(rec), class = "dar_error_too_many_time_levels")
  expect_equal(as(phyloseq::sample_data(phy), "data.frame")$time, original)
})

test_that("NA handling uses one explicit complete-case cohort", {
  phy <- make_longitudinal_phy()
  meta <- as(phyloseq::sample_data(phy), "data.frame")
  meta$batch[[1]] <- NA
  phyloseq::sample_data(phy) <- phyloseq::sample_data(meta)

  error_rec <- recipe(phy) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time"
    )
  expect_error(dar:::resolve_model(error_rec), class = "dar_error_model_missing_values")

  drop_rec <- recipe(phy) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time",
      na_action = "drop"
    )
  applied <- dar:::apply_model_to_recipe(drop_rec)
  expect_length(applied$resolved$dropped_samples, 1)
  expect_equal(phyloseq::nsamples(applied$rec@phyloseq), 11)
})

test_that("non-estimable fixed designs and empty interaction cells fail", {
  phy <- make_longitudinal_phy()
  meta <- as(phyloseq::sample_data(phy), "data.frame")
  meta$duplicate_batch <- meta$batch
  phyloseq::sample_data(phy) <- phyloseq::sample_data(meta)
  collinear <- recipe(phy) |>
    add_model(
      ~ condition + batch + duplicate_batch,
      targets = "condition", tax_level = "Species"
    )
  expect_error(dar:::resolve_model(collinear), class = "dar_error_non_estimable_model")

  keep <- !(meta$condition == "treated" & meta$time == 2)
  missing_cell_phy <- phyloseq::prune_samples(rownames(meta)[keep], phy)
  missing_cell <- recipe(missing_cell_phy) |>
    add_model(
      ~ condition * time, targets = "condition", tax_level = "Species",
      time = "time"
    )
  expect_error(dar:::resolve_model(missing_cell), class = "dar_error_non_estimable_model")
})

test_that("capability checks skip only unsupported statistical designs", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch + (1 | subject),
      targets = "condition", tax_level = "Species", time = "time"
    )
  expect_warning(
    rec <- dar:::add_step(rec, dar:::step("wilcox", id = "wilcox__model")),
    class = "dar_warning_model_step_skipped"
  )
  rec <- dar:::add_step(rec, dar:::step("ancom", id = "ancom__model"))
  status <- suppressWarnings(dar:::model_steps_status(rec))

  expect_false(status$compatible[status$engine == "wilcox"])
  expect_true(status$compatible[status$engine == "ancom"])
})

test_that("legacy recipes remain model-free", {
  expect_warning(
    rec <- recipe(make_longitudinal_phy(), "condition", "Species"),
    class = "dar_warning_deprecated_recipe_argument"
  )
  expect_null(get_model(rec))
  expect_equal(names(sample_data(rec)), c("sample_id", "condition"))
  expect_true(methods::validObject(rec, test = TRUE))
})

test_that("model-free DA execution emits one classified migration warning", {
  rec <- suppressWarnings(
    recipe(make_longitudinal_phy(), "condition", "Species")
  )
  rec <- dar:::add_step(rec, dar:::step("mock", id = "mock__legacy"))

  expect_warning(
    dar:::warn_model_free_da(rec),
    class = "dar_warning_deprecated_legacy_da"
  )
  expect_silent(dar:::warn_model_free_da(recipe(make_longitudinal_phy())))
})

test_that("models support multiple targets and are independent of step order", {
  rec <- recipe(make_longitudinal_phy())
  before <- add_model(
    rec, ~ condition + batch,
    targets = c("condition", "batch"), tax_level = "Species"
  )
  after <- suppressWarnings(
    dar:::add_step(rec, dar:::step("wilcox", id = "wilcox__ordered")) |>
      add_model(
        ~ condition + batch,
        targets = c("condition", "batch"), tax_level = "Species"
      )
  )

  expect_equal(get_model(before)$targets, c("condition", "batch"))
  expect_equal(
    names(sample_data(before)),
    c("sample_id", "condition", "batch")
  )
  expect_equal(get_model(before)$formula, get_model(after)$formula)
})

test_that("random effects are rechecked after cohort changes", {
  phy <- make_longitudinal_phy()
  metadata <- as(phyloseq::sample_data(phy), "data.frame")
  metadata$unique_subject <- rownames(metadata)
  phyloseq::sample_data(phy) <- phyloseq::sample_data(metadata)
  rec <- recipe(phy) |>
    add_model(
      ~ condition + (1 | unique_subject),
      targets = "condition", tax_level = "Species"
    )

  expect_error(
    dar:::resolve_model(rec),
    class = "dar_error_non_estimable_model"
  )
})

test_that("all engine compilers receive the same design and contrast plan", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time"
    )
  engines <- c("deseq", "aldex", "ancom", "corncob", "maaslin")
  compiled <- lapply(engines, function(engine) dar:::compile_model(rec, engine))

  expect_true(all(vapply(compiled, function(x) {
    identical(colnames(x$matrix), colnames(compiled[[1]]$matrix)) &&
      identical(x$contrasts$contrast_id, compiled[[1]]$contrasts$contrast_id)
  }, logical(1))))
  parameterized <- dar:::reparameterize_model_contrast(
    compiled[[1]], compiled[[1]]$contrasts[5, , drop = FALSE]
  )
  expect_equal(ncol(parameterized$matrix), ncol(compiled[[1]]$matrix))
  expect_true(parameterized$coefficient %in% colnames(parameterized$matrix))
})

test_that("export and import preserve the central model", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time",
      reference = c(condition = "control"), contrasts = "all",
      na_action = "drop"
    )
  path <- tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- import_steps(
    recipe(make_longitudinal_phy()), path
  )

  expect_equal(deparse(get_model(imported)$formula), deparse(get_model(rec)$formula))
  expect_identical(get_model(imported)$reference, get_model(rec)$reference)
  expect_equal(get_model(imported)$contrasts, "all")
  expect_equal(get_model(imported)$na_action, "drop")
  expect_equal(get_model(imported)$targets, "condition")
  expect_equal(get_model(imported)$tax_level, "Species")
})

test_that("model consensus never mixes contrasts or effect directions", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step("wilcox", id = "wilcox__one"))
  rec <- dar:::add_step(rec, dar:::step("lefse", id = "lefse__two"))
  row <- tibble::tibble(
    taxa_id = "taxon_1", taxa = "Species_1",
    contrast_id = "condition[treated-control]",
    comparison = "condition[treated-control]", contrast_type = "main",
    var = "condition", signif = TRUE
  )
  results <- list(
    wilcox__one = list(model = dplyr::mutate(row, effect = 1)),
    lefse__two = list(model = dplyr::mutate(row, effect = -1))
  )
  prepared <- dar:::prep_recipe(
    rec, results, list(),
    execution = list(
      executed_steps = names(results),
      skipped_steps = tibble::tibble(
        step_id = character(), engine = character(),
        compatible = logical(), reason = character()
      )
    )
  )

  expect_equal(nrow(dar:::run_bake(prepared, 2, NULL, NULL, "bake__x")), 0)
  prepared@results$lefse__two$model$effect <- 2
  baked <- dar:::run_bake(prepared, 2, NULL, NULL, "bake__x")
  expect_equal(nrow(baked), 1)
  expect_equal(baked$contrast_id, "condition[treated-control]")
  expect_equal(baked$effect, "up")
  intersections <- find_intersections(prepared)
  expect_equal(nrow(intersections), 1)
  expect_equal(intersections$sum_methods, 2)
})

test_that("modeled exclusion plots preserve contrast and effect keys", {
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step("deseq", id = "deseq__one"))
  rec <- dar:::add_step(rec, dar:::step("maaslin", id = "maaslin__two"))
  contrasts <- dar:::resolve_model(rec)$contrast_plan$contrast_id[1:2]
  rows <- tibble::tibble(
    taxa_id = "taxon_1",
    taxa = "Species_1",
    contrast_id = contrasts,
    comparison = contrasts,
    contrast_type = "main",
    var = "condition",
    effect = c(1, -1),
    padj = 0.01,
    signif = TRUE
  )
  results <- list(
    deseq__one = list(model = rows),
    maaslin__two = list(model = dplyr::slice(rows, 1))
  )
  prepared <- dar:::prep_recipe(rec, results, list())

  plot <- exclusion_plt(prepared)
  expect_s3_class(plot, "ggplot")
  expect_equal(
    dplyr::filter(plot$data, .data$method == "deseq__one") |>
      dplyr::arrange(.data$sum) |>
      dplyr::select("sum", "n", "total"),
    tibble::tibble(sum = c(1L, 2L), n = c(1L, 1L), total = c(2L, 2L))
  )
  expect_equal(
    dplyr::filter(plot$data, .data$method == "maaslin__two") |>
      dplyr::select("sum", "n", "total"),
    tibble::tibble(sum = 2L, n = 1L, total = 1L)
  )

  subset_plot <- exclusion_plt(prepared, steps = "deseq__one")
  expect_equal(unique(subset_plot$data$method), "deseq__one")
  expect_equal(subset_plot$data$sum, 1L)
  expect_equal(subset_plot$data$n, 2L)
})

test_that("modeled mutual plots normalize adjusted p-value aliases", {
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  rec <- dar:::add_step(rec, dar:::step("lefse", id = "lefse__stats"))
  contrast <- dar:::resolve_model(rec)$contrast_plan$contrast_id[[1]]
  rows <- tibble::tibble(
    taxa_id = c("taxon_1", "taxon_2"),
    taxa = c("Species_1", "Species_2"),
    contrast_id = contrast,
    comparison = contrast,
    contrast_type = "main",
    var = "condition",
    effect = c(1, -1),
    adjp = c(0.02, 0.03),
    padj = c(0.01, 0.015),
    signif = TRUE
  )
  prepared <- dar:::prep_recipe(
    rec,
    list(lefse__stats = list(model = rows)),
    list()
  )

  stats <- dar:::.all_stats(prepared)
  expect_named(stats, c("taxa_id", "comparison", "effect_v", "padj", "method"))
  expect_equal(stats$padj, rows$padj)
  expect_s3_class(
    suppressMessages(mutual_plt(prepared, count_cutoff = 1)),
    "ggplot"
  )
})

test_that("model result harmonization is idempotent", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  result <- tibble::tibble(
    taxa_id = "taxon_1",
    contrast_id = "condition[treated-control]",
    comparison = "condition[treated-control]",
    contrast_type = "main",
    var = "condition",
    effect = 1
  )

  once <- dar:::harmonize_model_result(result, rec)
  twice <- dar:::harmonize_model_result(once, rec)
  expect_identical(twice, once)
  expect_false(any(grepl("\\.[xy]$", names(twice))))
})

test_that("prep records skipped steps and errors when none can run", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time + batch,
      targets = "condition", tax_level = "Species", time = "time"
    )
  rec <- suppressWarnings(
    dar:::add_step(rec, dar:::step("wilcox", id = "arbitrary-wilcox-id"))
  )

  expect_error(
    suppressWarnings(prep(rec, parallel = FALSE)),
    class = "dar_error_no_compatible_steps"
  )

  prepared <- dar:::prep_recipe(
    rec, list(), list(),
    execution = list(
      executed_steps = character(),
      skipped_steps = dar:::model_steps_status(rec)
    )
  )
  expect_length(steps_ids(prepared, "da"), 0)
  expect_equal(
    steps_ids(prepared, "da", include_skipped = TRUE),
    "arbitrary-wilcox-id"
  )
})

test_that("sequential and parallel model execution preserve contrast IDs", {
  skip_if_not_installed("microbiome")
  skip_if_not_installed("rstatix")
  skip_if_not_installed("furrr")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_wilcox(id = "wilcox-model")

  sequential <- prep(rec, parallel = FALSE)
  parallel <- prep(rec, parallel = TRUE, workers = 1)
  expect_identical(
    unique(sequential@results[[1]]$model$contrast_id),
    unique(parallel@results[[1]]$model$contrast_id)
  )
})
