aldex3_test_fit <- function(rec, ...) {
  prep(step_aldex3(rec, nsample = 8, id = "aldex3__test", ...), parallel = FALSE)
}

expect_aldex3_native <- function(observed, native, coefficient, direction = 1) {
  expected <- summary(native, ignore.intercept = FALSE)
  expected <- expected[expected$parameter == coefficient, , drop = FALSE]
  observed <- observed[match(expected$entity, observed$taxa_id), , drop = FALSE]
  expect_identical(observed$taxa_id, expected$entity)
  expect_equal(observed$effect_size, expected$estimate * direction, tolerance = 1e-8)
  expect_equal(observed$p_value, as.numeric(native$p.val[coefficient, ]))
  expect_equal(observed$adj_p_value, expected$p.val.adj)
  expect_identical(observed$significant, !is.na(expected$p.val.adj) & expected$p.val.adj < 0.05)
  expect_true(all(observed$method == "aldex3"))
  expect_true(all(observed$effect_metric == "aldex3_log2_abundance_coefficient"))
}

test_that("ALDEx2 has an explicit constructor and no ambiguous alias", {
  expect_true("step_aldex2" %in% getNamespaceExports("dar"))
  expect_false("step_aldex" %in% getNamespaceExports("dar"))
  expect_false(exists("step_aldex", envir = asNamespace("dar"), inherits = FALSE))
  expect_false("step_aldex" %in% names(dar:::recipe_step_registry))
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_aldex2(id = "aldex2__explicit")
  expect_identical(dar:::step_engine(rec@steps[[1]]), "aldex2")
  path <- withr::local_tempfile(fileext = ".json")
  export_steps(rec, path)
  document <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  document$steps[[1]]$constructor <- "step_aldex"
  jsonlite::write_json(document, path, auto_unbox = TRUE, null = "null")
  expect_error(import_steps(recipe(make_multilevel_phy()), path), "not supported",
               class = "dar_error_invalid_recipe_record")
})

test_that("ALDEx3 agrees with native adjusted fixed-model inference", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species")
  compiled <- dar:::compile_model(rec, "aldex3")
  reads <- dar:::recipe_count_matrix(rec)
  native <- withr::with_seed(130, ALDEx3::aldex(
    reads, ~ condition + age, compiled$data, nsample = 8,
    scale = ALDEx3::clr.sm, gamma = 0.5, n.cores = 1,
    return.pars = c("estimate", "std.error", "p.val", "p.val.adj")
  ))
  withr::local_seed(909)
  rng_before <- .Random.seed
  prepared <- aldex3_test_fit(rec)
  expect_identical(.Random.seed, rng_before)
  expect_aldex3_native(tidy_results(prepared), native, "conditiontreated")
  details <- get_execution(prepared)$engine_details$aldex3__test
  expect_identical(details$version, as.character(getNamespaceVersion("ALDEx3")))
  expect_identical(details$method, "lm")
  expect_identical(details$seed, 130L)
  expect_equal(details$nsample, 8)
  expect_identical(details$scale_args, list(gamma = 0.5))
  expect_identical(details$scale_basis, colnames(compiled$matrix))
  expect_true(is.finite(details$elapsed_seconds) && details$elapsed_seconds >= 0)
  expect_identical(tidy_results(aldex3_test_fit(rec)), tidy_results(prepared))
  withr::local_rng_version(getRversion())
  RNGkind("L'Ecuyer-CMRG")
  expect_identical(tidy_results(aldex3_test_fit(rec)), tidy_results(prepared))
})

test_that("contrast releveling does not redefine coefficient-scale uncertainty", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  compiled <- dar:::compile_model(rec, "aldex3")
  original_X <- t(compiled$matrix)
  reads <- dar:::recipe_count_matrix(rec)
  prepared <- aldex3_test_fit(rec, scale_args = list(gamma = 0.7))
  observed <- tidy_results(prepared)
  native_scale <- function(X, logComp) {
    ALDEx3::clr.sm(original_X, logComp, gamma = 0.7)
  }
  for (index in seq_len(nrow(compiled$contrasts))) {
    contrast <- compiled$contrasts[index, ]
    metadata <- compiled$data
    metadata$condition <- relevel(metadata$condition, contrast$denominator)
    coefficient <- paste0("condition", contrast$numerator)
    native <- withr::with_seed(130, ALDEx3::aldex(
      reads, ~ condition, metadata, nsample = 8, scale = native_scale, n.cores = 1
    ))
    rows <- observed[observed$contrast_id == contrast$contrast_id, ]
    expect_aldex3_native(rows, native, coefficient)
    # Check the estimand independently using the original coefficient basis.
    base_fit <- withr::with_seed(130, ALDEx3::aldex(
      reads, ~ condition, compiled$data, nsample = 8,
      scale = ALDEx3::clr.sm, gamma = 0.7, n.cores = 1
    ))
    coefficients <- rowMeans(base_fit$estimate, dims = 2)
    expected <- as.numeric(crossprod(contrast$weights[[1]], coefficients))
    expect_equal(rows$effect_size, expected, tolerance = 1e-8)
  }
})

test_that("ALDEx3 handles exact mixed models and interaction contrasts", {
  skip_if_not_installed("ALDEx3")
  phy <- make_longitudinal_phy()
  phy <- phyloseq::prune_samples(phyloseq::get_variable(phy, "time") != 2, phy)
  counts <- as(phyloseq::otu_table(phy), "matrix")
  counts[1, phyloseq::get_variable(phy, "condition") == "control"] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  rec <- recipe(phy) |>
    add_model(
      ~ condition * time + (1 | subject), targets = "condition",
      time = "time", tax_level = "Species", contrasts = "all"
    )
  compiled <- dar:::compile_model(rec, "aldex3")
  prepared <- suppressWarnings(aldex3_test_fit(rec, scale_args = list(gamma = 0)))
  observed <- tidy_results(prepared)
  expect_setequal(unique(observed$contrast_id), compiled$contrasts$contrast_id)
  expect_false(anyDuplicated(observed[c("taxa_id", "contrast_id", "step_id")]) > 0)
  native <- suppressWarnings(withr::with_seed(130, ALDEx3::aldex(
    counts, ~ condition * time + (1 | subject), compiled$data,
    nsample = 8, method = "lme4", scale = ALDEx3::clr.sm, gamma = 0, n.cores = 1
  )))
  did <- compiled$contrasts$contrast_type == "difference_in_differences"
  expect_equal(sum(did), 1L)
  rows <- observed[observed$contrast_id == compiled$contrasts$contrast_id[did], ]
  expect_aldex3_native(rows, native, "conditiontreated:time1")
  baseline <- observed[observed$contrast_id == "condition[treated-control]@time[0]", ]
  expect_aldex3_native(baseline, native, "conditiontreated")
  expect_identical(get_execution(prepared)$engine_details$aldex3__test$method, "lme4")
  expect_true(all(is.finite(observed$effect_size)))
})

test_that("external scale vectors are aligned by sample ID after exclusions", {
  skip_if_not_installed("ALDEx3")
  phy <- make_longitudinal_phy()
  phyloseq::sample_data(phy)$age[2] <- NA_real_
  ids <- rev(phyloseq::sample_names(phy))
  args <- list(s.mu = setNames(seq_along(ids) / 10, ids),
               s.sd = setNames(rep(0.2, length(ids)), ids))
  rec <- recipe(phy) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species",
              na_action = "drop")
  prepared <- suppressWarnings(aldex3_test_fit(rec, scale_model = "sample", scale_args = args))
  compiled <- dar:::compile_model(prepared, "aldex3")
  reads <- dar:::recipe_count_matrix(prepared)
  native <- withr::with_seed(130, ALDEx3::aldex(
    reads, ~ condition + age, compiled$data, nsample = 8,
    scale = ALDEx3::sample.sm, s.mu = unname(args$s.mu[colnames(reads)]),
    s.sd = unname(args$s.sd[colnames(reads)]), n.cores = 1
  ))
  expect_aldex3_native(tidy_results(prepared), native, "conditiontreated")
  expect_identical(get_execution(prepared)$dropped_samples, phyloseq::sample_names(phy)[2])
  path <- withr::local_tempfile(fileext = ".json")
  configured <- step_aldex3(rec, nsample = 8, scale_model = "sample", scale_args = args)
  export_steps(configured, path)
  imported <- import_steps(recipe(phy), path)
  expect_identical(imported@steps, configured@steps)
})

test_that("TSS and coefficient-scale models retain their native assumptions", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  compiled <- dar:::compile_model(rec, "aldex3")
  ids <- colnames(compiled$matrix)
  args <- list(c.mu = setNames(c(0, 0.3, -0.2), ids),
               c.sd = setNames(c(0.1, 0.2, 0.4), ids))
  for (model in c("tss", "coefficient")) {
    scale_args <- if (model == "tss") list(gamma = 0.4) else args
    fitted <- aldex3_test_fit(rec, scale_model = model, scale_args = scale_args)
    native_args <- lapply(scale_args, unname)
    native <- withr::with_seed(130, do.call(ALDEx3::aldex, c(list(
      Y = dar:::recipe_count_matrix(rec), X = ~ condition, data = compiled$data,
      nsample = 8, n.cores = 1,
      scale = getExportedValue("ALDEx3", paste0(model, ".sm"))
    ), native_args)))
    observed <- tidy_results(fitted)
    expect_aldex3_native(
      observed[observed$contrast_id == "condition[treated-control]", ],
      native, "conditiontreated"
    )
    path <- withr::local_tempfile(fileext = ".json")
    export_steps(fitted, path)
    imported <- import_steps(recipe(make_multilevel_phy()), path)
    expect_identical(imported@steps, fitted@steps)
  }
})

test_that("ALDEx3 missing statistics are retained, not converted to zero", {
  skip_if_not_installed("ALDEx3")
  native_stage <- dar:::exec_engine_stage
  local_mocked_bindings(exec_engine_stage = function(...) {
    fit <- native_stage(...)
    fit$estimate[, 1, ] <- NA_real_
    fit$p.val[, 1] <- NA_real_
    fit$p.val.adj[, 1] <- NA_real_
    fit
  }, .package = "dar")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  observed <- tidy_results(aldex3_test_fit(rec))
  rows <- observed$taxa_id == "taxon_1"
  expect_true(all(is.na(observed$effect_size[rows])))
  expect_true(all(is.na(observed$p_value[rows])))
  expect_true(all(is.na(observed$adj_p_value[rows])))
  expect_false(any(observed$significant[rows]))
  expect_true(all(is.finite(observed$effect_size[!rows])))
})

test_that("ALDEx3 uses aggregated lineage identities in either matrix orientation", {
  skip_if_not_installed("ALDEx3")
  phy <- make_multilevel_phy()
  taxonomy <- cbind(
    Family = c("first", "first", "second", "third", "fourth", "fifth"),
    Genus = c("Shared", "Shared", "Shared", "Other", "Another", NA_character_),
    Species = paste0("Species_", seq_len(6))
  )
  rownames(taxonomy) <- phyloseq::taxa_names(phy)
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  previous <- NULL
  for (taxa_are_rows in c(TRUE, FALSE)) {
    if (!taxa_are_rows) {
      phyloseq::otu_table(phy) <- phyloseq::otu_table(
        t(as(phyloseq::otu_table(phy), "matrix")), taxa_are_rows = FALSE
      )
    }
    rec <- recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Genus")
    observed <- tidy_results(aldex3_test_fit(rec))
    abundances <- tidy_abundance(rec)
    expect_equal(length(unique(observed$taxa_id)), 5L)
    expect_false("taxon_2" %in% observed$taxa_id)
    expect_equal(length(unique(observed$taxa_id[observed$taxa %in% "Shared"])), 2L)
    expect_setequal(observed$taxa_id, abundances$taxa_id)
    if (!is.null(previous)) expect_identical(observed, previous)
    previous <- observed
  }
})

test_that("ALDEx3 exposes native fixed-model tests and adjustment without post-hoc BH", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species")
  compiled <- dar:::compile_model(rec, "aldex3")
  prepared <- aldex3_test_fit(
    rec, p_adj_method = "holm", scale_args = list(gamma = 0),
    engine_args = list(fit = list(test = "t.HC0", streamsize = 1))
  )
  native <- withr::with_seed(130, ALDEx3::aldex(
    dar:::recipe_count_matrix(rec), ~ condition + age, compiled$data,
    nsample = 8, scale = ALDEx3::clr.sm, gamma = 0, n.cores = 1,
    p.adjust.method = "holm", test = "t.HC0", streamsize = 1
  ))
  expect_aldex3_native(tidy_results(prepared), native, "conditiontreated")
  details <- get_execution(prepared)$engine_details$aldex3__test
  expect_identical(details$p_adj_method, "holm")
  expect_identical(details$test, "t.HC0")
})

test_that("ALDEx3 checks dar configuration and model compatibility", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  expect_error(step_aldex3(rec, scale_args = list(typo = 0.5)))
  expect_error(step_aldex3(rec, scale_model = "unknown"))
  expect_error(
    aldex3_test_fit(
      rec, scale_model = "coefficient", scale_args = list(c.mu = 1:3, c.sd = 1:3)
    ),
    class = "dar_error_invalid_model"
  )
  for (argument in c("method", "scale", "onesided", "nsample", "random", "typo")) {
    expect_error(step_aldex3(rec, engine_args = list(fit = setNames(list(TRUE), argument))))
  }
  expect_error(suppressWarnings(aldex3_test_fit(recipe(make_multilevel_phy()))), "requires")
  phy <- make_multilevel_phy()
  phyloseq::sample_data(phy)$duplicate <- phyloseq::get_variable(phy, "condition")
  singular <- recipe(phy) |>
    add_model(~ condition + duplicate, targets = "condition", tax_level = "Species")
  expect_error(aldex3_test_fit(singular), class = "dar_error_non_estimable_model")
  mixed <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + (1 | subject), targets = "condition", tax_level = "Species"
    )
  expect_error(
    aldex3_test_fit(mixed, engine_args = list(fit = list(test = "t"))),
    class = "dar_error_invalid_engine_args"
  )
})

test_that("ALDEx3 receives native settings and input without dar validation", {
  skip_if_not_installed("ALDEx3")
  observed <- NULL
  native_error <- errorCondition("ALDEx3 rejected a setting", class = "native_error")
  local_mocked_bindings(
    aldex = function(...) {
      observed <<- list(...)
      stop(native_error)
    },
    .package = "ALDEx3"
  )
  phy <- make_multilevel_phy()
  counts <- as(phyloseq::otu_table(phy), "matrix")
  counts[1, 1] <- 0.5
  counts[, 2] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species") |>
    step_aldex3(
      nsample = 0,
      n_cores = 0,
      p_adj_method = "native-choice",
      scale_args = list(gamma = -1),
      engine_args = list(fit = list(test = "native-test", streamsize = Inf))
    )
  error <- tryCatch(prep(rec, parallel = FALSE), error = identity)
  while (!is.null(error$parent)) error <- error$parent
  expect_identical(error, native_error)
  expect_equal(observed$Y, counts)
  expect_identical(observed$nsample, 0)
  expect_identical(observed$n.cores, 0)
  expect_identical(observed$p.adjust.method, "native-choice")
  expect_identical(observed$test, "native-test")
  expect_identical(observed$streamsize, Inf)
  expect_identical(environment(observed$scale)$native_scale_args$gamma, -1)
  expect_s4_class(step_aldex3(recipe(phy), nsample = 1, seed = NA_integer_), "Recipe")
})

test_that("ALDEx3 reports its own invalid numerical scale settings", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  withr::local_seed(42)
  original_seed <- .Random.seed
  error <- tryCatch(
    aldex3_test_fit(rec, scale_args = list(gamma = -1)), error = identity
  )
  while (!is.null(error$parent)) error <- error$parent
  expect_match(conditionMessage(error), "gamma.*non-negative")
  expect_false(inherits(error, "dar_error_invalid_engine_args"))
  expect_identical(.Random.seed, original_seed)
  expect_error(
    aldex3_test_fit(rec, scale_model = "sample", scale_args = list()),
    "s.mu cannot be NULL"
  )
  ids <- phyloseq::sample_names(get_phy(rec))
  expect_error(
    aldex3_test_fit(
      rec,
      scale_model = "sample",
      scale_args = list(
        s.mu = setNames(rep(0, length(ids)), ids),
        s.sd = setNames(rep(-1, length(ids)), ids)
      )
    ),
    "s.sd.*non-negative"
  )
})

test_that("ALDEx3 scale identifiers remain a dar integration requirement", {
  skip_if_not_installed("ALDEx3")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  ids <- phyloseq::sample_names(get_phy(rec))
  means <- setNames(rep(0, length(ids)), ids)
  for (value_ids in list(NULL, rep(ids[1], length(ids)), rev(ids)[-1])) {
    invalid <- means
    if (length(value_ids) && length(value_ids) < length(invalid)) invalid <- invalid[-1]
    names(invalid) <- value_ids
    expect_error(
      aldex3_test_fit(
        rec, scale_model = "sample", scale_args = list(s.mu = invalid)
      ),
      class = "dar_error_invalid_model"
    )
  }
})

test_that("ALDEx3 output identities must match the dar analysis units", {
  local_mocked_bindings(
    exec_engine_stage = function(...) {
      list(p.val.adj = matrix(
        0.1,
        nrow = 2,
        ncol = 6,
        dimnames = list(
          c("(Intercept)", "conditiontreated"), paste0("unexpected_", 1:6)
        )
      ))
    },
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  expect_error(aldex3_test_fit(rec), class = "dar_error_invalid_result_contract")
})
