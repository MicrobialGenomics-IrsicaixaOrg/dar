rademu_test_fit <- function(rec, ...) {
  prep(step_rademu(rec, id = "rademu__test", ...), parallel = FALSE)
}

rademu_native_fit <- function(rec, metadata, coefficient, test = "score",
                              cluster = NULL, constraint_param = 0.1,
                              reads = NULL, ...) {
  metadata <- data.frame(metadata, row.names = metadata$sample_id, check.names = FALSE)
  if (is.null(reads)) {
    reads <- t(dar:::recipe_count_matrix(rec, phy = dar:::analysis_phyloseq(rec)))
  }
  reads <- reads[rownames(metadata), , drop = FALSE]
  withr::with_seed(
    130,
    radEmu::emuFit(
      Y = reads,
      formula = dar:::resolve_model(rec)$fixed_formula,
      data = metadata,
      cluster = cluster,
      constraint_param = constraint_param,
      test_kj = data.frame(k = coefficient, j = seq_len(ncol(reads))),
      run_score_tests = identical(test, "score"),
      return_wald_p = identical(test, "wald"),
      ...
    ),
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
}

expect_rademu_native <- function(observed, native, coefficient,
                                 p_adj_method = "BH", max_significance = 0.05,
                                 test = "score") {
  expected <- native$coef[native$coef$covariate == coefficient, , drop = FALSE]
  expected <- expected[match(observed$taxa_id, expected$category), , drop = FALSE]
  p_value <- expected[[if (test == "score") "pval" else "wald_p"]]
  adjusted <- p.adjust(p_value, method = p_adj_method)
  expect_equal(observed$effect_size, expected$estimate, tolerance = 1e-8)
  expect_equal(observed$p_value, p_value, tolerance = 1e-8)
  expect_equal(observed$adj_p_value, adjusted, tolerance = 1e-8)
  expect_identical(observed$significant, !is.na(adjusted) & adjusted < max_significance)
  expect_true(all(observed$method == "rademu"))
  expect_true(all(observed$effect_metric == "rademu_centered_log_fold_difference"))
}

test_that("radEmu score results agree with an adjusted native fit", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition + age, targets = "condition", tax_level = "Species")
  metadata <- dar:::resolve_model(rec)$data
  native <- rademu_native_fit(rec, metadata, "conditiontreated")
  withr::local_seed(909)
  original_rng <- .Random.seed
  prepared <- rademu_test_fit(rec)
  expect_identical(.Random.seed, original_rng)
  expect_rademu_native(tidy_results(prepared), native, "conditiontreated")
  details <- get_execution(prepared)$engine_details$rademu__test
  expect_identical(details$package, "radEmu")
  expect_identical(details$version, as.character(getNamespaceVersion("radEmu")))
  expect_identical(details$test, "score")
  expect_identical(details$constraint, "pseudohuber_median")
  expect_identical(details$constraint_param, 0.1)
  expect_identical(details$seed, 130L)
  expect_null(details$cluster)
  expect_equal(details$native_calls, 1L)
  expect_true(details$elapsed_seconds >= 0)
  diagnostic <- details$diagnostics[[1]]
  expect_equal(diagnostic$score_test_hyperparams, native$score_test_hyperparams)
  expect_identical(diagnostic$estimation_converged, native$estimation_converged)
  expect_identical(tidy_results(rademu_test_fit(rec)), tidy_results(prepared))
})

test_that("radEmu Wald tests and factor references retain native inference", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(
      ~ condition, targets = "condition", tax_level = "Species",
      reference = c(condition = "treated")
    )
  prepared <- rademu_test_fit(
    rec, test = "wald", p_adj_method = "holm", max_significance = 0.1,
    constraint_param = 0.2, engine_args = list(fit = list(tolerance = 1e-5))
  )
  observed <- tidy_results(prepared)
  plan <- dar:::resolve_model(rec)$contrast_plan
  for (index in seq_len(nrow(plan))) {
    contrast <- plan[index, ]
    metadata <- dar:::resolve_model(rec)$data
    metadata$condition <- relevel(metadata$condition, contrast$denominator)
    coefficient <- paste0("condition", contrast$numerator)
    native <- rademu_native_fit(
      rec, metadata, coefficient, test = "wald", constraint_param = 0.2,
      tolerance = 1e-5
    )
    rows <- observed[observed$contrast_id == contrast$contrast_id, ]
    expect_rademu_native(rows, native, coefficient, "holm", 0.1, "wald")
    raw <- prepared@results$rademu__test$model
    raw <- raw[raw$contrast_id == contrast$contrast_id, ]
    expect_equal(raw$pvalue, raw$wald_p)
    expect_true(all(is.na(raw$pval)))
  }
})

test_that("radEmu clustered interactions are explicit native GEE fits", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time, targets = "condition", time = "time",
      tax_level = "Species", contrasts = "all"
    )
  prepared <- rademu_test_fit(rec, cluster = "subject", test = "wald")
  observed <- tidy_results(prepared)
  resolved <- dar:::resolve_model(rec)
  expect_setequal(unique(observed$contrast_id), resolved$contrast_plan$contrast_id)
  native <- rademu_native_fit(
    rec, resolved$data, "conditiontreated:time1", test = "wald",
    cluster = resolved$data$subject
  )
  contrast_id <- resolved$contrast_plan$contrast_id[
    resolved$contrast_plan$contrast_type == "difference_in_differences" &
      resolved$contrast_plan$numerator == "treated:1"
  ]
  expect_length(contrast_id, 1L)
  rows <- observed[observed$contrast_id == contrast_id, ]
  expect_rademu_native(rows, native, "conditiontreated:time1", test = "wald")
  details <- get_execution(prepared)$engine_details$rademu__test
  expect_identical(details$cluster, "subject")
  expect_identical(details$sample_ids, resolved$data$sample_id)
  expect_identical(details$requested_cluster_ids, resolved$data$subject)
  expect_true(all(vapply(details$diagnostics, function(item) {
    identical(item$cluster_ids, native$cluster)
  }, logical(1))))
})

test_that("radEmu score tests preserve native clustered covariance", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  resolved <- dar:::resolve_model(rec)
  native <- rademu_native_fit(
    rec, resolved$data, "conditiontreated", cluster = resolved$data$subject
  )
  observed <- tidy_results(rademu_test_fit(rec, cluster = "subject"))
  expect_rademu_native(observed, native, "conditiontreated")
})

test_that("radEmu score difference-in-differences matches its native hypothesis", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  phy <- make_longitudinal_phy()
  phy <- phyloseq::prune_samples(phyloseq::get_variable(phy, "time") != 2, phy)
  rec <- recipe(phy) |>
    add_model(
      ~ condition * time, targets = "condition", time = "time",
      tax_level = "Species", contrasts = "all"
    )
  prepared <- rademu_test_fit(rec, cluster = "subject")
  observed <- tidy_results(prepared)
  resolved <- dar:::resolve_model(rec)
  native <- rademu_native_fit(
    rec, resolved$data, "conditiontreated:time1", cluster = resolved$data$subject
  )
  rows <- observed[observed$contrast_type == "difference_in_differences", ]
  expect_rademu_native(rows, native, "conditiontreated:time1")
  expect_equal(nrow(observed), 5L * phyloseq::ntaxa(phy))
})

test_that("radEmu does not replace random effects with clustering", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition + (1 | subject), targets = "condition", tax_level = "Species"
    )
  expect_warning(
    configured <- step_rademu(rec, cluster = "subject"),
    "does not support random effects", class = "dar_warning_model_step_skipped"
  )
  expect_false(dar:::model_steps_status(configured)$compatible)
  expect_warning(
    expect_error(
      prep(configured, parallel = FALSE),
      class = "dar_error_no_compatible_steps"
    ),
    class = "dar_warning_model_step_skipped"
  )
})

test_that("radEmu uses the common cohort and aligns cluster IDs by sample", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  phy <- make_longitudinal_phy()
  phyloseq::sample_data(phy)$age[2] <- NA_real_
  rec <- recipe(phy) |>
    add_model(
      ~ condition + age, targets = "condition", tax_level = "Species",
      na_action = "drop"
    )
  prepared <- suppressWarnings(rademu_test_fit(rec, cluster = "subject", test = "wald"))
  resolved <- dar:::resolve_model(rec)
  native <- rademu_native_fit(
    rec, resolved$data, "conditiontreated", test = "wald",
    cluster = resolved$data$subject
  )
  expect_rademu_native(tidy_results(prepared), native, "conditiontreated", test = "wald")
  details <- get_execution(prepared)$engine_details$rademu__test
  expect_identical(details$sample_ids, resolved$data$sample_id)
  expect_identical(details$requested_cluster_ids, resolved$data$subject)
  expect_identical(
    get_execution(prepared)$dropped_samples,
    phyloseq::sample_names(phy)[2]
  )
})

test_that("radEmu preserves sparse native zero-comparison decisions", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  phy <- make_multilevel_phy()
  counts <- as(phyloseq::otu_table(phy), "matrix")
  group <- phyloseq::get_variable(phy, "condition")
  counts[1, group != "treated"] <- 0
  phyloseq::otu_table(phy) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  rec <- recipe(phy) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  prepared <- rademu_test_fit(
    rec, engine_args = list(fit = list(remove_zero_comparison_pvals = TRUE))
  )
  observed <- tidy_results(prepared)
  metadata <- dar:::resolve_model(rec)$data
  native <- rademu_native_fit(
    rec, metadata, "conditionrescue", remove_zero_comparison_pvals = TRUE
  )
  rows <- observed[observed$contrast_id == "condition[rescue-control]", ]
  expect_rademu_native(rows, native, "conditionrescue")
  expect_true(is.na(rows$p_value[rows$taxa_id == "taxon_1"]))
  expect_false(rows$significant[rows$taxa_id == "taxon_1"])
  expect_identical(get_phy(rec), phy)
  expect_identical(get_phy(prepared), phy)
})

test_that("radEmu fits aggregated counts and keeps homonymous lineages distinct", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  phy <- make_multilevel_phy()
  taxonomy <- cbind(
    Family = c("first", "first", "second", "third", "fourth", "fifth"),
    Genus = c("Shared", "Shared", "Shared", "Other", "Another", NA_character_),
    Species = paste0("Species_", seq_len(6))
  )
  rownames(taxonomy) <- phyloseq::taxa_names(phy)
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  original <- as(phyloseq::otu_table(phy), "matrix")
  aggregated <- rbind(
    taxon_1 = original[1, ] + original[2, ],
    original[3:6, , drop = FALSE]
  )
  previous <- NULL
  for (taxa_are_rows in c(TRUE, FALSE)) {
    if (!taxa_are_rows) {
      phyloseq::otu_table(phy) <- phyloseq::otu_table(
        t(as(phyloseq::otu_table(phy), "matrix")), taxa_are_rows = FALSE
      )
    }
    rec <- recipe(phy) |>
      add_model(~ condition, targets = "condition", tax_level = "Genus")
    prepared <- rademu_test_fit(rec, test = "wald")
    observed <- tidy_results(prepared)
    native <- rademu_native_fit(
      rec, dar:::resolve_model(rec)$data, "conditionrescue", test = "wald",
      reads = t(aggregated)
    )
    rows <- observed[observed$contrast_id == "condition[rescue-control]", ]
    expect_rademu_native(rows, native, "conditionrescue", test = "wald")
    expect_equal(length(unique(observed$taxa_id)), 5L)
    expect_false("taxon_2" %in% observed$taxa_id)
    expect_equal(length(unique(observed$taxa_id[observed$taxa %in% "Shared"])), 2L)
    if (!is.null(previous)) { expect_identical(observed, previous) }
    expect_identical(get_phy(prepared), phy)
    previous <- observed
  }
})

test_that("radEmu missing effects and P values remain missing", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  native_stage <- dar:::exec_engine_stage
  local_mocked_bindings(
    exec_engine_stage = function(...) {
      fit <- native_stage(...)
      fit$coef$estimate[fit$coef$category == "taxon_1"] <- NA_real_
      fit$coef$wald_p[fit$coef$category == "taxon_1"] <- NA_real_
      fit
    },
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  observed <- tidy_results(rademu_test_fit(rec, test = "wald"))
  missing <- observed$taxa_id == "taxon_1"
  expect_true(all(is.na(observed$effect_size[missing])))
  expect_true(all(is.na(observed$p_value[missing])))
  expect_true(all(is.na(observed$adj_p_value[missing])))
  expect_false(any(observed$significant[missing]))
})

test_that("radEmu receives one complete rarefied matrix and supports exploration", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  skip_if_not_installed("microbiome")
  native_stage <- dar:::exec_engine_stage
  received <- list()
  local_mocked_bindings(
    exec_engine_stage = function(engine, stage, engine_args, fixed, defaults = list()) {
      received[[length(received) + 1L]] <<- fixed$Y
      native_stage(engine, stage, engine_args, fixed, defaults)
    },
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  original <- get_phy(rec)
  withr::local_seed(130)
  expect_message(
    prepared <- rademu_test_fit(rec, test = "wald", rarefy = TRUE),
    "rarefaction"
  )
  minimum_depth <- min(colSums(dar:::recipe_count_matrix(rec)))
  expect_length(received, 3L)
  for (reads in received) {
    expect_identical(rownames(reads), phyloseq::sample_names(original))
    expect_identical(colnames(reads), phyloseq::taxa_names(original))
    expect_true(all(rowSums(reads) == minimum_depth))
    expect_identical(reads, received[[1]])
  }
  expect_identical(get_phy(prepared), original)
  rows <- tidy_results(prepared)
  contrast <- rows$contrast_id[[1]]
  plot <- effect_plt(prepared, contrast_id = contrast, top_n = 3)
  expect_s3_class(plot, "ggplot")
  expect_true(all(plot$data$method == "rademu"))
  expect_equal(nrow(plot$data), 3L)
  abundance <- tidy_abundance(prepared, taxa_ids = rows$taxa_id[[1]])
  expect_equal(nrow(abundance), phyloseq::nsamples(original))
  expect_setequal(intersection_df(prepared)$taxa_id, rows$taxa_id)
})

test_that("radEmu checks integration settings but leaves native parameters native", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  expect_error(step_rademu(rec, test = "unknown"))
  expect_error(step_rademu(rec, cluster = c("first", "second")))
  expect_error(
    rademu_test_fit(rec, cluster = "missing"),
    class = "dar_error_missing_vars"
  )
  phyloseq::sample_data(rec@phyloseq)$cluster <- c(NA, rep("one", 11))
  expect_error(
    rademu_test_fit(rec, cluster = "cluster"), class = "dar_error_model_missing_values"
  )
  protected <- dar:::engine_args_registry()$rademu$stages$fit$protected
  for (argument in protected) {
    expect_error(
      step_rademu(rec, engine_args = list(fit = setNames(list(TRUE), argument))),
      class = "dar_error_invalid_engine_args"
    )
  }
  native_error <- structure(list(message = "native parameter error"),
                            class = c("native_error", "error", "condition"))
  received <- NULL
  local_mocked_bindings(
    emuFit = function(...) {
      received <<- list(...)
      stop(native_error)
    },
    .package = "radEmu"
  )
  withr::local_seed(809)
  rng_before <- .Random.seed
  error <- tryCatch(
    rademu_test_fit(
      rec, constraint_param = -1, engine_args = list(fit = list(maxit = -10))
    ),
    error = identity
  )
  while (!is.null(error$parent)) { error <- error$parent }
  expect_identical(error, native_error)
  expect_identical(received$constraint_param, -1)
  expect_identical(received$maxit, -10)
  expect_identical(.Random.seed, rng_before)
})

test_that("radEmu output identities are part of the result contract", {
  skip_if_not_installed("radEmu", "2.3.2.0")
  native_stage <- dar:::exec_engine_stage
  local_mocked_bindings(
    exec_engine_stage = function(...) {
      fit <- native_stage(...)
      fit$coef$category <- "unrelated_taxon"
      fit
    },
    .package = "dar"
  )
  rec <- recipe(make_multilevel_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  expect_error(
    rademu_test_fit(rec, test = "wald"),
    class = "dar_error_invalid_result_contract"
  )
})

test_that("radEmu configuration can be exported without executing models", {
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(
      ~ condition * time, targets = "condition", time = "time", tax_level = "Species"
    ) |>
    step_rademu(
      cluster = "subject", constraint_param = 0.2,
      engine_args = list(fit = list(tolerance = 1e-5))
    )
  path <- withr::local_tempfile(fileext = ".json")
  export_steps(rec, path)
  imported <- import_steps(recipe(make_longitudinal_phy()), path)
  expect_identical(imported@steps, rec@steps)
  expect_identical(
    dar:::recipe_serialization_document(imported),
    dar:::recipe_serialization_document(rec)
  )
})
