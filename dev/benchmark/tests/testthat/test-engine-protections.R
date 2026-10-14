# These fixtures inspect controls and returned matrices, never statistical fits.
protection_simulation_fixture <- function() {
  list(
    specification = list(mode = "pilot"),
    contrast_plan = data.frame(numerator = "treated", denominator = "control")
  )
}

protection_ancom_fixture <- function() {
  raw <- data.frame(
    taxon = c("a", "b", "c", "d"),
    diff_conditiontreated = c(TRUE, TRUE, FALSE, NA),
    passed_ss_conditiontreated = c(FALSE, TRUE, TRUE, TRUE),
    diff_robust_conditiontreated = c(FALSE, TRUE, FALSE, NA),
    q_conditiontreated = c(0.001, 0.01, 0.5, 0.001)
  )
  list(
    stages = list(list(stage = "fit", completed = TRUE,
                       arguments = list(pseudo_sens = TRUE), result = list(res = raw))),
    results = data.frame(taxa_id = raw$taxon, significant = raw$diff_robust_conditiontreated)
  )
}

protection_linda_fixture <- function(mode = "imputation", zeros = TRUE) {
  original <- rbind(
    a = c(if (zeros) 0 else 1, 2, 8, 16),
    b = c(6, 12, 24, 48), c = c(3, 6, 12, 24), d = c(2, 4, 8, 16)
  )
  colnames(original) <- paste0("s", seq_len(ncol(original)))
  before <- getFromNamespace("winsor.fun", "MicrobiomeStat")(original, 0.97, "count")
  observed <- before
  if (zeros && mode == "imputation") {
    depth <- colSums(before)
    zero_depth <- matrix(rep(depth, nrow(before)), nrow = nrow(before), byrow = TRUE)
    zero_depth[before > 0] <- 0
    observed <- before + zero_depth / depth[max.col(zero_depth, ties.method = "first")]
  } else if (zeros && mode == "pseudo_count") {
    observed <- before + 0.5
  }
  list(
    stage = "fit", completed = TRUE,
    arguments = list(feature.dat = original, feature.dat.type = "count",
                     is.winsor = TRUE, outlier.pct = 0.03, adaptive = TRUE, pseudo.cnt = 0.5),
    effective_function = benchmark_engine_implementation("linda", "fit"),
    result = list(feature.dat.use = observed)
  )
}

test_that("the pilot explicitly requests robust ANCOM without fitting it", {
  simulation <- benchmark_simulate("PF_G1", preflight = TRUE)
  rec <- benchmark_add_engine(simulation$recipe, "ancom_pilot")
  expect_true(rec@steps[[1L]][["engine_args"]]$fit$pseudo_sens)
  expect_equal(length(rec@steps), 1L)
})

test_that("ANCOM audit preserves robust native decisions including missing values", {
  fixture <- protection_ancom_fixture()
  checks <- benchmark_native_checks("ancom_pilot", fixture$stages, fixture$results,
                                     protection_simulation_fixture())
  expect_true(all(checks$passed))
  expect_true(all(c("ancom_robust_calls", "ancom_robust_conjunction",
                    "ancom_native_significance") %in% checks$check))
  # Small q-values alone do not admit a failed sensitivity flag or an NA decision.
  expect_identical(fixture$results$significant, c(FALSE, TRUE, FALSE, NA))
  fixture$results$significant[1L] <- TRUE
  checks <- benchmark_native_checks("ancom_pilot", fixture$stages, fixture$results,
                                     protection_simulation_fixture())
  expect_false(checks$passed[checks$check == "ancom_native_significance"])
})

test_that("ANCOM audit rejects basic mode and inconsistent or absent sensitivity decisions", {
  fixture <- protection_ancom_fixture()
  fixture$stages[[1L]]$arguments$pseudo_sens <- FALSE
  fixture$stages[[1L]]$result$res$diff_robust_conditiontreated[1L] <- TRUE
  checks <- benchmark_native_checks("ancom_pilot", fixture$stages, fixture$results,
                                     protection_simulation_fixture())
  expect_false(checks$passed[checks$check == "ancom_robust_calls"])
  expect_false(checks$passed[checks$check == "ancom_robust_conjunction"])
  fixture$stages[[1L]]$result$res$passed_ss_conditiontreated <- NULL
  expect_error(benchmark_native_checks("ancom_pilot", fixture$stages, fixture$results,
                                        protection_simulation_fixture()), "sensitivity flag")
})

test_that("LinDA implementation metadata verifies exactly the local repair", {
  native <- MicrobiomeStat::linda
  original_body <- body(native)
  effective <- dar:::engine_stage_function("linda", "fit")
  identity <- benchmark_engine_implementation("linda", "fit", effective)
  expect_identical(identity$implementation, "dar_local_adaptive_zero_fix")
  expect_equal(identity$linda_correction_count, 1)
  expect_true(identity$formals_unchanged)
  expect_true(identity$namespace_unchanged)
  expect_true(identity$linda_adaptive_corrected)
  expect_false(identity$linda_adaptive_bug_present)
  expect_false(identical(identity$effective_sha256, identity$native_sha256))
  expect_identical(body(MicrobiomeStat::linda), original_body)
  old <- benchmark_engine_implementation("linda", "fit", native)
  expect_identical(old$implementation, "native")
  expect_true(old$linda_adaptive_bug_present)
  expect_false(old$linda_adaptive_corrected)
})

test_that("LinDA AST audit detects unrelated edits including structural NULL deletion", {
  native <- function(x, y) {
    value <- list(NULL, x[, y])
    zero.handling <- "Imputation"
    if (zero.handling == "imputation") value
  }
  effective <- dar:::linda_engine_function(native, version = "1.4")
  expect_equal(benchmark_linda_correction_count(native, effective), 1)
  expect_identical(body(effective)[[2L]], body(native)[[2L]])
  body(effective)[[2L]][[3L]][[2L]] <- NULL
  expect_identical(benchmark_linda_correction_count(native, effective), Inf)
})

test_that("LinDA audit verifies executed imputation or pseudocount without refitting", {
  set.seed(173)
  saved_rng <- .Random.seed
  imputed <- benchmark_linda_zero_handling(protection_linda_fixture("imputation"))
  pseudo_counted <- benchmark_linda_zero_handling(protection_linda_fixture("pseudo_count"))
  inactive <- benchmark_linda_zero_handling(protection_linda_fixture(zeros = FALSE))
  expect_true(imputed$passed)
  expect_identical(imputed$mode, "imputation")
  expect_true(pseudo_counted$passed)
  expect_identical(pseudo_counted$mode, "pseudo_count")
  expect_true(inactive$passed)
  expect_identical(inactive$mode, "no_zeros")
  expect_identical(.Random.seed, saved_rng)
  invalid <- protection_linda_fixture()
  invalid$result$feature.dat.use[2L, 2L] <- 999
  expect_false(benchmark_linda_zero_handling(invalid)$passed)
  invalid$result$feature.dat.use <- NULL
  expect_false(benchmark_linda_zero_handling(invalid)$passed)
})

test_that("LinDA native checks reject old unpatched dispatchers and invalid processed counts", {
  fit <- protection_linda_fixture()
  checks <- benchmark_native_checks("linda_pilot", list(fit), data.frame(),
                                     protection_simulation_fixture())
  expect_true(all(checks$passed))
  fit$effective_function <- benchmark_engine_implementation("linda", "fit", MicrobiomeStat::linda)
  checks <- benchmark_native_checks("linda_pilot", list(fit), data.frame(),
                                     protection_simulation_fixture())
  expect_false(checks$passed[checks$check == "linda_corrected_engine"])
  fit$result$feature.dat.use[1L, 1L] <- -1
  checks <- benchmark_native_checks("linda_pilot", list(fit), data.frame(),
                                     protection_simulation_fixture())
  expect_false(checks$passed[checks$check == "linda_executed_zero_handling"])
})
