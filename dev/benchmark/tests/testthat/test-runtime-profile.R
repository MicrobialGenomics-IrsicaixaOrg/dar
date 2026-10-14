test_that("the method profile identifies robust ANCOM and corrected LinDA", {
  skip_if_not_installed("ANCOMBC")
  skip_if_not_installed("MicrobiomeStat")
  profile <- benchmark_method_profile()

  expect_identical(profile$protocol_revision, "0.2")
  expect_identical(profile$ancom$pseudo_sens, TRUE)
  expect_identical(profile$ancom$decision, "diff_robust")
  expect_identical(profile$linda$adaptive, TRUE)
  expect_true(profile$linda$implementation$linda_adaptive_corrected)
  expect_identical(
    profile$linda$implementation$effective_sha256,
    benchmark_engine_implementation("linda", "fit")$effective_sha256
  )
})

test_that("an uncorrected LinDA implementation cannot identify a ready pilot", {
  testthat::local_mocked_bindings(
    benchmark_engine_implementation = function(...) {
      list(linda_adaptive_corrected = FALSE)
    },
    .env = .GlobalEnv
  )
  expect_error(benchmark_method_profile(), "requires the corrected LinDA")
})

test_that("method decisions and effective implementations affect the receipt hash", {
  profile <- list(
    protocol_revision = "0.2",
    ancom = list(pseudo_sens = TRUE, decision = "diff_robust"),
    linda = list(implementation = list(effective_sha256 = "corrected"))
  )
  testthat::local_mocked_bindings(
    benchmark_method_profile = function() profile,
    benchmark_version_manifest = function() data.frame(
      package = "unit-fixture", version = "1", build_identity = "test"
    ),
    .env = .GlobalEnv
  )
  reservoir <- list(checksum = "unit-reservoir")
  robust <- benchmark_contract(repository_root, benchmark_root, reservoir)
  expect_identical(robust$contract$schema_version, 2L)
  expect_identical(robust$contract$method_profile, profile)

  profile$ancom <- list(pseudo_sens = FALSE, decision = "diff")
  basic <- benchmark_contract(repository_root, benchmark_root, reservoir)
  expect_false(identical(robust$fingerprint, basic$fingerprint))
  profile$ancom <- robust$contract$method_profile$ancom
  profile$linda$implementation$effective_sha256 <- "uncorrected"
  uncorrected <- benchmark_contract(repository_root, benchmark_root, reservoir)
  expect_false(identical(robust$fingerprint, uncorrected$fingerprint))

  directory <- withr::local_tempdir()
  saveRDS(list(mode = "preflight", passed = TRUE, fingerprint = basic$fingerprint),
          file.path(directory, "preflight-receipt.rds"))
  arguments <- list(mode = "pilot", confirm_pilot = TRUE, preflight_dir = directory)
  expect_error(benchmark_check_pilot_authority(arguments, robust$fingerprint),
               "different source/configuration/runtime")
})

test_that("future reports label the amended method profile explicitly", {
  directory <- withr::local_tempdir()
  profile <- list(protocol_revision = "0.2", linda = list(
    implementation = list(implementation = "dar_local_adaptive_zero_fix")
  ))
  runs <- data.frame(procedure = character(), reason = character(),
                     status = character())
  checks <- data.frame(check = character(), passed = logical(), detail = character())
  benchmark_write_report(directory, runs, checks, "preflight", method_profile = profile)
  report <- paste(readLines(file.path(directory, "report.md")), collapse = "\n")

  expect_match(report, "Protocol revision: **0.2**", fixed = TRUE)
  expect_match(report, "native `diff_robust`", fixed = TRUE)
  expect_match(report, "dar_local_adaptive_zero_fix", fixed = TRUE)
})
