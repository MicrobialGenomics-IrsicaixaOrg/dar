test_that("CLI defaults to preflight and rejects ambiguous/full-pilot calls", {
  args <- benchmark_parse_arguments("--output_dir=new-run")
  expect_identical(args$mode, "preflight")
  expect_false(args$confirm_pilot)
  expect_error(benchmark_parse_arguments(character()), "output_dir")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--output_dir=b")), "duplicate")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--mode=full")), "Invalid mode")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--confirm_pilot=yes")), "Invalid mode")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--mode=pilot")), "locked")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--base_seed=1")), "Unknown")
  expect_error(benchmark_parse_arguments(c("--output_dir=a", "--confirm_pilot=true")), "not applicable")
})

test_that("native job budgets and repeat seeds are frozen without running fits", {
  preflight <- benchmark_job_plan(TRUE)
  pilot <- benchmark_job_plan(FALSE)
  expect_equal(nrow(preflight), 20L)
  expect_equal(nrow(pilot), 860L)
  expect_setequal(preflight$cell_id, c("PF_G1", "PF_G2"))
  expect_true(all(preflight$mode == "preflight"))
  expect_true(all(pilot$mode == "pilot"))
  expect_false(any(pilot$cell_id %in% c("P5", "P6") & pilot$engine == "wilcox_pilot"))
  keys <- c("cell_id", "replicate", "attempt", "engine")
  expect_equal(anyDuplicated(pilot[keys]), 0L)
  one <- preflight[preflight$attempt == 1L, ]
  two <- preflight[preflight$attempt == 2L, ]
  expect_identical(one$engine_seed, two$engine_seed)
  expect_identical(one$data_seed, two$data_seed)
})

test_that("only a passing matching receipt unlocks the full pilot", {
  directory <- withr::local_tempdir()
  arguments <- list(mode = "pilot", confirm_pilot = TRUE, preflight_dir = directory)
  expect_error(benchmark_check_pilot_authority(arguments, "current"), "receipt")
  file <- file.path(directory, "preflight-receipt.rds")
  saveRDS(list(mode = "preflight", passed = FALSE, fingerprint = "current"), file)
  expect_error(benchmark_check_pilot_authority(arguments, "current"), "failed")
  saveRDS(list(mode = "preflight", passed = TRUE, fingerprint = "old"), file)
  expect_error(benchmark_check_pilot_authority(arguments, "current"), "different")
  saveRDS(list(mode = "preflight", passed = TRUE, fingerprint = "current"), file)
  expect_true(benchmark_check_pilot_authority(arguments, "current"))
})

test_that("an install timestamp is not mistaken for scientific runtime identity", {
  expect_identical(
    benchmark_build_identity("R 4.6.1; aarch64; 2026-10-14 12:00:00 UTC; unix"),
    benchmark_build_identity("R 4.6.1; aarch64; 2026-10-14 13:00:00 UTC; unix")
  )
  expect_false(identical(
    benchmark_build_identity("R 4.6.1; aarch64; time; unix"),
    benchmark_build_identity("R 4.5.1; aarch64; time; unix")
  ))
})

test_that("combining isolated result fixtures preserves actual directional bake", {
  simulation <- benchmark_simulate("PF_G1", preflight = TRUE)
  outcomes <- stats::setNames(lapply(benchmark_engines(), function(id) {
    benchmark_prepared_fixture(simulation, id)
  }), benchmark_engines())
  prepared <- benchmark_combine_prepared(outcomes, simulation)
  expect_setequal(prepared@execution$executed_steps, benchmark_engines())
  majority <- benchmark_consensus(prepared, 3L, simulation)
  unanimity <- benchmark_consensus(prepared, 4L, simulation)
  expect_equal(nrow(majority$results), 1L)
  expect_equal(nrow(unanimity$results), 1L)
  expect_true(majority$checked)
  expect_true(all(majority$selected$method_count == 4))
  outcomes$maaslin_pilot <- benchmark_prepared_fixture(simulation, "maaslin_pilot", direction = -1)
  prepared <- benchmark_combine_prepared(outcomes, simulation)
  expect_equal(nrow(benchmark_consensus(prepared, 3L, simulation)$results), 1L)
  expect_equal(nrow(benchmark_consensus(prepared, 4L, simulation)$results), 0L)
  outcomes$maaslin_pilot$status <- "engine_error"
  expect_error(benchmark_combine_prepared(outcomes, simulation), "unavailable")
})

test_that("an unavailable panel is never a successful empty consensus", {
  simulation <- benchmark_simulate("PF_G1", preflight = TRUE)
  outcomes <- stats::setNames(lapply(benchmark_engines(), function(id) {
    list(status = "engine_error")
  }), benchmark_engines())
  consensus <- benchmark_dataset_consensus(outcomes, simulation, attempt = 1L)
  expect_true(all(consensus$runs$status == "unavailable"))
  expect_true(all(is.na(consensus$metrics$fdp)))
  expect_true(all(consensus$metrics$conservative_directional_power == 0))
})

test_that("native dispatcher observation does not run or modify a statistical fit", {
  capture <- new.env(parent = emptyenv())
  capture$stages <- list()
  assign(".benchmark_capture", capture, envir = .GlobalEnv)
  withr::defer(rm(".benchmark_capture", envir = .GlobalEnv))
  testthat::local_mocked_bindings(
    engine_stage_function = function(...) function(x, exact) list(value = x, exact = exact),
    .package = "dar"
  )
  suppressMessages(trace("exec_engine_stage", where = asNamespace("dar"), print = FALSE,
    tracer = quote(get("benchmark_trace_enter", envir = .GlobalEnv)(environment())),
    exit = quote(get("benchmark_trace_exit", envir = .GlobalEnv)(environment(), returnValue()))))
  withr::defer(suppressMessages(untrace("exec_engine_stage", where = asNamespace("dar"))))
  value <- dar:::exec_engine_stage("wilcox", "test", list(test = list(exact = FALSE)),
                                   fixed = list(x = 7))
  expect_identical(value, list(value = 7, exact = FALSE))
  expect_true(capture$stages[[1L]]$completed)
  expect_identical(capture$stages[[1L]]$result, value)
  expect_identical(capture$stages[[1L]]$arguments, list(x = 7, exact = FALSE))
})
