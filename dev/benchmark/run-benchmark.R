#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Unable to locate benchmark entrypoint.", call. = FALSE)
benchmark_root <- normalizePath(dirname(sub("^--file=", "", script_argument)), mustWork = TRUE)
repository_root <- normalizePath(file.path(benchmark_root, "..", ".."), mustWork = TRUE)
source(file.path(benchmark_root, "R", "runtime.R"), local = .GlobalEnv)
benchmark_load(benchmark_root)
arguments <- benchmark_parse_arguments(commandArgs(trailingOnly = TRUE))
outcome <- benchmark_run(arguments, benchmark_root, repository_root)
message("Benchmark artifacts: ", outcome$output_dir)
message("Mode: ", arguments$mode, "; execution status: ", outcome$status)
quit(status = outcome$status, save = "no")
