# The checkout is installed into a run-local library, never the user's library.
benchmark_load <- function(root) {
  for (file in c("design", "generate", "metrics", "decision", "engines", "runtime", "runner")) {
    source(file.path(root, "R", paste0(file, ".R")), local = .GlobalEnv)
  }
  invisible(root)
}

benchmark_dependencies <- function() {
  c("ANCOMBC", "MicrobiomeStat", "ALDEx2", "maaslin3", "phyloseq", "mia",
    "microbiome", "callr", "ps", "digest", "jsonlite")
}

benchmark_prepare_runtime <- function(repository_root, output_dir) {
  missing <- benchmark_dependencies()[!vapply(benchmark_dependencies(), requireNamespace,
                                              logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing dependencies (nothing will be installed automatically): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if ("dar" %in% loadedNamespaces()) {
    stop("Start in a fresh vanilla R process so the checkout can be pinned.", call. = FALSE)
  }
  library <- file.path(output_dir, "runtime", "library")
  dir.create(library, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(output_dir, "runtime", "install.log")
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
      shQuote(paste0("--library=", library)), shQuote(repository_root)),
    stdout = log, stderr = log,
    env = c("R_PROFILE_USER=/dev/null", "R_ENVIRON_USER=/dev/null",
            "BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE"))
  if (status != 0L) stop("Installing the checkout failed; see runtime/install.log.", call. = FALSE)
  .libPaths(c(normalizePath(library), .libPaths()))
  if (!requireNamespace("dar", quietly = TRUE) ||
      normalizePath(find.package("dar")) != normalizePath(file.path(library, "dar"))) {
    stop("Refusing a fallback to a different installed dar.", call. = FALSE)
  }
  normalizePath(library)
}

benchmark_file_hashes <- function(files, root) {
  files <- sort(normalizePath(files, mustWork = TRUE))
  data.frame(file = substring(files, nchar(normalizePath(root)) + 2L),
             sha256 = vapply(files, function(file) digest::digest(file = file, algo = "sha256"),
                             character(1)), stringsAsFactors = FALSE)
}

benchmark_version_manifest <- function() {
  installed <- utils::installed.packages()
  packages <- unique(c("dar", benchmark_dependencies()))
  dependencies <- tools::package_dependencies(packages, db = installed,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
  packages <- sort(intersect(unique(c(packages, unlist(dependencies))), rownames(installed)))
  data.frame(package = packages, version = installed[packages, "Version"],
             build_identity = vapply(installed[packages, "Built"], benchmark_build_identity,
                                      character(1)), stringsAsFactors = FALSE)
}

benchmark_build_identity <- function(built) {
  parts <- strsplit(built, ";", fixed = TRUE)[[1L]]
  # R version, architecture and OS affect identity; INSTALL timestamp does not.
  if (length(parts) >= 4L) parts <- parts[c(1L, 2L, length(parts))]
  paste(trimws(parts), collapse = ";")
}

benchmark_method_profile <- function() {
  linda <- benchmark_engine_implementation("linda", "fit")
  if (!isTRUE(linda$linda_adaptive_corrected) ||
      !linda$implementation %in% c("native", "dar_local_adaptive_zero_fix") ||
      !isTRUE(linda$formals_unchanged) || !isTRUE(linda$namespace_unchanged)) {
    stop("The pilot requires the corrected LinDA adaptive-zero implementation.",
         call. = FALSE)
  }
  list(
    protocol_revision = "0.2",
    ancom = list(pseudo_sens = TRUE, decision = "diff_robust"),
    linda = list(adaptive = TRUE, implementation = linda),
    ancom_implementation = benchmark_engine_implementation("ancom", "fit")
  )
}

benchmark_contract <- function(repository_root, code_root, reservoir) {
  source_files <- c(file.path(repository_root, c("DESCRIPTION", "NAMESPACE")),
                    list.files(file.path(repository_root, "R"), full.names = TRUE,
                               pattern = "[.]R$"))
  code_files <- c(file.path(code_root, "run-benchmark.R"),
                  list.files(file.path(code_root, "R"), full.names = TRUE, pattern = "[.]R$"))
  contract <- list(
    format = "dar-consensus-pilot", schema_version = 2L,
    method_profile = benchmark_method_profile(),
    R_version = R.version.string, platform = R.version$platform,
    dar_source = benchmark_file_hashes(source_files, repository_root),
    benchmark_source = benchmark_file_hashes(code_files, code_root),
    protocol_sha256 = digest::digest(file = file.path(code_root, "pilot-protocol.md"), algo = "sha256"),
    packages = benchmark_version_manifest(), reservoir_checksum = reservoir$checksum
  )
  list(contract = contract, fingerprint = digest::digest(contract, algo = "sha256"))
}

benchmark_write_json <- function(value, file) {
  jsonlite::write_json(value, file, auto_unbox = TRUE, pretty = TRUE,
                        null = "null", na = "null", digits = NA)
}

benchmark_ram_limit <- function() {
  system <- Sys.info()[["sysname"]]
  bytes <- if (identical(system, "Darwin")) {
    suppressWarnings(as.numeric(system2("/usr/sbin/sysctl", c("-n", "hw.memsize"),
                                          stdout = TRUE, stderr = FALSE)))
  } else if (file.exists("/proc/meminfo")) {
    line <- grep("^MemTotal:", readLines("/proc/meminfo"), value = TRUE)
    as.numeric(sub("^MemTotal:[[:space:]]+([0-9]+).*", "\\1", line)) * 1024
  } else NA_real_
  if (length(bytes) != 1L || !is.finite(bytes) || bytes <= 0) {
    stop("Cannot verify physical RAM; no native fits will be launched.", call. = FALSE)
  }
  0.75 * bytes / 1024^2
}

benchmark_cost_projection <- function(runs, plan) {
  cost <- ifelse(is.finite(runs$user_seconds + runs$system_seconds),
                  runs$user_seconds + runs$system_seconds, runs$worker_elapsed_seconds)
  # Within each measured cell, use its observed per-engine mean. Unmeasured
  # cells use the same engine's measured mean and remain explicitly provisional.
  expected <- vapply(seq_len(nrow(plan)), function(i) {
    rows <- runs$procedure == plan$engine[[i]] & runs$cell_id == plan$cell_id[[i]]
    if (!any(rows)) rows <- runs$procedure == plan$engine[[i]]
    available <- cost[rows & is.finite(cost)]
    if (length(available)) mean(available) else NA_real_
  }, numeric(1))
  if (any(!is.finite(expected))) return(NA_real_)
  sum(expected) / 3600
}

benchmark_check_pilot_authority <- function(arguments, fingerprint) {
  if (arguments$mode != "pilot") return(invisible(TRUE))
  if (!isTRUE(arguments$confirm_pilot) || !nzchar(arguments$preflight_dir)) {
    stop("Full pilot requires explicit confirmation and a passing preflight directory.", call. = FALSE)
  }
  receipt_file <- file.path(arguments$preflight_dir, "preflight-receipt.rds")
  if (!file.exists(receipt_file)) stop("No preflight receipt found.", call. = FALSE)
  receipt <- readRDS(receipt_file)
  if (!isTRUE(receipt$passed) || !identical(receipt$mode, "preflight") ||
      !identical(receipt$fingerprint, fingerprint)) {
    stop("Preflight failed or belongs to a different source/configuration/runtime.", call. = FALSE)
  }
  invisible(TRUE)
}

benchmark_preflight_checks <- function(outcomes, consensus_records) {
  checks <- list()
  append <- function(check, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
    )
  }
  append("fit_budget", length(outcomes) == 20L, "Exactly two datasets x five engines x two attempts.")
  append("native_fits", length(outcomes) == 20L &&
           all(vapply(outcomes, function(x) x$status == "success", logical(1))),
         "Every planned engine attempt must succeed; failures remain visible.")
  for (cell in c("PF_G1", "PF_G2")) {
    for (engine in benchmark_engines()) {
      one <- outcomes[[paste(cell, 1L, engine, 1L, sep = "__")]]
      two <- outcomes[[paste(cell, 1L, engine, 2L, sep = "__")]]
      ok <- !is.null(one) && !is.null(two) && one$status == "success" && two$status == "success"
      if (ok) {
        a <- one$results; b <- two$results
        order_rows <- function(x) x[order(x$taxa_id, x$contrast_id), , drop = FALSE]
        ok <- isTRUE(all.equal(order_rows(a), order_rows(b), tolerance = 1e-10,
                               check.attributes = TRUE))
      }
      append(paste("repeat", cell, engine, sep = "__"), ok,
             "Same dataset and seed must reproduce native results, not runtime measurements.")
    }
  }
  append("actual_consensus", nrow(consensus_records) == 8L &&
           all(consensus_records$status == "success"),
         "Actual majority and unanimity must match explicit voting for both attempts/datasets.")
  do.call(rbind, checks)
}
