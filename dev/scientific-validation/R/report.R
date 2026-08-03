read_validation_artifact <- function(input_dir, filename) {
  paths <- list.files(
    input_dir,
    pattern = paste0("^", gsub("[.]", "[.]", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  tables <- lapply(paths, function(path) {
    tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(cnd) data.frame()
    )
  })
  dplyr::bind_rows(tables)
}

aggregate_validation_artifacts <- function(input_dir, output_dir,
                                           render = TRUE, strict = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  artifact_names <- c(
    "results", "runs", "truth", "manifests", "metrics", "summary", "gates"
  )
  artifacts <- stats::setNames(lapply(artifact_names, function(name) {
    read_validation_artifact(input_dir, paste0(name, ".csv"))
  }), artifact_names)
  artifacts$truth <- dplyr::distinct(artifacts$truth)
  artifacts$manifests <- dplyr::distinct(artifacts$manifests)
  for (name in artifact_names) {
    utils::write.csv(
      artifacts[[name]], file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE, na = ""
    )
  }
  saveRDS(artifacts, file.path(output_dir, "validation-results.rds"))

  if (render) {
    rmarkdown::render(
      file.path(validation_root_path(), "report", "scientific-validation.Rmd"),
      output_file = "scientific-validation.html",
      output_dir = normalizePath(output_dir),
      params = list(data_dir = normalizePath(output_dir)),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
  }

  gate_failure <- nrow(artifacts$gates) > 0L &&
    any(artifacts$gates$status == "fail", na.rm = TRUE)
  engine_failure <- nrow(artifacts$runs) > 0L &&
    any(artifacts$runs$status == "engine_error", na.rm = TRUE)
  list(
    status = if (strict && (gate_failure || engine_failure)) 1L else 0L,
    gate_failure = gate_failure,
    engine_failure = engine_failure,
    artifacts = artifacts
  )
}
