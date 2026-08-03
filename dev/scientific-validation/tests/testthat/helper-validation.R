validation_root <- normalizePath(
  file.path("..", ".."),
  mustWork = TRUE
)
options(dar.validation.root = validation_root)

source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "simulate.R"), local = FALSE)
source(file.path(validation_root, "R", "engines.R"), local = FALSE)
source(file.path(validation_root, "R", "metrics.R"), local = FALSE)
source(file.path(validation_root, "R", "cli.R"), local = FALSE)
