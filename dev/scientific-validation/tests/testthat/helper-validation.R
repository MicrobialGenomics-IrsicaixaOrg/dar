validation_root <- normalizePath(
  file.path("..", ".."),
  mustWork = TRUE
)

source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "simulate.R"), local = FALSE)
