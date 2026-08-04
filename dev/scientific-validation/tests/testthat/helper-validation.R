validation_root <- normalizePath(
  file.path("..", ".."),
  mustWork = TRUE
)
options(dar.validation.root = validation_root)

source(file.path(validation_root, "R", "scenarios.R"), local = FALSE)
source(file.path(validation_root, "R", "simulate.R"), local = FALSE)
source(file.path(validation_root, "R", "engines.R"), local = FALSE)
source(file.path(validation_root, "R", "metrics.R"), local = FALSE)
source(file.path(validation_root, "R", "baseline.R"), local = FALSE)
source(file.path(validation_root, "R", "cli.R"), local = FALSE)
source(file.path(validation_root, "R", "report.R"), local = FALSE)

validation_metric_fixture <- function(reverse = FALSE) {
  truth <- expand.grid(
    taxa_id = sprintf("taxon_%02d", 1:6),
    contrast_id = "condition[treated-control]",
    stringsAsFactors = FALSE
  )
  truth$scenario <- "cross_sectional_signal"
  truth$replicate <- 1L
  truth$true_effect <- c(1, -1, 0, 0, 0, 0)
  truth$is_null <- truth$true_effect == 0
  truth$expected_direction <- sign(truth$true_effect)
  truth$truth_scale <- "centered_log_abundance"
  effects <- c(1.4, -1.4, 0.1, 0.1, 0.1, 0.1)
  if (reverse) effects[1:2] <- -effects[1:2]
  results <- data.frame(
    engine = "deseq",
    scenario = truth$scenario,
    replicate = truth$replicate,
    taxa_id = truth$taxa_id,
    contrast_id = truth$contrast_id,
    effect = effects,
    padj = c(0.01, 0.02, 0.01, 0.5, 0.5, 0.5),
    std_error = 0.2,
    effect_metric = "log2_fold_change",
    truth_multiplier = 1 / log(2),
    comparable = TRUE
  )
  list(results = results, truth = truth)
}
