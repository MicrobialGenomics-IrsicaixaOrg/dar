# Frozen pilot design. These helpers do not run statistical engines.
benchmark_design <- function() {
  data.frame(
    cell_id = paste0("P", seq_len(9L)),
    cell_index = seq_len(9L),
    generator = c(rep("G1", 6L), rep("G2", 3L)),
    n_per_group = c(30L, 30L, 15L, rep(30L, 6L)),
    depth_median = c(15000L, 15000L, 1500L, rep(15000L, 6L)),
    condition_signal = c(
      "null", "balanced", "balanced", "unidirectional", "null",
      "balanced", "null", "balanced", "unidirectional"
    ),
    effect_magnitude = c(
      0, log(2), log(2), log(4), 0, log(2), 0, log(2), log(4)
    ),
    adjust_batch = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
    n_taxa = 100L,
    n_replicates = 20L,
    mode = "pilot",
    stringsAsFactors = FALSE
  )
}

benchmark_preflight_design <- function() {
  data.frame(
    cell_id = c("PF_G1", "PF_G2"),
    cell_index = seq_len(2L),
    generator = c("G1", "G2"),
    n_per_group = 30L,
    depth_median = 15000L,
    condition_signal = "balanced",
    effect_magnitude = log(4),
    adjust_batch = FALSE,
    n_taxa = 100L,
    n_replicates = 1L,
    mode = "preflight",
    stringsAsFactors = FALSE
  )
}

benchmark_case <- function(cell_id, preflight = FALSE) {
  if (!is.logical(preflight) || length(preflight) != 1L || is.na(preflight)) {
    stop("preflight must be a single non-missing logical value.", call. = FALSE)
  }
  design <- if (preflight) benchmark_preflight_design() else benchmark_design()
  if (!is.character(cell_id) || length(cell_id) != 1L ||
      is.na(cell_id) || !cell_id %in% design$cell_id) {
    stop("Unknown cell for the selected benchmark mode.", call. = FALSE)
  }
  design[match(cell_id, design$cell_id), , drop = FALSE]
}

benchmark_integer <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(name, " must be a valid integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

benchmark_data_seed <- function(cell_id, replicate = 1L, base_seed = NULL,
                                preflight = FALSE) {
  spec <- benchmark_case(cell_id, preflight = preflight)
  replicate <- benchmark_integer(replicate, "replicate")
  if (replicate > spec$n_replicates) {
    stop("Replicate exceeds the frozen design for this mode.", call. = FALSE)
  }
  if (is.null(base_seed)) {
    base_seed <- if (preflight) 19800310L else 20261014L
  }
  base_seed <- benchmark_integer(base_seed, "base_seed", minimum = 0L)
  benchmark_integer(base_seed + 10000 * spec$cell_index + replicate, "seed")
}

benchmark_signal_positions <- function() {
  list(
    condition = c(21L, 27L, 33L, 39L, 45L, 51L, 57L, 63L, 69L, 75L),
    batch = c(22L, 28L, 34L, 40L, 46L, 52L, 58L, 64L, 70L, 76L)
  )
}

benchmark_seed_ledger <- function(preflight = FALSE, base_seed = NULL) {
  design <- if (preflight) benchmark_preflight_design() else benchmark_design()
  rows <- lapply(seq_len(nrow(design)), function(index) {
    replicates <- seq_len(design$n_replicates[[index]])
    data.frame(
      mode = design$mode[[index]],
      cell_id = design$cell_id[[index]],
      replicate = replicates,
      data_seed = vapply(replicates, function(replicate) {
        benchmark_data_seed(
          design$cell_id[[index]], replicate, base_seed, preflight
        )
      }, integer(1L)),
      stringsAsFactors = FALSE
    )
  })
  ledger <- do.call(rbind, rows)
  for (index in seq_len(5L)) {
    ledger[[c("ancom", "linda", "aldex", "maaslin", "wilcox")[[index]]]] <-
      ledger$data_seed + 100000L * index
  }
  ledger
}
