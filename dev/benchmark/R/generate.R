# Dataset generation and the empirical reservoir are isolated from engine fits.
benchmark_metadata <- function(spec) {
  n <- spec$n_per_group
  metadata <- data.frame(
    sample_id = sprintf("sample_%03d", seq_len(2L * n)),
    condition = factor(rep(c("control", "treated"), each = n),
                       levels = c("control", "treated")),
    stringsAsFactors = FALSE
  )
  if (spec$adjust_batch) {
    if (n != 30L) {
      stop("The frozen confounded design requires 30 samples per group.",
           call. = FALSE)
    }
    metadata$batch <- factor(c(
      sample(c(rep("A", 24L), rep("B", 6L))),
      sample(c(rep("A", 6L), rep("B", 24L)))
    ), levels = c("A", "B"))
  }
  rownames(metadata) <- metadata$sample_id
  metadata
}

benchmark_effects <- function(spec, taxa_ids) {
  positions <- benchmark_signal_positions()
  beta <- gamma <- stats::setNames(rep(0, spec$n_taxa), taxa_ids)
  if (identical(spec$condition_signal, "balanced")) {
    beta[positions$condition] <- rep(c(1, -1), 5L) * spec$effect_magnitude
  } else if (identical(spec$condition_signal, "unidirectional")) {
    beta[positions$condition] <- spec$effect_magnitude
  } else if (!identical(spec$condition_signal, "null")) {
    stop("Unknown signal specification.", call. = FALSE)
  }
  if (spec$adjust_batch) {
    gamma[positions$batch] <- rep(c(1, -1), 5L) * log(2)
  }
  list(beta = beta, gamma = gamma)
}

benchmark_make_recipe <- function(counts, metadata, adjust_batch = FALSE) {
  if (nrow(counts) < 2L || anyDuplicated(rownames(counts))) {
    stop("The common universe needs at least two unique taxon IDs.", call. = FALSE)
  }
  # The source lineage is recorded separately. Unique single-rank labels stop
  # a second aggregation of distinct selected lineages during synthetic fits.
  taxonomy <- matrix(
    paste0("benchmark_", rownames(counts)), ncol = 1L,
    dimnames = list(rownames(counts), "Species")
  )
  samples <- metadata[, setdiff(names(metadata), c("sample_id", "source_profile")),
                      drop = FALSE]
  phy <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy), phyloseq::sample_data(samples)
  )
  rec <- dar::add_model(
    dar::recipe(phy),
    formula = if (adjust_batch) ~ condition + batch else ~ condition,
    targets = "condition", tax_level = "Species",
    reference = c(condition = "control"), contrasts = "simple",
    na_action = "error"
  )
  list(phy = phy, recipe = rec)
}

benchmark_truth <- function(beta, retained, contrast_id) {
  if (length(retained) != length(beta) || !any(retained)) {
    stop("The common universe must retain at least one original taxon.",
         call. = FALSE)
  }
  clr <- median <- rep(NA_real_, length(beta))
  clr[retained] <- beta[retained] - mean(beta[retained])
  median[retained] <- beta[retained] - stats::median(beta[retained])
  data.frame(
    taxa_id = names(beta), contrast_id = contrast_id,
    beta_latent = as.numeric(beta),
    beta_clr_original = as.numeric(beta - mean(beta)),
    beta_median_original = as.numeric(beta - stats::median(beta)),
    beta_clr_U = clr, beta_median_U = median,
    retained = retained, stringsAsFactors = FALSE
  )
}

benchmark_prepare_reservoir <- function(output_dir = NULL) {
  source <- new.env(parent = emptyenv())
  utils::data("metaHIV_phy", package = "dar", envir = source)
  phy <- source$metaHIV_phy
  if (is.null(phy)) stop("The bundled metaHIV_phy is unavailable.", call. = FALSE)
  metadata <- as(phyloseq::sample_data(phy), "data.frame")
  # Clinical categories and generic sample IDs are not donor identifiers. The
  # currently bundled data have no unequivocal subject/donor column.
  donor_candidates <- intersect(
    c("donor_id", "subject_id", "participant_id", "patient_id"), names(metadata)
  )
  donor_column <- NULL
  if (length(donor_candidates) == 1L) {
    donor_values <- as.character(metadata[[donor_candidates]])
    if (!anyNA(donor_values) && all(nzchar(trimws(donor_values)))) {
      donor_column <- donor_candidates
    }
  }
  sample_ids <- sort(rownames(metadata), method = "radix")
  if (!is.null(donor_column)) {
    sample_ids <- sample_ids[!duplicated(metadata[sample_ids, donor_column])]
  }
  phy <- phyloseq::prune_samples(sample_ids, phy)
  metadata <- metadata[sample_ids, , drop = FALSE]
  # This model exists solely to invoke dar's canonical modeled analysis view.
  metadata$benchmark_reservoir_group <- factor(
    rep(c("A", "B"), length.out = nrow(metadata)), levels = c("A", "B")
  )
  phyloseq::sample_data(phy) <- phyloseq::sample_data(metadata)
  source_recipe <- dar::add_model(
    dar::recipe(phy), ~ benchmark_reservoir_group,
    targets = "benchmark_reservoir_group", tax_level = "Species",
    reference = c(benchmark_reservoir_group = "A"), na_action = "error"
  )
  aggregated <- dar:::analysis_phyloseq(source_recipe)
  counts <- dar:::recipe_count_matrix(source_recipe, phy = aggregated)
  counts <- counts[, sample_ids, drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || any(colSums(counts) <= 0)) {
    stop("The empirical source must contain finite non-negative integer counts.",
         call. = FALSE)
  }
  relative <- sweep(counts, 2L, colSums(counts), "/")
  statistics <- data.frame(
    taxa_id = rownames(counts), prevalence = rowMeans(counts > 0),
    mean_relative_abundance = rowMeans(relative), stringsAsFactors = FALSE
  )
  eligible <- statistics[statistics$prevalence >= 0.1, , drop = FALSE]
  if (nrow(eligible) < 100L) {
    stop("Fewer than 100 empirical taxa pass pooled prevalence >= 0.1; ",
         "a prospective protocol amendment is required.", call. = FALSE)
  }
  choose <- order(-eligible$prevalence, -eligible$mean_relative_abundance,
                  eligible$taxa_id, method = "radix")
  selected <- eligible[choose[seq_len(100L)], , drop = FALSE]
  selected <- selected[order(selected$mean_relative_abundance, selected$taxa_id,
                             method = "radix"), , drop = FALSE]
  rownames(selected) <- NULL
  counts <- counts[selected$taxa_id, , drop = FALSE]
  taxonomy <- as(phyloseq::tax_table(aggregated), "matrix")
  taxonomy <- taxonomy[selected$taxa_id, , drop = FALSE]
  metadata$benchmark_reservoir_group <- NULL
  reservoir <- list(counts = counts, taxonomy = taxonomy, metadata = metadata,
                    taxa = selected)
  reservoir$checksum <- digest::digest(reservoir, algo = "sha256")
  reservoir$manifest <- list(
    source = "dar::metaHIV_phy", samples = ncol(counts), taxa = nrow(counts),
    aggregated_source_taxa = nrow(statistics), eligible_taxa = nrow(eligible),
    donor_column = donor_column,
    donor_independence = if (is.null(donor_column)) {
      "independencia de donantes no verificable"
    } else "one sample per donor; lowest sample ID retained",
    pooled_prevalence_minimum = 0.1,
    abundance_definition = "mean relative abundance before selecting taxa",
    selection = "descending prevalence, descending mean relative abundance, ascending ID",
    signal_order = "ascending mean relative abundance, ascending ID",
    smoothing = 0.5, checksum = reservoir$checksum
  )
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(reservoir, file.path(output_dir, "reservoir.rds"))
    utils::write.csv(selected, file.path(output_dir, "reservoir-taxa.csv"),
                     row.names = FALSE)
    utils::write.csv(metadata, file.path(output_dir, "reservoir-samples.csv"),
                     row.names = TRUE)
    jsonlite::write_json(reservoir$manifest, file.path(output_dir, "reservoir-manifest.json"),
                        auto_unbox = TRUE, pretty = TRUE, null = "null")
  }
  reservoir
}

benchmark_validate_reservoir <- function(reservoir) {
  if (!is.list(reservoir) || !is.matrix(reservoir$counts) ||
      nrow(reservoir$counts) != 100L || ncol(reservoir$counts) < 1L ||
      any(!is.finite(reservoir$counts)) || any(reservoir$counts < 0) ||
      any(reservoir$counts != floor(reservoir$counts)) ||
      is.null(rownames(reservoir$counts)) || is.null(colnames(reservoir$counts)) ||
      anyDuplicated(rownames(reservoir$counts)) ||
      anyDuplicated(colnames(reservoir$counts))) {
    stop("G2 requires a frozen reservoir of exactly 100 uniquely identified count profiles.",
         call. = FALSE)
  }
  invisible(TRUE)
}

benchmark_simulate <- function(cell_id, replicate = 1L, base_seed = NULL,
                               reservoir = NULL, preflight = FALSE) {
  spec <- benchmark_case(cell_id, preflight = preflight)
  seed <- benchmark_data_seed(cell_id, replicate, base_seed, preflight)
  prior_rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(prior_rng)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else assign(".Random.seed", prior_rng, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  rng_start <- .Random.seed
  metadata <- benchmark_metadata(spec)
  n_samples <- nrow(metadata)
  if (identical(spec$generator, "G1")) {
    taxa_ids <- sprintf("taxon_%03d", seq_len(spec$n_taxa))
  } else {
    benchmark_validate_reservoir(reservoir)
    taxa_ids <- rownames(reservoir$counts)
  }
  coefficients <- benchmark_effects(spec, taxa_ids)
  condition <- as.integer(metadata$condition == "treated")
  if (identical(spec$generator, "G1")) {
    mu <- seq(log(0.02), log(20), length.out = spec$n_taxa)
    blocks <- matrix(stats::rnorm(n_samples * 10L), nrow = n_samples)
    individual <- matrix(stats::rnorm(n_samples * spec$n_taxa), nrow = n_samples)
    noise <- 0.8 * (sqrt(0.3) * blocks[, rep(seq_len(10L), each = 10L)] +
                     sqrt(0.7) * individual)
    latent <- matrix(rep(mu, each = n_samples), nrow = n_samples) + noise
    latent <- latent + outer(condition, coefficients$beta)
    if (spec$adjust_batch) {
      latent <- latent + outer(as.integer(metadata$batch == "B"), coefficients$gamma)
    }
  } else {
    profile_index <- sample.int(ncol(reservoir$counts), n_samples, replace = TRUE)
    metadata$source_profile <- colnames(reservoir$counts)[profile_index]
    smoothed <- reservoir$counts[, profile_index, drop = FALSE] + 0.5
    baseline <- sweep(smoothed, 2L, colSums(smoothed), "/")
    latent <- log(t(baseline)) + outer(condition, coefficients$beta)
  }
  # Log-softmax avoids numerical overflow without changing multinomial probabilities.
  proportions <- t(apply(latent, 1L, function(x) {
    abundance <- exp(x - max(x))
    abundance / sum(abundance)
  }))
  depths <- as.integer(pmax(500, round(exp(log(spec$depth_median) +
                                           0.5 * stats::rnorm(n_samples)))))
  counts_original <- vapply(seq_len(n_samples), function(index) {
    as.integer(stats::rmultinom(1L, depths[[index]], proportions[index, ]))
  }, integer(spec$n_taxa))
  dimnames(counts_original) <- list(taxa_ids, metadata$sample_id)
  retained <- rowSums(counts_original) > 0
  counts <- counts_original[retained, , drop = FALSE]
  objects <- benchmark_make_recipe(counts, metadata, adjust_batch = spec$adjust_batch)
  resolved <- dar:::resolve_model(objects$recipe)
  plan <- resolved$contrast_plan
  if (nrow(plan) != 1L || !identical(as.character(plan$var[[1L]]), "condition")) {
    stop("The frozen pilot must resolve to one condition contrast.", call. = FALSE)
  }
  truth <- benchmark_truth(coefficients$beta, retained, as.character(plan$contrast_id[[1L]]))
  truth$cell_id <- cell_id
  truth$replicate <- as.integer(replicate)
  manifest <- list(
    mode = spec$mode, cell_id = cell_id, replicate = as.integer(replicate),
    seed = seed, samples = n_samples, original_taxa = spec$n_taxa,
    retained_taxa = sum(retained), removed_all_zero_taxa = taxa_ids[!retained],
    specification = spec, contrast_id = as.character(plan$contrast_id[[1L]]),
    rng_kind = RNGkind(), rng_start = rng_start, rng_end = .Random.seed,
    counts_checksum = digest::digest(counts_original, algo = "sha256"),
    metadata_checksum = digest::digest(metadata, algo = "sha256"),
    reservoir_checksum = if (identical(spec$generator, "G2")) reservoir$checksum else NULL,
    clr_shift_original = mean(coefficients$beta),
    clr_shift_U = mean(coefficients$beta[retained]),
    median_shift_original = stats::median(coefficients$beta),
    median_shift_U = stats::median(coefficients$beta[retained])
  )
  c(objects, list(
    truth = truth, manifest = manifest, counts = counts,
    counts_original = counts_original, metadata = metadata,
    coefficients = coefficients, specification = spec, contrast_plan = plan
  ))
}
