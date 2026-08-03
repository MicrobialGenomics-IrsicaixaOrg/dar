validation_metadata <- function(spec) {
  if (identical(spec$design, "cross_sectional")) {
    metadata <- expand.grid(
      replicate = seq_len(24L),
      condition = c("control", "treated"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    metadata$batch <- ifelse(metadata$replicate %% 2L == 0L, "B", "A")
  } else if (identical(spec$design, "confounded")) {
    metadata <- data.frame(
      condition = rep(c("control", "treated"), each = 24L),
      batch = c(rep("A", 18L), rep("B", 6L), rep("A", 6L), rep("B", 18L)),
      stringsAsFactors = FALSE
    )
    metadata$replicate <- ave(
      seq_len(nrow(metadata)), metadata$condition, FUN = seq_along
    )
  } else {
    subjects <- expand.grid(
      subject_index = seq_len(12L),
      condition = c("control", "treated"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    subjects$subject <- paste0(substr(subjects$condition, 1L, 1L), subjects$subject_index)
    subjects$batch <- ifelse(subjects$subject_index %% 2L == 0L, "B", "A")
    metadata <- merge(
      subjects,
      data.frame(time = 0:2),
      all = TRUE,
      sort = FALSE
    )
    metadata <- metadata[order(metadata$condition, metadata$subject_index, metadata$time), ]
    if (identical(spec$design, "longitudinal_unbalanced")) {
      remove <- (metadata$subject %in% c("c2", "t3") & metadata$time == 2L) |
        (metadata$subject == "t10" & metadata$time == 1L)
      metadata <- metadata[!remove, , drop = FALSE]
    }
  }

  metadata$sample_id <- sprintf("sample_%03d", seq_len(nrow(metadata)))
  metadata$condition <- factor(metadata$condition, levels = c("control", "treated"))
  metadata$batch <- factor(metadata$batch, levels = c("A", "B"))
  if ("time" %in% names(metadata)) {
    metadata$time <- factor(metadata$time, levels = c(0, 1, 2))
  }
  metadata
}

validation_fixed_formula <- function(formula) {
  reformulas::nobars(formula)
}

validation_coefficients <- function(design, spec, n_taxa) {
  coefficients <- matrix(
    0,
    nrow = n_taxa,
    ncol = ncol(design),
    dimnames = list(sprintf("taxon_%02d", seq_len(n_taxa)), colnames(design))
  )
  coefficients[, "(Intercept)"] <- seq(log(0.5), log(20), length.out = n_taxa)

  set_balanced <- function(column, magnitude, positive, negative) {
    if (column %in% colnames(coefficients)) {
      coefficients[positive, column] <<- magnitude
      coefficients[negative, column] <<- -magnitude
    }
  }

  if (!identical(spec$signal, "null")) {
    set_balanced("conditiontreated", 0.9, 1:4, 5:8)
  }
  if (spec$confounders) {
    set_balanced("batchB", 0.55, 9:10, 11:12)
  }
  if (identical(spec$signal, "longitudinal")) {
    set_balanced("time1", 0.30, 13:14, 15:16)
    set_balanced("time2", 0.55, 13:14, 15:16)
    set_balanced("conditiontreated:time1", 0.45, 1:2, 5:6)
    set_balanced("conditiontreated:time2", 0.80, 1:2, 5:6)
  }
  coefficients
}

validation_softmax <- function(x) {
  shifted <- x - max(x)
  exp(shifted) / sum(exp(shifted))
}

simulate_validation_scenario <- function(
    scenario_id,
    replicate = 1L,
    base_seed = 140L,
    n_taxa = 60L
) {
  spec <- validation_scenario(scenario_id)
  seed <- validation_seed(scenario_id, replicate, base_seed)
  set.seed(seed)

  metadata_complete <- validation_metadata(spec)
  fixed_formula <- validation_fixed_formula(spec$formula)
  design <- stats::model.matrix(fixed_formula, data = metadata_complete)
  coefficients <- validation_coefficients(design, spec, n_taxa)
  latent <- design %*% t(coefficients)

  if (spec$random) {
    subject_levels <- unique(metadata_complete$subject)
    random_intercepts <- matrix(
      stats::rnorm(length(subject_levels) * n_taxa, sd = 0.30),
      nrow = length(subject_levels),
      dimnames = list(subject_levels, rownames(coefficients))
    )
    subject_index <- match(metadata_complete$subject, subject_levels)
    latent <- latent + random_intercepts[subject_index, , drop = FALSE]
  }

  latent <- latent + matrix(
    stats::rnorm(length(latent), sd = 0.35),
    nrow = nrow(latent)
  )
  library_sizes <- 8000L + stats::rpois(nrow(metadata_complete), lambda = 4000)
  counts <- vapply(seq_len(nrow(metadata_complete)), function(index) {
    as.integer(stats::rmultinom(
      1L,
      size = library_sizes[[index]],
      prob = validation_softmax(latent[index, ])
    ))
  }, integer(n_taxa))
  rownames(counts) <- rownames(coefficients)
  colnames(counts) <- metadata_complete$sample_id

  if (spec$sparse) {
    sparse_taxa <- seq.int(n_taxa - 9L, n_taxa)
    dropout <- matrix(
      stats::runif(length(sparse_taxa) * ncol(counts)) < 0.65,
      nrow = length(sparse_taxa)
    )
    sparse_counts <- counts[sparse_taxa, , drop = FALSE]
    sparse_counts[dropout] <- 0L
    counts[sparse_taxa, ] <- sparse_counts
  }

  metadata_observed <- metadata_complete
  if (spec$missing) {
    missing_rows <- c(2L, floor(nrow(metadata_observed) / 2L) + 1L)
    metadata_observed$batch[missing_rows] <- NA
  }
  rownames(metadata_observed) <- metadata_observed$sample_id
  sample_metadata <- metadata_observed
  sample_metadata$sample_id <- NULL
  taxonomy <- matrix(
    sprintf("Species_%02d", seq_len(n_taxa)),
    ncol = 1L,
    dimnames = list(rownames(counts), "Species")
  )
  phy <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(sample_metadata)
  )
  rec <- dar::recipe(phy)
  rec <- dar::add_model(
    rec,
    formula = spec$formula,
    targets = "condition",
    tax_level = "Species",
    time = spec$time,
    reference = c(condition = "control"),
    contrasts = spec$contrasts,
    na_action = spec$na_action
  )
  resolved <- dar:::resolve_model(rec)
  plan <- resolved$contrast_plan
  truth <- do.call(rbind, lapply(seq_len(nrow(plan)), function(index) {
    effect <- as.numeric(coefficients %*% plan$weights[[index]])
    data.frame(
      scenario = scenario_id,
      replicate = as.integer(replicate),
      taxa_id = rownames(coefficients),
      contrast_id = as.character(plan$contrast_id[[index]]),
      contrast_type = as.character(plan$contrast_type[[index]]),
      true_effect = effect,
      is_null = abs(effect) <= sqrt(.Machine$double.eps),
      expected_direction = sign(effect),
      truth_scale = "centered_log_abundance",
      stringsAsFactors = FALSE
    )
  }))

  list(
    phy = phy,
    recipe = rec,
    truth = truth,
    manifest = data.frame(
      scenario = scenario_id,
      replicate = as.integer(replicate),
      seed = seed,
      samples = phyloseq::nsamples(phy),
      complete_samples = nrow(resolved$data),
      dropped_samples = length(resolved$dropped_samples),
      taxa = phyloseq::ntaxa(phy),
      contrasts = nrow(plan),
      stringsAsFactors = FALSE
    ),
    coefficients = coefficients,
    specification = spec
  )
}
