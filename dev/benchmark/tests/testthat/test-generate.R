benchmark_test_reservoir <- function() {
  counts <- matrix(rep(seq_len(100L), 4L), nrow = 100L,
                   dimnames = list(sprintf("source_%03d", seq_len(100L)),
                                   paste0("profile_", seq_len(4L))))
  list(counts = counts, checksum = digest::digest(counts, algo = "sha256"))
}

test_that("G1 generation preserves exactly 100 original taxa and count depths", {
  first <- benchmark_simulate("P2")
  second <- benchmark_simulate("P2")
  expect_identical(first$counts_original, second$counts_original)
  expect_identical(first$metadata, second$metadata)
  expect_equal(dim(first$counts_original), c(100L, 60L))
  expect_true(all(colSums(first$counts_original) >= 500L))
  expect_equal(rownames(first$counts), first$truth$taxa_id[first$truth$retained])
  expect_equal(phyloseq::taxa_names(dar:::analysis_phyloseq(first$recipe)),
               rownames(first$counts))
  expect_equal(first$truth$beta_latent[benchmark_signal_positions()$condition],
               rep(c(1, -1), 5L) * log(2))
  expect_equal(first$truth$beta_clr_original, first$truth$beta_latent,
               tolerance = 1e-10)
  expect_equal(first$truth$beta_median_original, first$truth$beta_latent)
  expect_equal(unique(first$truth$contrast_id), first$contrast_plan$contrast_id)
  expect_false(identical(first$counts_original,
                         benchmark_simulate("P2", replicate = 2L)$counts_original))
})

test_that("generation does not consume the caller's RNG state", {
  set.seed(71L)
  prior <- .Random.seed
  benchmark_simulate("P1")
  expect_identical(.Random.seed, prior)
})

test_that("nuls and adjusted confounding follow the known conditional truth", {
  for (cell_id in c("P1", "P5")) {
    data <- benchmark_simulate(cell_id)
    expect_true(all(data$truth$beta_latent == 0))
    expect_true(all(data$truth$beta_clr_original == 0))
  }
  adjusted <- benchmark_simulate("P6")
  cells <- table(adjusted$metadata$condition, adjusted$metadata$batch)
  expect_equal(unname(cells), matrix(c(24L, 6L, 6L, 24L), nrow = 2L))
  expect_equal(sum(adjusted$coefficients$gamma > 0), 5L)
  expect_equal(sum(adjusted$coefficients$gamma < 0), 5L)
  expect_length(intersect(which(adjusted$coefficients$beta != 0),
                          which(adjusted$coefficients$gamma != 0)), 0L)
  expect_equal(deparse(dar::get_model(adjusted$recipe)$formula), "~condition + batch")
})

test_that("unidirectional truth distinguishes latent and CLR objectives", {
  data <- benchmark_simulate("P4")
  expect_equal(sum(abs(data$truth$beta_latent) > 1e-10), 10L)
  expect_equal(sum(abs(data$truth$beta_clr_original) > 1e-10), 100L)
  expect_equal(data$manifest$clr_shift_original, 0.1 * log(4))
  expect_equal(data$truth$beta_median_original, data$truth$beta_latent)
  removed <- c(FALSE, rep(TRUE, 99L))
  beta <- stats::setNames(c(log(2), -log(2), rep(0, 98L)), paste0("taxon", 1:100))
  truth <- benchmark_truth(beta, removed, "condition[treated-control]")
  expect_true(is.na(truth$beta_clr_U[[1L]]))
  expect_equal(truth$beta_clr_original, as.numeric(beta))
  expect_equal(mean(truth$beta_clr_U[truth$retained]), 0, tolerance = 1e-10)
  expect_false(all(abs(truth$beta_clr_U[truth$retained] -
                        truth$beta_clr_original[truth$retained]) <= 1e-10))
})

test_that("G2 uses complete profiles and keeps their canonical taxon IDs", {
  reservoir <- benchmark_test_reservoir()
  first <- benchmark_simulate("P8", reservoir = reservoir)
  second <- benchmark_simulate("P8", reservoir = reservoir)
  expect_identical(first$counts_original, second$counts_original)
  expect_identical(rownames(first$counts_original), rownames(reservoir$counts))
  expect_true(all(first$metadata$source_profile %in% colnames(reservoir$counts)))
  expect_equal(first$manifest$reservoir_checksum, reservoir$checksum)
  expect_equal(dim(first$counts_original), c(100L, 60L))
  expect_error(benchmark_simulate("P8"), "frozen reservoir")
  reservoir$counts <- reservoir$counts[-1L, , drop = FALSE]
  expect_error(benchmark_simulate("P8", reservoir = reservoir), "exactly 100")
})

test_that("preflight data are deterministic and never pilot replicates", {
  for (cell_id in c("PF_G1", "PF_G2")) {
    data <- benchmark_simulate(cell_id, reservoir = benchmark_test_reservoir(),
                               preflight = TRUE)
    repeated <- benchmark_simulate(cell_id, reservoir = benchmark_test_reservoir(),
                                   preflight = TRUE)
    expect_identical(data$counts_original, repeated$counts_original)
    expect_identical(data$manifest$mode, "preflight")
    expect_equal(sum(data$truth$beta_latent > 0), 5L)
    expect_equal(sum(data$truth$beta_latent < 0), 5L)
    expect_equal(max(abs(data$truth$beta_latent)), log(4))
  }
})

test_that("the empirical reservoir is frozen deterministically from bundled data", {
  first <- benchmark_prepare_reservoir()
  second <- benchmark_prepare_reservoir()
  expect_identical(first$checksum, second$checksum)
  expect_equal(nrow(first$counts), 100L)
  expect_true(all(first$taxa$prevalence >= 0.1))
  expect_true(all(diff(first$taxa$mean_relative_abundance) >= 0))
  expect_equal(first$taxa$taxa_id, rownames(first$counts))
  expect_equal(rownames(first$taxonomy), rownames(first$counts))
  expect_match(first$manifest$donor_independence, "no verificable")
  data <- benchmark_simulate("PF_G2", reservoir = first, preflight = TRUE)
  expect_equal(phyloseq::taxa_names(dar:::analysis_phyloseq(data$recipe)),
               rownames(data$counts))
})
test_that("empirical metadata are plain exportable data frames", {
  directory <- withr::local_tempdir()
  reservoir <- benchmark_prepare_reservoir(output_dir = directory)
  expect_identical(class(reservoir$metadata), "data.frame")
  samples <- utils::read.csv(file.path(directory, "reservoir-samples.csv"), row.names = 1L)
  expect_equal(nrow(samples), ncol(reservoir$counts))
  expect_setequal(rownames(samples), colnames(reservoir$counts))
  expect_true(file.exists(file.path(directory, "reservoir-manifest.json")))
})
