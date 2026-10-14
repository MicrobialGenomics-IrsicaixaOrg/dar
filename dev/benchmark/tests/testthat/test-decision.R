benchmark_fixture_decision <- function() {
  panel <- c("ancom_pilot", "linda_pilot", "aldex_pilot", "maaslin_pilot")
  cells <- paste0("P", 1:9)
  runs <- rbind(
    expand.grid(cell_id = cells, replicate = 1:20,
                procedure = c(panel, "majority_3of4", "unanimity_4of4"),
                stringsAsFactors = FALSE),
    expand.grid(cell_id = setdiff(cells, c("P5", "P6")), replicate = 1:20,
                procedure = "wilcox_pilot", stringsAsFactors = FALSE)
  )
  runs$status <- "success"
  summary <- merge(unique(runs[c("cell_id", "procedure")]), data.frame(
    truth = c("truth_latent", "truth_clr", "truth_median"), stringsAsFactors = FALSE
  ), by = NULL)
  summary$n_planned <- summary$n_recorded <- summary$n_valid <- 20L
  summary$is_null <- summary$cell_id %in% c("P1", "P5", "P7")
  summary$signal_retention <- 1
  summary$fdr <- summary$paired_procedure_fdp <- 0.1
  summary$directional_power <- summary$paired_procedure_directional_power <- 0.6
  summary$paired_reference_directional_power <- 0.6
  summary$delta_fdp <- summary$delta_directional_power <- 0
  summary$n_paired_fdp <- summary$n_paired_power <- 20L
  summary$null_events <- ifelse(summary$cell_id %in% c("P1", "P5", "P7"), 0L, NA_integer_)
  list(summary = summary, runs = runs, diagnostics = list(
    mode = "pilot", unexplained_discrepancies = 0L,
    hypotheses_interpretable = TRUE, failures_documented = TRUE,
    projected_cpu_hours = 2, uncertainty_too_wide = FALSE,
    discordant_taxa = data.frame(cell_id = c("P4", "P9"), mean = c(0, 0),
                                 explained = c(TRUE, TRUE))
  ))
}

testthat::test_that("preflight and partial pilots never pass technical gates", {
  fixture <- benchmark_fixture_decision()
  fixture$diagnostics$mode <- "preflight"
  result <- do.call(benchmark_decide, fixture)
  testthat::expect_equal(result$decision, "not_applicable")
  testthat::expect_false(result$eligible_for_expansion)
  fixture <- benchmark_fixture_decision()
  fixture$runs <- fixture$runs[-1L, ]
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "incomplete")
  fixture <- benchmark_fixture_decision()
  fixture$summary <- fixture$summary[-1L, ]
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "incomplete")
})

testthat::test_that("A enforces fit completeness, retention and documented cost", {
  fixture <- benchmark_fixture_decision()
  fixture$runs$status[fixture$runs$cell_id == "P2" & fixture$runs$procedure == "linda_pilot" &
                        fixture$runs$replicate <= 2] <- "engine_error"
  fixture$summary$n_valid[fixture$summary$cell_id == "P2" &
                            fixture$summary$procedure == "linda_pilot"] <- 18L
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "A_revise_technical")
  fixture <- benchmark_fixture_decision()
  fixture$summary$signal_retention[fixture$summary$cell_id == "P8"] <- 0.89
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "A_revise_technical")
  fixture <- benchmark_fixture_decision()
  fixture$diagnostics$projected_cpu_hours <- 25
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "A_revise_technical")
})

testthat::test_that("B requires the same pathway in G1 and G2", {
  fixture <- benchmark_fixture_decision()
  selected <- fixture$summary$procedure == "majority_3of4" &
    fixture$summary$cell_id %in% c("P2", "P8") & fixture$summary$truth == "truth_latent"
  fixture$summary$delta_fdp[selected] <- -0.06
  fixture$summary$paired_procedure_directional_power[selected] <- 0.5
  result <- do.call(benchmark_decide, fixture)
  testthat::expect_equal(result$decision, "B_propose_advantage_study")
  testthat::expect_false(result$fdr_control_established)
  fixture$summary$delta_fdp[fixture$summary$cell_id == "P8"] <- 0
  fixture$summary$delta_directional_power[fixture$summary$cell_id == "P8"] <- 0.11
  testthat::expect_false(do.call(benchmark_decide, fixture)$B)
})

testthat::test_that("a null inflation alert blocks the advantage pathway", {
  fixture <- benchmark_fixture_decision()
  selected <- fixture$summary$procedure == "majority_3of4" &
    fixture$summary$cell_id %in% c("P2", "P8") & fixture$summary$truth == "truth_latent"
  fixture$summary$delta_fdp[selected] <- -0.06
  fixture$summary$null_events[fixture$summary$cell_id == "P1" &
                               fixture$summary$procedure == "majority_3of4"] <- 4L
  testthat::expect_false(do.call(benchmark_decide, fixture)$B)
})

testthat::test_that("C requires a practical effect beyond algebraic truth changes", {
  fixture <- benchmark_fixture_decision()
  selected <- fixture$summary$cell_id %in% c("P4", "P9") &
    fixture$summary$procedure == "majority_3of4" & fixture$summary$truth == "truth_clr"
  fixture$summary$fdr[selected] <- 0.4
  testthat::expect_false(do.call(benchmark_decide, fixture)$C)
  fixture$diagnostics$discordant_taxa$mean <- c(5, 5)
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "C_propose_limit_study")
  fixture$diagnostics$discordant_taxa$explained[[1L]] <- FALSE
  testthat::expect_false(do.call(benchmark_decide, fixture)$C)
})

testthat::test_that("D is not automatic when power or uncertainty is uninformative", {
  fixture <- benchmark_fixture_decision()
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "D_do_not_expand")
  fixture$diagnostics$uncertainty_too_wide <- TRUE
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "inconclusive")
  fixture$diagnostics$uncertainty_too_wide <- FALSE
  fixture$summary$paired_reference_directional_power[fixture$summary$cell_id == "P8"] <- 0.19
  testthat::expect_equal(do.call(benchmark_decide, fixture)$decision, "inconclusive")
})
