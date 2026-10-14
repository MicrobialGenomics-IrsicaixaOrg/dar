# This decision is for a completed pilot, never for the preflight. It describes
# practical continuation criteria, not an FDR guarantee or a significance test.
benchmark_decide <- function(summary, runs, diagnostics) {
  answer <- function(decision, reasons, A = FALSE, B = FALSE, C = FALSE) {
    list(decision = decision, eligible_for_expansion = B || C, A = A, B = B,
         C = C, reasons = reasons, provisional = TRUE,
         fdr_control_established = FALSE)
  }
  if (!is.list(diagnostics) || !identical(diagnostics$mode, "pilot") ||
      ("mode" %in% names(runs) && any(runs$mode != "pilot"))) {
    return(answer("not_applicable", "Preflight is not a scientific pilot decision."))
  }
  cells <- paste0("P", 1:9)
  panel <- c("ancom_pilot", "linda_pilot", "aldex_pilot", "maaslin_pilot")
  consensus <- c("majority_3of4", "unanimity_4of4")
  compatible_wilcox <- setdiff(cells, c("P5", "P6"))
  expected <- expand.grid(cell_id = cells, replicate = 1:20,
                          procedure = c(panel, consensus), stringsAsFactors = FALSE)
  expected <- rbind(expected, expand.grid(
    cell_id = compatible_wilcox, replicate = 1:20,
    procedure = "wilcox_pilot", stringsAsFactors = FALSE
  ))
  run_keys <- c("cell_id", "replicate", "procedure")
  required_summary <- c(
    "cell_id", "procedure", "truth", "n_planned", "n_recorded", "n_valid",
    "signal_retention", "fdr", "directional_power", "null_events", "is_null",
    "delta_fdp", "delta_directional_power", "n_paired_fdp", "n_paired_power",
    "paired_reference_directional_power", "paired_procedure_directional_power",
    "paired_procedure_fdp"
  )
  if (!is.data.frame(runs) || !all(c(run_keys, "status") %in% names(runs)) ||
      !is.data.frame(summary) || !all(required_summary %in% names(summary)) ||
      anyNA(runs[c(run_keys, "status")]) || anyDuplicated(runs[run_keys]) ||
      anyDuplicated(summary[c("cell_id", "procedure", "truth")])) {
    return(answer("incomplete", "Missing or duplicated execution/summary records."))
  }
  key <- function(data, columns) do.call(paste, c(data[columns], sep = "\r"))
  if (nrow(runs) != nrow(expected) ||
      !setequal(key(runs, run_keys), key(expected, run_keys)) ||
      any(!runs$status %in% c("success", "engine_error", "timeout", "contract_invalid", "unavailable"))) {
    return(answer("incomplete", "All planned 9 x 20 dataset/procedure records are required."))
  }
  summary_expected <- unique(expected[c("cell_id", "procedure")])
  summary_expected <- merge(summary_expected, data.frame(
    truth = c("truth_latent", "truth_clr", "truth_median"), stringsAsFactors = FALSE
  ), by = NULL)
  summary_keys <- c("cell_id", "procedure", "truth")
  if (nrow(summary) != nrow(summary_expected) ||
      !setequal(key(summary, summary_keys), key(summary_expected, summary_keys)) ||
      anyNA(summary[c(summary_keys, "n_planned", "n_recorded", "n_valid")]) ||
      any(summary$n_planned != 20L | summary$n_recorded != 20L)) {
    return(answer("incomplete", "Each planned procedure requires all three truth summaries."))
  }
  numeric_counts <- c("n_planned", "n_recorded", "n_valid", "n_paired_fdp", "n_paired_power")
  bounded <- c("signal_retention", "fdr", "directional_power",
               "paired_reference_directional_power", "paired_procedure_directional_power",
               "paired_procedure_fdp")
  if (!is.logical(summary$is_null) || anyNA(summary$is_null) ||
      any(vapply(summary[numeric_counts], function(x) {
        !is.numeric(x) || any(!is.finite(x) | x != floor(x) | x < 0 | x > 20)
      }, logical(1))) ||
      any(vapply(summary[bounded], function(x) {
        !is.numeric(x) || any(is.finite(x) & (x < 0 | x > 1))
      }, logical(1))) ||
      any(summary$n_valid > 0 & !is.finite(summary$fdr)) ||
      any(summary$n_valid > 0 & !summary$is_null & !is.finite(summary$directional_power)) ||
      any(summary$n_paired_fdp > 0 & !is.finite(summary$delta_fdp)) ||
      any(summary$n_paired_power > 0 & !is.finite(summary$delta_directional_power))) {
    return(answer("incomplete", "Incomplete or invalid scientific metric denominators/estimates."))
  }
  reasons <- character()
  native <- runs[runs$procedure %in% c(panel, "wilcox_pilot"), ]
  completion <- aggregate(as.integer(native$status == "success"),
                          native[c("cell_id", "procedure")], sum)
  if (any(completion$x < 19L)) reasons <- c(reasons, "Fewer than 19/20 successful native fits in a cell.")
  for (cell in cells) {
    data <- runs[runs$cell_id == cell & runs$procedure %in% panel, ]
    complete <- vapply(split(data$status, data$replicate), function(x) all(x == "success"), logical(1))
    if (sum(complete) < 18L) reasons <- c(reasons, paste0(cell, ": fewer than 18/20 complete panels."))
  }
  signal_cells <- c("P2", "P3", "P4", "P6", "P8", "P9")
  retention <- summary[summary$procedure == "linda_pilot" & summary$truth == "truth_latent" &
                         summary$cell_id %in% signal_cells, ]
  if (any(!is.finite(retention$signal_retention)) || any(retention$signal_retention < 0.9)) {
    reasons <- c(reasons, "Original signal retention is below 90% or undocumented.")
  }
  for (i in seq_len(nrow(summary))) {
    matching <- runs$cell_id == summary$cell_id[[i]] & runs$procedure == summary$procedure[[i]]
    if (summary$n_valid[[i]] != sum(runs$status[matching] == "success")) {
      reasons <- c(reasons, "Summary availability does not match execution statuses.")
      break
    }
  }
  if (!isTRUE(diagnostics$unexplained_discrepancies == 0L)) {
    reasons <- c(reasons, "Unexplained ID, orientation or consensus discrepancies remain.")
  }
  if (!isTRUE(diagnostics$hypotheses_interpretable)) {
    reasons <- c(reasons, "Common hypotheses/universes are not interpretable.")
  }
  if (!isTRUE(diagnostics$failures_documented)) reasons <- c(reasons, "Failure/absence causes are not documented.")
  if (length(diagnostics$projected_cpu_hours) != 1L ||
      !is.finite(diagnostics$projected_cpu_hours) || diagnostics$projected_cpu_hours > 24) {
    reasons <- c(reasons, "Projected CPU cost exceeds 24 hours or is undocumented.")
  }
  if (length(reasons)) return(answer("A_revise_technical", unique(reasons)))
  get_row <- function(procedure, cell, truth = "truth_latent") {
    summary[summary$procedure == procedure & summary$cell_id == cell & summary$truth == truth, ]
  }
  primary <- c("P2", "P6", "P8")
  majority <- do.call(rbind, lapply(primary, function(cell) get_row(consensus[[1L]], cell)))
  informative <- is.finite(majority$paired_reference_directional_power) &
    majority$paired_reference_directional_power >= 0.2 &
    majority$n_paired_fdp >= 18L & majority$n_paired_power >= 18L
  pathway_1 <- informative & majority$delta_fdp <= -0.05 &
    majority$paired_procedure_directional_power >= 0.8 * majority$paired_reference_directional_power
  pathway_2 <- informative & majority$delta_directional_power >= 0.1 &
    majority$delta_fdp <= 0.02 & majority$paired_procedure_fdp <= 0.1
  repeats <- function(pattern) isTRUE(pattern[[3L]]) && any(pattern[1:2] %in% TRUE)
  harm <- any(majority$delta_directional_power < -0.2 |
                majority$delta_fdp > 0.05, na.rm = TRUE)
  null_alerts <- vapply(c("P1", "P5", "P7"), function(cell) {
    events <- get_row(consensus[[1L]], cell)$null_events
    is.finite(events) && events >= 4L
  }, logical(1))
  B <- (repeats(pathway_1) || repeats(pathway_2)) && !harm && !any(null_alerts)
  if (B) return(answer("B_propose_advantage_study", "Same majority pathway in G2 and G1, without specified harm/null alerts.", A = TRUE, B = TRUE))
  # A truth-dependent change must concern the same metric in both generators.
  truth_change <- lapply(c("P4", "P9"), function(cell) {
    latent <- get_row(consensus[[1L]], cell)
    clr <- get_row(consensus[[1L]], cell, "truth_clr")
    c(fdp = abs(latent$fdr - clr$fdr) >= 0.2,
      directional_power = abs(latent$directional_power - clr$directional_power) >= 0.2)
  })
  discordance <- diagnostics$discordant_taxa
  disagreement_explained <- is.data.frame(discordance) &&
    all(c("cell_id", "mean", "explained") %in% names(discordance)) &&
    all(vapply(c("P4", "P9"), function(cell) {
      data <- discordance[discordance$cell_id == cell, ]
      nrow(data) == 1L && isTRUE(data$mean >= 5) && isTRUE(data$explained)
    }, logical(1)))
  estimand_limit <- any((truth_change[[1L]] & truth_change[[2L]]) %in% TRUE) && disagreement_explained
  unanimity_loss <- vapply(primary, function(cell) {
    row <- get_row(consensus[[2L]], cell)
    isTRUE(row$delta_directional_power <= -0.3) && isTRUE(row$delta_fdp > -0.05)
  }, logical(1))
  C <- estimand_limit || repeats(unanimity_loss) || (null_alerts[[1L]] && null_alerts[[3L]])
  if (C) return(answer("C_propose_limit_study", "A prespecified consensus limitation warrants a separately approved focused study.", A = TRUE, C = TRUE))
  if (!all(informative)) {
    return(answer("inconclusive", "Some primary cells lack informative LinDA power or paired results; any observability amendment must be prospective.", A = TRUE))
  }
  if (!identical(diagnostics$uncertainty_too_wide, FALSE)) {
    return(answer("inconclusive", "Uncertainty has not been judged narrow enough for a no-expansion conclusion.", A = TRUE))
  }
  answer("D_do_not_expand", "Technical viability passed without a prespecified useful advantage or limitation; prioritize the software article.", A = TRUE)
}
