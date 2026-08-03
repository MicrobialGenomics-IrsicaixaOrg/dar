#' LinDA differential-abundance analysis
#'
#' `step_linda()` fits bias-corrected linear or linear mixed-effects models with
#' [MicrobiomeStat::linda()]. It requires a centralized model defined by
#' [add_model()] and executes every contrast in that model's contrast plan.
#'
#' @param rec A [Recipe-class] object with a centralized model.
#' @param prev_filter Minimum feature prevalence retained for analysis.
#' @param mean_abund_filter Minimum mean relative abundance retained.
#' @param max_abund_filter Minimum maximum relative abundance retained.
#' @param winsorize Whether high-abundance outliers are winsorized.
#' @param outlier_pct Expected fraction of outliers used for winsorization.
#' @param adaptive Whether LinDA chooses zero handling from depth-covariate
#'   correlations.
#' @param zero_handling Zero treatment used when `adaptive = FALSE`: either
#'   `"pseudo-count"` or `"imputation"`.
#' @param pseudo_count Positive pseudo-count used by pseudo-count zero handling.
#' @param corr_cut Significance cutoff used by adaptive zero handling.
#' @param p_adj_method Multiple-testing correction accepted by [p.adjust()].
#' @param alpha Adjusted p-value threshold used to classify significance.
#' @param n_cpus Number of cores used by LinDA mixed-effects models.
#' @param rarefy Whether counts are rarefied for this step.
#' @param id Unique identifier for this configured step.
#' @param engine_args Named lists of advanced arguments for the native `fit`
#'   stage. Arguments managed by dar or exposed above cannot be overridden.
#'
#' @return A [Recipe-class] object.
#' @family Diff taxa steps
#' @aliases step_linda
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' expect_s3_class(step_linda(rec)@steps[[1]], "step_linda")
#' data(test_prep_rec)
#' expect_error(step_linda(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(
#'     ~ RiskGroup2,
#'     targets = "RiskGroup2",
#'     tax_level = "Species"
#'   ) |>
#'   step_linda(prev_filter = 0.1, adaptive = TRUE)
#'
#' rec
step_linda <- function(rec,
                       prev_filter = 0,
                       mean_abund_filter = 0,
                       max_abund_filter = 0,
                       winsorize = TRUE,
                       outlier_pct = 0.03,
                       adaptive = TRUE,
                       zero_handling = c("pseudo-count", "imputation"),
                       pseudo_count = 0.5,
                       corr_cut = 0.1,
                       p_adj_method = "BH",
                       alpha = 0.05,
                       n_cpus = 1L,
                       rarefy = FALSE,
                       id = rand_id("linda"),
                       engine_args = list()) {
  check_recipe(rec)
  if (is.null(get_model(rec))) {
    cli::cli_abort(
      c(
        "x" = "{.fun step_linda} requires a centralized statistical model.",
        "i" = "Define the design, targets and taxonomic level with {.fun add_model}."
      ),
      class = "dar_error_model_required"
    )
  }
  recipes_pkg_check(required_pkgs_linda(), "step_linda()")

  checkmate::assert_number(prev_filter, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(mean_abund_filter, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_number(max_abund_filter, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_flag(winsorize)
  checkmate::assert_number(outlier_pct, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_flag(adaptive)
  zero_handling <- match.arg(zero_handling)
  checkmate::assert_number(pseudo_count, lower = 0, finite = TRUE)
  if (pseudo_count == 0) {
    cli::cli_abort("{.arg pseudo_count} must be positive.")
  }
  checkmate::assert_number(corr_cut, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_choice(p_adj_method, stats::p.adjust.methods)
  checkmate::assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_count(n_cpus, positive = TRUE)

  add_step(
    rec,
    step(
      subclass = "linda",
      prev_filter = prev_filter,
      mean_abund_filter = mean_abund_filter,
      max_abund_filter = max_abund_filter,
      winsorize = winsorize,
      outlier_pct = outlier_pct,
      adaptive = adaptive,
      zero_handling = zero_handling,
      pseudo_count = pseudo_count,
      corr_cut = corr_cut,
      p_adj_method = p_adj_method,
      alpha = alpha,
      n_cpus = as.integer(n_cpus),
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("linda", engine_args)
    )
  )
}

#' @noRd
#' @autoglobal
run_linda <- function(rec,
                      prev_filter,
                      mean_abund_filter,
                      max_abund_filter,
                      winsorize,
                      outlier_pct,
                      adaptive,
                      zero_handling,
                      pseudo_count,
                      corr_cut,
                      p_adj_method,
                      alpha,
                      n_cpus,
                      rarefy,
                      id,
                      engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "linda", engine_args, id)
  if (is.null(get_model(rec))) {
    cli::cli_abort(
      "{.fun step_linda} cannot execute without {.fun add_model}.",
      class = "dar_error_model_required"
    )
  }
  run_linda_model(
    rec = rec,
    prev_filter = prev_filter,
    mean_abund_filter = mean_abund_filter,
    max_abund_filter = max_abund_filter,
    winsorize = winsorize,
    outlier_pct = outlier_pct,
    adaptive = adaptive,
    zero_handling = zero_handling,
    pseudo_count = pseudo_count,
    corr_cut = corr_cut,
    p_adj_method = p_adj_method,
    alpha = alpha,
    n_cpus = n_cpus,
    rarefy = rarefy,
    engine_args = engine_args
  )
}

#' @noRd
required_pkgs_linda <- function(x, ...) {
  "MicrobiomeStat"
}
