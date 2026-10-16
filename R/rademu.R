#' radEmu differential-abundance analysis
#'
#' Fit planned contrasts with [radEmu::emuFit()]. A centralized [add_model()]
#' defines the fixed design and taxonomic analysis units. The default robust
#' score test is recommended by radEmu; Wald testing is available explicitly.
#'
#' @details Effects are natural-log fold differences relative to the native
#'   pseudo-Huber centre across analysis taxa. They are not unqualified
#'   absolute-abundance changes or log2 fold changes. Each planned contrast
#'   is reparameterized and fitted with its own native centring constraint.
#'   This retains the native hypothesis rather than treating the nonlinear
#'   centre as a linear transformation of previously centred coefficients.
#'
#'   Supply `cluster` to identify a metadata column containing cluster labels.
#'   Labels are aligned with the common model cohort. Clustered inference uses
#'   radEmu's GEE framework, not a mixed model. Random-effect formulas such as
#'   `(1 | subject)` are incompatible and are never silently translated. For
#'   repeated observations, use a fixed formula and explicitly set
#'   `cluster = "subject"`. Missing cluster IDs are not dropped privately.
#'
#'   dar adjusts the selected native P values across returned taxa separately
#'   for each planned contrast. It does not adjust across all contrasts jointly.
#'   Missing native statistics remain missing, including native zero-comparison
#'   exclusions. Native convergence warnings are retained; diagnostics and
#'   effective settings are available through [get_execution()].
#'
#'   The complete taxonomic analysis matrix is used in every fit. Counts are
#'   not split by taxon or restricted to the two groups named in a contrast.
#'   Score tests can be expensive because each taxon requires a restricted fit.
#'   Native numerical validation and errors are left to radEmu.
#'
#' @param rec An unprepared [Recipe-class] object.
#' @param test Native test to report: `"score"` (default) or `"wald"`.
#' @param cluster Optional single metadata column identifying dependent
#'   observations. `NULL` uses native independent-observation inference.
#' @param max_significance Adjusted P value cutoff.
#' @param p_adj_method Multiple-testing adjustment passed to [p.adjust()].
#' @param constraint_param Native pseudo-Huber smoothing parameter, default 0.1.
#' @param seed Integer seed used for native fitting. Restores the caller's RNG.
#' @param rarefy Whether to rarefy before fitting; default `FALSE`.
#' @param id Unique step identifier.
#' @param engine_args Named list with an optional `fit` list of native
#'   arguments, for example `tolerance`, `maxit`, `maxit_null`, `alpha` and
#'   `remove_zero_comparison_pvals`. Solver controls can be passed directly in
#'   `fit`. Data, design, clustering, centring functions, fitted-model reuse,
#'   tested coordinates, test selection and result shape are managed by dar;
#'   `control` and `return_both_score_pvals` cannot be overridden.
#' @return A Recipe with the radEmu step added.
#' @family Diff taxa steps
#' @include recipe-class.R
#' @export
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus") |>
#'   step_rademu()
#' # prep(rec) executes native score tests for every planned contrast and taxon.
step_rademu <- function(rec, test = c("score", "wald"), cluster = NULL,
                        max_significance = 0.05, p_adj_method = "BH",
                        constraint_param = 0.1, seed = 130L, rarefy = FALSE,
                        id = rand_id("rademu"), engine_args = list()) {
  check_recipe(rec)
  test <- match.arg(test)
  checkmate::assert_string(cluster, null.ok = TRUE, min.chars = 1)
  recipes_pkg_check("radEmu", "step_rademu()")
  add_step(
    rec,
    step(
      "rademu",
      test = test,
      cluster = cluster,
      max_significance = max_significance,
      p_adj_method = p_adj_method,
      constraint_param = constraint_param,
      seed = seed,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("rademu", engine_args)
    )
  )
}

#' @noRd
run_rademu <- function(rec, test, cluster, max_significance, p_adj_method,
                       constraint_param, seed, rarefy, id, engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "rademu", engine_args)
  run_rademu_model(
    rec,
    test,
    cluster,
    max_significance,
    p_adj_method,
    constraint_param,
    seed,
    rarefy,
    engine_args
  )
}

#' @noRd
required_pkgs_rademu <- function(x, ...) "radEmu"
