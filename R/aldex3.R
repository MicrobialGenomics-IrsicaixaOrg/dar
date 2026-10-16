#' ALDEx3 fixed- and mixed-effects analysis
#'
#' Fit planned contrasts with the official ALDEx3 engine. Fixed models use
#' native linear regression; mixed models use exact `lme4` fits, not the
#' approximate `blmm` engine. A centralized [add_model()] is required.
#'
#' @details Supply sequence counts for inference. Taxonomic aggregation and
#'   common preprocessing precede native Dirichlet sampling. No extra
#'   prevalence filter or implicit "other" feature is added. Filtering changes
#'   the composition and the features defining CLR scale assumptions.
#'   External scale information must describe the retained community; removed
#'   counts are not reconstructed.
#'   Native numerical parameters and input values are passed to ALDEx3.
#'   Engine-specific validation and errors are left to that package. dar
#'   checks its configuration and taxon, sample and contrast identities.
#'
#'   Scale models are specified by name, not executable functions. `"clr"`
#'   and `"tss"` accept `scale_args = list(gamma = 0.5)`, the native default
#'   SD on the log2 coefficient scale. `gamma = 0` fixes the corresponding
#'   scale assumption. `"sample"` requires named `s.mu` and `s.sd` vectors
#'   of external log2-scale means and SDs, one per input sample ID.
#'   `"coefficient"` requires named `c.mu` and `c.sd` vectors, one per column
#'   of the resolved fixed-effects design. This interface supports independent
#'   SDs, not covariance matrices or arbitrary scale functions.
#'   The default is a modeling choice, not an estimate of scale uncertainty.
#'   Coefficient coding and the units of numeric covariates affect these
#'   priors; choose scale information appropriate to the study.
#'
#'   Scale uncertainty and priors always refer to the original resolved
#'   design, even when a planned contrast requires releveling. Each contrast
#'   uses the same seed and scale-model basis. The native model is refitted
#'   for each contrast; mixed models can therefore be expensive.
#'
#'   Effects are mean Monte Carlo coefficients on the log2 abundance scale
#'   conditional on the specified scale model (`aldex3_log2_abundance_coefficient`).
#'   Sequence counts alone do not measure absolute microbial load. P values
#'   and adjusted P values retain ALDEx3's native directional Monte Carlo
#'   procedure, including adjustment across taxa separately for each
#'   coefficient. They are not BH applied to averaged P values, nor a joint
#'   adjustment over all planned contrasts. Missing statistics remain missing.
#'   Execution settings and the effective package version are recorded by
#'   [get_execution()].
#'
#' @param rec An unprepared Recipe object.
#' @param max_significance Native adjusted P value cutoff.
#' @param nsample Number of Monte Carlo draws passed to ALDEx3. Default 2000.
#' @param scale_model One of `"clr"`, `"tss"`, `"sample"`, `"coefficient"`.
#' @param scale_args Named list of scale-model parameters; see Details.
#' @param p_adj_method Native multiple-testing adjustment, default `"BH"`.
#' @param n_cores Number of cores for native mixed-model fitting, default one.
#' @param seed Integer Monte Carlo seed. Restores the caller's RNG state.
#' @param rarefy Whether to rarefy before fitting; default `FALSE`.
#' @param id Unique step identifier.
#' @param engine_args Named list with an optional `fit` list. Supports native
#'   `test` (`"t.HC3"`, `"t.HC0"`, `"t"`) for fixed models and `streamsize`.
#'   Validation is left to ALDEx3. The native test setting does not apply
#'   to mixed models.
#'   Data, model, scale, engine, seed and two-sided inference are managed by dar.
#' @return A Recipe with the ALDEx3 step added.
#' @family Diff taxa steps
#' @include recipe-class.R
#' @export
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus") |>
#'   step_aldex3(scale_model = "clr", scale_args = list(gamma = 0.5))
#' # Use verified sequence counts for inference; prep() executes the recipe.
step_aldex3 <- function(rec, max_significance = 0.05, nsample = 2000,
                        scale_model = "clr", scale_args = list(gamma = 0.5),
                        p_adj_method = "BH", n_cores = 1, seed = 130L,
                        rarefy = FALSE, id = rand_id("aldex3"),
                        engine_args = list()) {
  check_recipe(rec)
  scale_model <- match.arg(scale_model, c("clr", "tss", "sample", "coefficient"))
  scale_args <- normalize_aldex3_scale_args(scale_model, scale_args)
  recipes_pkg_check("ALDEx3", "step_aldex3()")
  add_step(
    rec,
    step(
      "aldex3",
      max_significance = max_significance,
      nsample = nsample,
      scale_model = scale_model,
      scale_args = scale_args,
      p_adj_method = p_adj_method,
      n_cores = n_cores,
      seed = seed,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("aldex3", engine_args)
    )
  )
}

#' @noRd
normalize_aldex3_scale_args <- function(scale_model, scale_args) {
  allowed <- switch(
    scale_model,
    clr = "gamma",
    tss = "gamma",
    sample = c("s.mu", "s.sd"),
    coefficient = c("c.mu", "c.sd")
  )
  if (is.list(scale_args) && !length(scale_args) && length(allowed) == 1L) {
    scale_args <- list(gamma = 0.5)
  }
  if (!is.list(scale_args) || (length(scale_args) > 0L &&
      (is.null(names(scale_args)) || anyDuplicated(names(scale_args)) ||
       any(!nzchar(names(scale_args))) || anyNA(names(scale_args)) ||
       length(setdiff(names(scale_args), allowed)) > 0L))) {
    cli::cli_abort(
      paste(
        "Scale model {.val {scale_model}} supports only",
        "{.arg {allowed}} in {.arg scale_args}."
      ),
      class = "dar_error_invalid_engine_args"
    )
  }
  scale_args
}

#' @noRd
align_aldex3_scale_args <- function(scale_model, scale_args, ids, allow_extra = FALSE) {
  if (scale_model %in% c("clr", "tss")) return(scale_args)
  purrr::map(scale_args, function(value) {
    if (is.null(value)) { return(NULL) }
    value_ids <- names(value)
    if (is.null(value_ids) || anyNA(value_ids) || any(!nzchar(value_ids)) ||
        anyDuplicated(value_ids) || !all(ids %in% value_ids) ||
        (!allow_extra && !setequal(value_ids, ids))) {
      cli::cli_abort(
        "Scale vectors must have unique names covering the resolved {.val {scale_model}} IDs.",
        class = "dar_error_invalid_model"
      )
    }
    unname(value[ids])
  })
}

#' @noRd
run_aldex3 <- function(rec, max_significance, nsample, scale_model, scale_args,
                       p_adj_method, n_cores, seed, rarefy, id,
                       engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "aldex3", engine_args)
  run_aldex3_model(
    rec,
    max_significance,
    nsample,
    scale_model,
    scale_args,
    p_adj_method,
    n_cores,
    seed,
    rarefy,
    engine_args
  )
}

#' @noRd
required_pkgs_aldex3 <- function(x, ...) "ALDEx3"
