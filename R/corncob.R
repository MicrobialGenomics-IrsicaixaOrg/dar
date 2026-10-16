#' corncob analysis
#'
#' Corncob an individual taxon regression model that uses abundance tables and
#' sample data. corncob is able to model differential abundance and differential
#' variability, and addresses each of the challenges presented below:
#' 
#' * different sequencing depth
#' * excessive zeros from unobserved taxa
#' * high variability of empirical relative abundances (overdispersion)
#' * within-taxon correlation
#' * hypothesis testing with categorical and continuous covariates
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe. Execution requires [add_model()] before
#'   [prep()]; the step can be configured before adding the model.
#' @param phi.formula An object of class formula without the response: a
#'   symbolic description of the model to be fitted to the dispersion.
#' @param formula_null Retained design argument, default `NULL`. Explicit
#'   values are overridden with a warning by the [add_model()] contrasts.
#' @param phi.formula_null Formula for overdispersion under null, without
#'   response. Default `NULL`. Explicit values are overridden with a warning
#'   by the [add_model()] contrasts.
#' @param link Link function for abundance covariates, defaults to "logit".
#' @param phi.link Link function for dispersion covariates, defaults to "logit".
#' @param test Hypothesis test, default `"Wald"`. Only Wald contrasts are
#'   supported; other tests make the step incompatible with the shared model.
#' @param boot Whether to request bootstrap inference, default `FALSE`.
#'   Bootstrap is unsupported and makes the step incompatible.
#' @param B Retained bootstrap iteration setting. Not used for Wald contrasts.
#' @param filter_discriminant Boolean. Defaults to TRUE. If FALSE, discriminant
#'   taxa will not be filtered out.
#' @param fdr_cutoff Integer. Defaults to 0.05. Desired type 1 error rate.
#' @param fdr Character. Defaults to "fdr". False discovery rate control method,
#'   see p.adjust for more options.
#' @param log2FC Minimum absolute effect on the reported coefficient scale.
#'   Despite the argument name, a logit-link coefficient is not a log2 fold
#'   change. Default zero.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @inheritParams step_wilcox
#' @param engine_args Named lists of advanced arguments for the native `fit`
#'   stage. Arguments managed by dar or exposed above cannot be overridden.
#' @details Centralized models use non-bootstrap Wald contrasts. With
#'   `engine_args = list(fit = list(robust = TRUE))`, inference uses corncob's
#'   native sandwich covariance rather than the inverse Hessian. This does not
#'   add support for random effects or repeated-measures dependence.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_corncob
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @md
#' @tests 
#' data(metaHIV_phy)
#' 
#' test <-
#'  recipe(metaHIV_phy) |>
#'  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_abundance() |> 
#'  step_corncob()
#' test <- suppressWarnings(prep(test, parallel = FALSE))
#'  
#' expect_s4_class(test, "PrepRecipe")
#' 
#' data(test_prep_rec)
#' expect_error(step_corncob(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")
#'
#' rec_default <- step_corncob(rec)
#'
#' # Use a stricter significance threshold.
#' rec_strict <- step_corncob(rec, fdr_cutoff = 0.01)
step_corncob <- function(rec,
                         phi.formula = stats::formula(~ 1),
                         formula_null = NULL,
                         phi.formula_null = NULL,
                         link = "logit",
                         phi.link = "logit",
                         test = "Wald",
                         boot = FALSE,
                         B = 1000,
                         filter_discriminant = TRUE,
                         fdr_cutoff = 0.05,
                         fdr = "fdr",
                         log2FC = 0,
                         rarefy = FALSE,
                         id = rand_id("corncob"),
                         engine_args = list()) {
  check_recipe(rec)
  recipes_pkg_check(required_pkgs_corncob(), "step_corncob()")

  add_step(
    rec,
    step(
      subclass = "corncob",
      phi.formula = phi.formula,
      formula_null = formula_null,
      phi.formula_null = phi.formula_null,
      link = link,
      phi.link = phi.link,
      test = test,
      boot = boot,
      B = B,
      filter_discriminant = filter_discriminant,
      fdr_cutoff = fdr_cutoff,
      fdr = fdr,
      log2FC = log2FC,
      rarefy = rarefy,
      id = id,
      engine_args = normalize_engine_args("corncob", engine_args)
    )
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
run_corncob <- function(rec,
                        phi.formula,
                        formula_null,
                        phi.formula_null,
                        link,
                        phi.link,
                        test,
                        boot,
                        B,
                        filter_discriminant,
                        fdr_cutoff,
                        fdr,
                        log2FC,
                        rarefy,
                        id,
                        engine_args = list()) {
  engine_args <- check_engine_args_execution(rec, "corncob", engine_args)

  run_corncob_model(
    rec = rec,
    phi.formula = phi.formula,
    link = link,
    phi.link = phi.link,
    filter_discriminant = filter_discriminant,
    fdr_cutoff = fdr_cutoff,
    fdr = fdr,
    log2FC = log2FC,
    rarefy = rarefy,
    engine_args = engine_args
  )
}

#' Extract mean-model covariance using the requested native inference
#' @noRd
corncob_model_covariance <- function(fit) {
  covariance <- if (isTRUE(fit$robust)) {
    utils::getFromNamespace("sand_vcov", "corncob")(fit)
  } else {
    solve(utils::getFromNamespace("hessian", "corncob")(fit))
  }
  mean_coefficients <- seq_len(fit$np.mu)
  covariance[mean_coefficients, mean_coefficients, drop = FALSE]
}

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_corncob <- function(x, ...) { c("corncob") }
