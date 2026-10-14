test_that("corncob selects the requested native mean-model covariance", {
  skip_if_not_installed("corncob")
  calls <- character()
  testthat::local_mocked_bindings(
    hessian = function(fit) {
      calls <<- c(calls, "hessian")
      diag(c(4, 25, 100))
    },
    sand_vcov = function(fit) {
      calls <<- c(calls, "sand_vcov")
      diag(c(1, 4, 9))
    },
    .package = "corncob"
  )
  fit <- list(np.mu = 2L, np.phi = 1L, robust = FALSE)
  basic <- dar:::corncob_model_covariance(fit)
  expect_equal(basic, diag(c(0.25, 0.04)))
  expect_identical(calls, "hessian")

  fit$robust <- TRUE
  robust <- dar:::corncob_model_covariance(fit)
  expect_equal(robust, diag(c(1, 4)))
  expect_identical(calls, c("hessian", "sand_vcov"))
  expect_identical(dim(robust), c(2L, 2L))
})

test_that("central corncob propagates robust covariance to its decisions", {
  skip_if_not_installed("corncob")
  testthat::local_mocked_bindings(
    hessian = function(fit) diag(c(4, 25, 100)),
    sand_vcov = function(fit) diag(c(1, 4, 9)),
    .package = "corncob"
  )
  testthat::local_mocked_bindings(
    engine_stage_function = function(engine, stage) {
      expect_identical(engine, "corncob")
      expect_identical(stage, "fit")
      function(..., robust = FALSE) {
        list(np.mu = 2L, np.phi = 1L, b.mu = c(0, 1),
             robust = robust, sep_da = FALSE)
      }
    },
    .package = "dar"
  )
  rec <- recipe(make_longitudinal_phy()) |>
    add_model(~ condition, targets = "condition", tax_level = "Species")
  run <- function(engine_args = list()) {
    dar:::run_corncob_model(
      rec, ~ 1, "logit", "logit", FALSE, 0.05, "BH", 0, FALSE,
      engine_args = engine_args
    )$model
  }
  default <- run()
  basic <- run(list(fit = list(robust = FALSE)))
  robust <- run(list(fit = list(robust = TRUE)))

  expect_identical(default, basic)
  expect_equal(basic$effect, rep(1, nrow(basic)))
  expect_identical(robust$effect, basic$effect)
  expect_equal(basic$stderr, rep(0.2, nrow(basic)))
  expect_equal(robust$stderr, rep(2, nrow(robust)))
  expect_equal(basic$pval, rep(2 * stats::pnorm(-5), nrow(basic)))
  expect_equal(robust$pval, rep(2 * stats::pnorm(-0.5), nrow(robust)))
  expect_equal(basic$padj, stats::p.adjust(basic$pval, method = "BH"))
  expect_equal(robust$padj, stats::p.adjust(robust$pval, method = "BH"))
  expect_true(all(basic$signif))
  expect_false(any(robust$signif))
  expect_identical(robust$taxa_id, basic$taxa_id)
  expect_identical(robust$contrast_id, basic$contrast_id)
})

test_that("legacy corncob forwards the requested bootstrap count", {
  skip_if_not_installed("corncob")
  calls <- list()
  testthat::local_mocked_bindings(
    differentialTest = function(formula, data, phi.formula, formula_null,
                                phi.formula_null, link, phi.link, test,
                                boot, B, filter_discriminant, fdr_cutoff, fdr) {
      calls[[length(calls) + 1L]] <<- list(B = B, boot = boot)
      ids <- phyloseq::taxa_names(data)
      models <- lapply(ids, function(id) {
        list(
          np.mu = 2L,
          coefficients = rbind(
            intercept = c(Estimate = 0, Std.Error = 0.1),
            conditiontreated = c(Estimate = 1, Std.Error = 0.1)
          )
        )
      })
      structure(list(
        data = data, significant_taxa = ids, significant_models = models,
        p_fdr = stats::setNames(rep(0.01, length(ids)), ids)
      ), class = "differentialTest")
    },
    .package = "corncob"
  )
  phy <- make_longitudinal_phy()
  taxonomy <- cbind(
    Genus = paste0("Genus_", seq_len(phyloseq::ntaxa(phy))),
    as(phyloseq::tax_table(phy), "matrix")
  )
  phyloseq::tax_table(phy) <- phyloseq::tax_table(taxonomy)
  rec <- suppressWarnings(recipe(phy, "condition", "Species"))

  custom <- rec |>
    step_corncob(boot = TRUE, B = 7L, id = "corncob__bootstrap_custom")
  default <- rec |>
    step_corncob(boot = TRUE, id = "corncob__bootstrap_default")
  custom_prepared <- suppressWarnings(prep(custom, parallel = FALSE))
  default_prepared <- suppressWarnings(prep(default, parallel = FALSE))
  expect_s4_class(custom_prepared, "PrepRecipe")
  expect_s4_class(default_prepared, "PrepRecipe")
  expect_identical(calls, list(
    list(B = 7L, boot = TRUE), list(B = 1000, boot = TRUE)
  ))
})

test_that("corncob robust covariance agrees with the native implementation", {
  skip_if_not_installed("corncob")
  withr::local_seed(130)
  data <- data.frame(
    condition = factor(rep(c("control", "treated"), each = 10)),
    M = rep(200L, 20)
  )
  proportion <- stats::rbeta(20, shape1 = 4 + 2 * (data$condition == "treated"),
                            shape2 = 10)
  data$W <- stats::rbinom(20, size = data$M, prob = proportion)
  fit <- suppressWarnings(corncob::bbdml(
    cbind(W, M - W) ~ condition, phi.formula = ~ 1,
    data = data, robust = TRUE
  ))
  native <- utils::getFromNamespace("sand_vcov", "corncob")(fit)
  expect_identical(dim(native), c(fit$np.total, fit$np.total))
  index <- seq_len(fit$np.mu)
  expect_equal(
    dar:::corncob_model_covariance(fit),
    native[index, index, drop = FALSE]
  )
  fit$robust <- FALSE
  expect_equal(
    dar:::corncob_model_covariance(fit),
    solve(utils::getFromNamespace("hessian", "corncob")(fit))[
      index, index, drop = FALSE
    ]
  )
})
