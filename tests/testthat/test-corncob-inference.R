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
