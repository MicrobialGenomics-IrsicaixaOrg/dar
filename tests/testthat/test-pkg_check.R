library(dar)

test_that("testrecipes_pkg_check", {
  expected <-
    "1 package is needed for step_aldex() and is not installed: (randompackage). Start a clean R session then run: pak::pkg_install(c(\"bioc::randompackage\")) "
  res <- capture.output(dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()"))
  expect_equal(res, expected)
  expect_invisible(dar:::recipes_pkg_check(dar:::required_pkgs_aldex(), "step_aldex()"))
  expect_invisible(dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()"))
})
