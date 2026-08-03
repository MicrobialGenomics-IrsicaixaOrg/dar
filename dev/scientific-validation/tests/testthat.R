library(testthat)

pkgload::load_all(".", quiet = TRUE)

test_dir(
  file.path("dev", "scientific-validation", "tests", "testthat"),
  reporter = "summary"
)
