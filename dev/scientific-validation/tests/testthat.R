library(testthat)

Sys.setenv(TESTTHAT_PARALLEL = "FALSE")
pkgload::load_all(".", quiet = TRUE)

test_dir(
  file.path("dev", "scientific-validation", "tests", "testthat"),
  reporter = "summary"
)
