test_that("optional plotting dependencies fail with a classified error", {
  expect_error(
    dar:::require_optional_package("darPackageThatCannotExist", "example_plot"),
    class = "dar_error_missing_package"
  )
})
