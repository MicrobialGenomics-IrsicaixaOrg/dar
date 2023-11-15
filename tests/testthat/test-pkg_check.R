# File created by roxut; edit the function definition file, not this file

# Test found in pkg_check.R:12 (file:line)
  
expected <- capture.output(
  glue::glue("1 package is needed for step_aldex() and is not installed: ", 
  "(randompackage). \n Start a clean R session then run: ", 
  "BiocManager::install(c(\"randompackage\")) ") %>% print()
)
  
res <- capture.output(
  dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()")
)
expect_equal(res, expected)
expect_invisible(
  dar:::recipes_pkg_check(dar:::required_pkgs_aldex(), "step_aldex()")
)
expect_invisible(
  dar:::recipes_pkg_check(dar:::required_pkgs_error(), "step_aldex()")
)