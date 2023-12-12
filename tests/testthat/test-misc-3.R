# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:114 (file:line)
  
data(test_prep_rec)
exprs <- test_prep_rec@steps |> purrr::map_chr(step_to_expr)
expect_length(exprs, 5)
expect_true(all(stringr::str_detect(exprs, "run_")))