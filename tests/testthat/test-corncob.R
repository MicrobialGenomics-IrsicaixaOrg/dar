# File created by roxut; edit the function definition file, not this file

# Test found in corncob.R:49 (file:line)
  
data(test_prep_rec)
expect_error(step_corncob(test_prep_rec, rarefy = TRUE))
