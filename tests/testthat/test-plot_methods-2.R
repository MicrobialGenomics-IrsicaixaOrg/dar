# File created by roxut; edit the function definition file, not this file

# Test found in plot_methods.R:207 (file:line)
  
data(test_rec)
data(test_prep_rec)
test_1 <- 
test_2 <- 
test_3 <- 

expect_error(exclusion_plt(test_rec))
expect_s3_class(exclusion_plt(test_prep_rec), "ggplot")
expect_s3_class(
  exclusion_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]),
  "ggplot"
)