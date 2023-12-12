# File created by roxut; edit the function definition file, not this file

# Test found in plot_methods.R:123 (file:line)
  
data(test_rec)
data(test_prep_rec)
expect_error(intersection_plt(test_rec))
expect_s3_class(intersection_plt(test_prep_rec), "upset")
expect_s3_class(
  intersection_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]), 
  "upset"
)
expect_s3_class(intersection_plt(test_prep_rec, ordered_by = "degree"), "upset")