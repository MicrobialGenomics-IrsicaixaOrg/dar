# File created by roxut; edit the function definition file, not this file

# Test found in plot_methods.R:16 (file:line)
  
data(test_rec)
data(test_prep_rec)
expect_error(corr_heatmap(test_rec))
expect_s3_class(corr_heatmap(test_prep_rec), "plotly")
expect_s3_class(
  corr_heatmap(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]), 
  "plotly"
)
expect_s3_class(corr_heatmap(test_prep_rec, type = "da"), "plotly")