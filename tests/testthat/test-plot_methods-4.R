# File created by roxut; edit the function definition file, not this file

# Test found in plot_methods.R:621 (file:line)
  
data(test_rec)
data(test_prep_rec)

expect_error(mutual_plt(test_rec))
expect_s3_class(mutual_plt(test_prep_rec), "ggplot")
expect_s3_class(mutual_plt(test_prep_rec, count_cutoff = 2), "ggplot")
expect_error(mutual_plt(test_prep_rec, count_cutoff = 0))
expect_error(mutual_plt(test_prep_rec, count_cutoff = 10))
expect_s3_class(mutual_plt(test_prep_rec, comparisons = "hts_msm"), "ggplot")
expect_error(mutual_plt(test_prep_rec, comparisons = "hts_sm"))
expect_s3_class(
  mutual_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]), 
  "ggplot"
)
expect_error(mutual_plt(test_prep_rec, steps = "da"))
expect_s3_class(mutual_plt(test_prep_rec, top_n = 10), "ggplot")
expect_error(mutual_plt(test_prep_rec, top_n = 0))