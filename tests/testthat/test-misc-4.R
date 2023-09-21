# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:159 (file:line)
  
data(test_prep_rec)
res_1 <- find_intersections(test_prep_rec, steps = steps_ids(test_prep_rec, type = "da"))
res_2 <- find_intersections(test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-1])
res_3 <- find_intersections(test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-2])

expect_equal(nrow(res_1), 88)
expect_equal(nrow(res_2), 36)
expect_equal(nrow(res_3), 86)
expect_s3_class(res_1, "tbl_df")
