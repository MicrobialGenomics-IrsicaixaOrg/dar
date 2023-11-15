# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:80 (file:line)
  
data("metaHIV_phy")
phy <- metaHIV_phy
df <- phyloseq::otu_table(phy)

test_1 <- to_tibble(df)
test_2 <- to_tibble(df, id_name = "random_id")

expect_length(test_1, 157)
expect_s3_class(test_1, "tbl_df")
expect_s3_class(test_2, "tbl_df")
expect_equal(names(test_1)[1], "otu_id")
expect_equal(names(test_2)[1], "random_id")
