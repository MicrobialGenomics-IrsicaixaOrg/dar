# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:30 (file:line)
  
data("metaHIV_phy")
phy <- metaHIV_phy
var <- "RiskGroup2"

is_list <- get_comparisons(var, phy)
is_tibb <- get_comparisons(var, phy, as_list = FALSE)
filtered <- get_comparisons(var, phy, n_cut = 30)

expect_equal(class(is_list), "list")
expect_s3_class(is_tibb, "tbl_df")
expect_length(is_list, 3)
expect_equal(nrow(is_tibb), 3)
expect_length(filtered, 1)
