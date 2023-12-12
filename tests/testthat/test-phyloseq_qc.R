# File created by roxut; edit the function definition file, not this file

# Test found in phyloseq_qc.R:23 (file:line)
  
data(test_rec)
test <- phy_qc(test_rec)

## Check columns
expected_cols <- c(
  "var_levels", "n", "n_zero", "pct_zero", "pct_all_zero", "pct_singletons", 
  "pct_doubletons", "count_mean", "count_min", "count_max"
)
expect_equal(colnames(test), expected_cols)

## Check column types
expect_true(all(sapply(test[,1], is.character)))
expect_true(all(sapply(test[,-1], is.numeric)))
