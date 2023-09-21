# File created by roxut; edit the function definition file, not this file

# Test found in filter_taxa.R:23 (file:line)
  
data(metaHIV_phy)
rec <- 
  recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0 * length(x))") %>% 
  step_metagenomeseq(rm_zeros = 0)
  
expect_error(prep(rec))
  