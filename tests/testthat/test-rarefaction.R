# File created by roxut; edit the function definition file, not this file

# Test found in rarefaction.R:34 (file:line)
  
data(metaHIV_phy)
data(test_prep_rec)
test <- 
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Class") |>
  step_rarefaction() 
  
expect_s4_class(prep(test), "PrepRecipe")
expect_error(step_rarefaction(test_prep_rec))