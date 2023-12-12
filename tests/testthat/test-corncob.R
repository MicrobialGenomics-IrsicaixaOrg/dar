# File created by roxut; edit the function definition file, not this file

# Test found in corncob.R:50 (file:line)
  
data(metaHIV_phy)

test <-
 recipe(metaHIV_phy, "RiskGroup2", "Class") |>
 step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
 step_filter_by_abundance() |> 
 step_corncob() |> 
 prep()
 
expect_s4_class(test, "PrepRecipe")

data(test_prep_rec)
expect_error(step_corncob(test_prep_rec))