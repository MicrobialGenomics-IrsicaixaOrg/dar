# File created by roxut; edit the function definition file, not this file

# Test found in aldex2.R:61 (file:line)
  
data(metaHIV_phy)

test <-
 recipe(metaHIV_phy, "RiskGroup2", "Class") |>
 step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
 step_filter_by_abundance() |> 
 step_aldex() |> 
 prep()
 
expect_s4_class(test, "PrepRecipe")

data(test_prep_rec)
expect_error(step_aldex(test_prep_rec))