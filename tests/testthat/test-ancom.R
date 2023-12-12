# File created by roxut; edit the function definition file, not this file

# Test found in ancom.R:75 (file:line)
  
data(metaHIV_phy)

test <-
 recipe(metaHIV_phy, "RiskGroup2", "Genus") |>
 step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
 step_filter_by_prevalence() |> 
 step_ancom()
 
expect_s4_class(prep(test), "PrepRecipe")

data(test_prep_rec)
expect_error(step_ancom(test_prep_rec, rarefy = TRUE))