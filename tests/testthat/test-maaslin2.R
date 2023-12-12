# File created by roxut; edit the function definition file, not this file

# Test found in maaslin2.R:45 (file:line)
  
data(metaHIV_phy)

test <-
 recipe(metaHIV_phy, "RiskGroup2", "Class") |>
 step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
 step_filter_by_abundance() |> 
 step_maaslin() |> 
 step_maaslin(rarefy = TRUE) |> 
 step_maaslin(rarefy = "no_seed") 
 
expect_s4_class(prep(test), "PrepRecipe")

data(test_prep_rec)
expect_error(step_maaslin(test_prep_rec))