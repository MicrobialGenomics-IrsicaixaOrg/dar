# File created by roxut; edit the function definition file, not this file

# Test found in lefse.R:41 (file:line)
  
data(metaHIV_phy)

test <-
 recipe(metaHIV_phy, "RiskGroup2", "Class") |>
 step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
 step_filter_by_prevalence() |> 
 step_lefse() |> 
 step_lefse(rarefy = FALSE) 
 
expect_match(
  capture_warnings(prep(test)),
  "Convert counts to relative abundances with 'relativeAb\\(\\)'", 
  all = TRUE
)  
expect_s4_class(suppressWarnings(prep(test)), "PrepRecipe")

data(test_prep_rec)
expect_error(step_lefse(test_prep_rec))