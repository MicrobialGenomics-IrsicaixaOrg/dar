# File created by roxut; edit the function definition file, not this file

# Test found in steps_and_checks.R:15 (file:line)
  
data(metaHIV_phy)
expect_s4_class(
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
    step_metagenomeseq() |>
    step_metagenomeseq(),
  "Recipe"
)

data(test_prep_rec)
expect_s4_class(
  test_prep_rec |> 
    bake() |>
    bake(),
 "PrepRecipe"
)