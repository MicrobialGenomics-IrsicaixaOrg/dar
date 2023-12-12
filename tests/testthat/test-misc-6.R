# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:366 (file:line)
  
data(metaHIV_phy)
data(test_prep_rec)
expect_snapshot_file(export_steps(test_prep_rec, "test.json"), "test.json")
expect_snapshot_output(
  recipe(metaHIV_phy, "RiskGroup2", "Species") |> 
    import_steps(system.file("extdata", "test_bake.json", package = "dar"))
)