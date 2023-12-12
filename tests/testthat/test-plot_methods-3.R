# File created by roxut; edit the function definition file, not this file

# Test found in plot_methods.R:326 (file:line)
  
data(test_rec)
data(test_prep_rec)
taxa_ids <- c("Otu_96", "Otu_78", "Otu_88", "Otu_35", "Otu_94", "Otu_34")

expect_error(abundance_plt(test_rec))
expect_s3_class(abundance_plt(test_prep_rec), "ggplot")
expect_s4_class(abundance_plt(test_prep_rec, type = "heatmap"), "HeatmapList")
expect_s3_class(abundance_plt(test_prep_rec, taxa_ids = taxa_ids), "ggplot")
expect_s3_class(abundance_plt(test_prep_rec, top_n = 10), "ggplot")
expect_s3_class(abundance_plt(test_prep_rec, transform = "clr"), "ggplot")
expect_s4_class(
  abundance_plt(test_prep_rec, taxa_ids = taxa_ids, type = "heatmap"), 
  "HeatmapList"
)
expect_s3_class(
  abundance_plt(test_prep_rec, transform = "scale", scale = 10), 
  "ggplot"
)