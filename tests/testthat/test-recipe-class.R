# File created by roxut; edit the function definition file, not this file

# Test found in recipe-class.R:71 (file:line)
   
data(metaHIV_phy) 
colnames(metaHIV_phy@tax_table) <-
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Sp")
  
expect_error( 
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species") 
)

data(GlobalPatterns, package = "mia")
expect_s4_class(recipe(GlobalPatterns), "Recipe")
