# File created by roxut; edit the function definition file, not this file

# Test found in read_data.R:310 (file:line)
  
## phyloseq validation
phy <-
  system.file("extdata", "metaHIV_phy.rds", package = "dar") %>%
  read_phyloseq()

no_tax_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::sample_data(phy))
no_sam_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::tax_table(phy))
only_otu_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy))

expect_error(validate_phyloseq(
  read_file(system.file("extdata", "metaHIV_taxas.txt", package = "dar")), 
  slots = ""
))
expect_s4_class(validate_phyloseq(phy), "phyloseq")
expect_s4_class(validate_phyloseq(no_sam_phy, "tax_table"), "phyloseq")
expect_s4_class(validate_phyloseq(no_tax_phy, "sample_data"), "phyloseq")
expect_s4_class(validate_phyloseq(only_otu_phy, ""), "otu_table")
expect_error(validate_phyloseq(no_tax_phy))
expect_error(validate_phyloseq(no_sam_phy))
expect_error(validate_phyloseq(only_otu_phy))
expect_error(validate_phyloseq(no_tax_phy, "tax_table"))
expect_error(validate_phyloseq(no_sam_phy, "sample_data"))