# File created by roxut; edit the function definition file, not this file

# Test found in read_data.R:257 (file:line)
  
df <- read_file(system.file("extdata", "metaHIV_taxas.txt", package = "dar"))

expect_s3_class(validate_tax_table(df), c("tbl_df", "tbl", "data.frame"))
expect_equal(dim(df), c(451, 8))
expect_error(validate_tax_table(dplyr::select(df, -1)))
expect_error(validate_tax_table(dplyr::mutate(
  df, otu_id = ifelse(otu_id == "Otu_1", NA_character_, otu_id)
)))
expect_error(validate_tax_table(dplyr::mutate(
  df, otu_id = ifelse(otu_id == "Otu_1", "Otu_2", otu_id)
)))
expect_error(validate_tax_table(dplyr::mutate(df, dplyr::across(2, as.factor))))
expect_error(validate_tax_table(dplyr::mutate(df, otu_id = 1:nrow(df))))