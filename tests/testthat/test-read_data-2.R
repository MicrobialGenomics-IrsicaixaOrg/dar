# File created by roxut; edit the function definition file, not this file

# Test found in read_data.R:214 (file:line)
  
df <- read_file(system.file("extdata", "metaHIV_metadata.txt", package = "dar"))

expect_s3_class(validate_sample_data(df), c("tbl_df", "tbl", "data.frame"))
expect_equal(dim(df), c(156, 4))
expect_error(validate_sample_data(dplyr::select(df, -1)))
expect_error(validate_sample_data(dplyr::mutate(
  df, sample_id = ifelse(sample_id == "Sample_169", NA_character_, sample_id)
)))
expect_error(validate_sample_data(dplyr::mutate(
  df, sample_id = ifelse(sample_id == "Sample_169", "Sample_162", sample_id)
)))
expect_error(validate_sample_data(dplyr::mutate(df, sample_id = 1:nrow(df))))