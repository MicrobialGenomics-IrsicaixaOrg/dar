# File created by roxut; edit the function definition file, not this file

# Test found in read_data.R:127 (file:line)
  
df <- read_file(system.file("extdata", "metaHIV_counts.txt", package = "dar"))

expect_invisible(validate_otu(df))
expect_error(validate_otu(dplyr::select(df, -1)))
expect_error(validate_otu(dplyr::mutate(
  df, otu_id = ifelse(otu_id == "Otu_1", NA_character_, otu_id)
)))
expect_error(validate_otu(dplyr::mutate(
  df, otu_id = ifelse(otu_id == "Otu_1", "Otu_2", otu_id)
)))
expect_error(validate_otu(dplyr::mutate(df, dplyr::across(2, as.character))))
expect_error(validate_otu(dplyr::mutate(df, dplyr::across(
 2, ~ ifelse(.x == 0, NA_integer_, .x)
))))
expect_error(validate_otu(dplyr::mutate(df, dplyr::across(2, ~ .x * -1))))