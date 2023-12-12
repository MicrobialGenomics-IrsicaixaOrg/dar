# File created by roxut; edit the function definition file, not this file

# Test found in read_data.R:16 (file:line)
  

## read phyloseq from rds file
good_file <- system.file("extdata", "metaHIV_phy.rds", package = "dar")
bad_file <- system.file("extdata", "invented_file.rds", package = "dar")
expect_s4_class(read_data(good_file), "phyloseq")
expect_error(read_data(bad_file))

## read phyloseq from multiple files
good_files <- c(
  system.file("extdata", "metaHIV_counts.txt", package = "dar"),
  system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
  system.file("extdata", "metaHIV_taxas.txt", package = "dar")
)
expect_s4_class(read_data(good_files), "phyloseq")

only_two_files <- c(
  system.file("extdata", "metaHIV_counts.txt", package = "dar"),
  system.file("extdata", "metaHIV_metadata.txt", package = "dar")
)
expect_error(read_data(only_two_files))

duplicated_files <- c(
  system.file("extdata", "metaHIV_counts.txt", package = "dar"),
  system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
  system.file("extdata", "metaHIV_metadata.txt", package = "dar")
)
expect_error(read_data(duplicated_files))
