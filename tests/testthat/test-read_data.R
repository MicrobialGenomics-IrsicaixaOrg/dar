test_that("read phyloseq from rds file", {
  good_file <- system.file("extdata", "metaHIV_phy.rds", package = "dar")
  bad_file <- system.file("extdata", "invented_file.rds", package = "dar")

  expect_s4_class(read_data(good_file), "phyloseq")
  expect_error(read_data(bad_file))
})

test_that("read phyloseq from multiple files", {
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
  expect_error(read_data(only_otu_phy))

  duplicated_files <- c(
    system.file("extdata", "metaHIV_counts.txt", package = "dar"),
    system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
    system.file("extdata", "metaHIV_metadata.txt", package = "dar")
  )
  expect_error(read_data(duplicated_files))
})

test_that("validate phyloseq content", {
  phy <- read_data(system.file("extdata", "metaHIV_phy.rds", package = "dar"))

  expect_equal(dim(phy@otu_table), c(451, 156))
  expect_equal(dim(phy@sam_data), c(156, 3))
  expect_equal(dim(phy@tax_table), c(451, 7))
})

test_that("phyloseq validation", {
  phy <-
    system.file("extdata", "metaHIV_phy.rds", package = "dar") %>%
    read_phyloseq()

  no_tax_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::sample_data(phy))
  no_sam_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::tax_table(phy))
  only_otu_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy))

  expect_s4_class(validate_phyloseq(phy), "phyloseq")
  expect_s4_class(validate_phyloseq(no_sam_phy, "tax_table"), "phyloseq")
  expect_s4_class(validate_phyloseq(no_tax_phy, "sample_data"), "phyloseq")
  expect_s4_class(validate_phyloseq(only_otu_phy, ""), "otu_table")
  expect_error(validate_phyloseq(no_tax_phy))
  expect_error(validate_phyloseq(no_sam_phy))
  expect_error(validate_phyloseq(only_otu_phy))
  expect_error(validate_phyloseq(no_tax_phy, "tax_table"))
  expect_error(validate_phyloseq(no_sam_phy, "sample_data"))
})

test_that("validate otu counts", {
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
})

test_that("validate sample data", {
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
})

test_that("validate taxa table", {
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
})
