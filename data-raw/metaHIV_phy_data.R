## code to prepare `metaHIV_phy` dataset goes here

# 1/ Phyloseq object as internal data -----------------------------------------------
set.seed(123)
metaHIV_phy <-
    aws.s3::s3readRDS(
        object = "TestMetaHIV/WMGS/SpeciesQuantification/Metaphlan3/Metaphlan3Counts_psObject.rds",
        bucket = "s3://datatoshare-mnoguera/"
    )

phyloseq::sample_data(metaHIV_phy) <-
    phyloseq::sample_data(metaHIV_phy) %>%
    .[, c("RiskGroup2", "Cluster", "HIV_Status")]

usethis::use_data(metaHIV_phy, overwrite = TRUE)

# 2/ Phyloseq and subfiles as raw rds -----------------------------------------------
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)

saveRDS(metaHIV_phy, "inst/extdata/metaHIV_phy.rds")
phyloseq::sample_data(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    saveRDS("inst/extdata/metaHIV_metadata.rds")

phyloseq::get_taxa(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "otu_id") %>%
    saveRDS("inst/extdata/metaHIV_counts.rds")

phyloseq::tax_table(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "otu_id") %>%
    saveRDS("inst/extdata/metaHIV_taxas.rds")

# 3/ Write phyloseq object inputs as txt files --------------------------------------
phyloseq::sample_data(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    data.table::fwrite("inst/extdata/metaHIV_metadata.txt")

phyloseq::get_taxa(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "otu_id") %>%
    data.table::fwrite("inst/extdata/metaHIV_counts.txt")

phyloseq::tax_table(metaHIV_phy) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "otu_id") %>%
    data.table::fwrite("inst/extdata/metaHIV_taxas.txt")
