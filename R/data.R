#' Phyloseq object from metaHIV project
#'
#' A Phyloseq object containing abundance counts and sample_data for metaHIV
#' project. Count reads were annotated with Metaphlan3.
#'
#' @autoglobal
#' @usage data("metaHIV_phy")
#' @source s3://fcatala-09142020-eu-west-1/cloud_test/SpeciesQuantification/Kraken2
#' @format A phyloseq object with 451 taxas, 30 samples, 3 sample variables and
#'   7 taxonomic ranks.
"metaHIV_phy"

#' Pastery data for step id generation
#'
#' Tibble contain
#' 
#' @name pastry_df
#' @keywords internal
#' @autoglobal
#' @format A tbl_df object with 228 unique pasteries.
#' @source https://raw.githubusercontent.com/prasertcbs/basic-dataset/master/pastry.csv
NULL

#' Recipe for metaHIV_phy data
#'
#' A recipe created for a metaHIV_phy object uning "Riskgroup2" as a `var_info`
#' and "Genus" as a `tax_info`.
#'
#' @format A recipe object.
#' @autoglobal
#' @usage data("test_rec")
"test_rec"

#' Prep-recipe for metaHIV_phy data
#'
#' A recipe created for a metaHIV_phy object uning "Riskgroup2" as a `var_info`
#' and "Genus" as a `tax_info`. Also includes `step_deseq`, `step_maaslin` and
#' `step_metagenomeSeq`.
#'
#' @format A prep-recipe object.
#' @autoglobal
#' @usage data("test_prep_rec")
"test_prep_rec"
