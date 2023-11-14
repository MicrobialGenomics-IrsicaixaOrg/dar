#' Phyloseq object from metaHIV project
#'
#' A Phyloseq object containing abundance counts and sample_data for metaHIV
#' project. Count reads were annotated with Metaphlan3.
#'
#' @format A phyloseq object with 451 taxas, 30 samples, 3 sample variables and
#'   7 taxonomic ranks.
#'
#' @source
#'   s3://fcatala-09142020-eu-west-1/cloud_test/SpeciesQuantification/Kraken2
#' @name data-metaHIV_phy
#' @aliases metaHIV_phy
#' @docType data
#' @keywords data
#' @return The function does not return a value explicitly.
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' metaHIV_phy
"metaHIV_phy"

#' Pastery data for step id generation
#'
#' Tibble contain
#'
#' @format A tbl_df object with 228 unique pasteries.
#'
#' @source
#' https://raw.githubusercontent.com/prasertcbs/basic-dataset/master/pastry.csv
#'
#' @name data-pastry_df
#' @aliases pastry_df
#' @docType data
#' @keywords data internal
#' @autoglobal
#' @return The function does not return a value explicitly.
NULL

#' Recipe for metaHIV_phy data
#'
#' A recipe created for a metaHIV_phy object uning "Riskgroup2" as a `var_info`
#' and "Genus" as a `tax_info`.
#'
#' @format A recipe object.
#'
#' @name data-test_rec
#' @aliases test_rec
#' @docType data
#' @keywords data
#' @autoglobal
#' @return The function does not return a value explicitly.
"test_rec"

#' Prep-recipe for metaHIV_phy data
#'
#' A recipe created for a metaHIV_phy object uning "Riskgroup2" as a `var_info`
#' and "Genus" as a `tax_info`. Also includes `step_deseq`, `step_maaslin` and
#' `step_metagenomeSeq`.
#'
#' @format A prep-recipe object.
#'
#' @name data-test_prep_rec
#' @aliases test_prep_rec
#' @docType data
#' @keywords data
#' @autoglobal
#' @return The function does not return a value explicitly.
"test_prep_rec"
