#' Phyloseq object from metaHIV project
#'
#' A phyloseq object containing the abundance table, sample metadata and
#' taxonomic annotations supplied for the metaHIV project associated with
#' Noguera-Julian et al. (2016). It is provided as a software example. The
#' profiling pipeline and the measurement units of this derived table have
#' not been independently verified; integer values alone do not establish
#' that a table contains raw read counts.
#'
#' @autoglobal
#' @usage data("metaHIV_phy")
#' @source Noguera-Julian et al. (2016),
#'   \doi{10.1016/j.ebiom.2016.01.032}.
#' @format A phyloseq object with 451 taxa, 156 samples, 3 sample variables and
#'   7 taxonomic ranks.
"metaHIV_phy"

#' Pastery data for step id generation
#'
#' Tibble contain
#' 
#' @name pastry_df
#' @keywords internal
#' @autoglobal
#' @return tibble with pastry names
#' @format A tbl_df object with 228 unique pasteries.
#' @source https://raw.githubusercontent.com/prasertcbs/basic-dataset/master/pastry.csv
NULL

#' Recipe for metaHIV_phy data
#'
#' An unprepared recipe for `metaHIV_phy`, with a centralized `~ RiskGroup2`
#' model at genus level and `hts` as its reference. It includes microbial
#' subsetting, feature-prevalence filtering, MaAsLin3 and DESeq2 steps.
#'
#' @format A Recipe object.
#' @autoglobal
#' @usage data("test_rec")
"test_rec"

#' PrepRecipe for metaHIV_phy data
#'
#' A prepared version of [test_rec], fitted through the centralized model API
#' at genus level. Results retain planned pairwise `RiskGroup2` contrasts and
#' method-specific effects. Original stored features remain available through
#' [otu_table()]; [tidy_abundance()] returns their shared genus-level view.
#' Unavailable native estimates remain missing; [effect_plt()] omits them
#' rather than treating them as zero effects.
#' Generated reproducibly by `data-raw/metaHIV_recipes.R`.
#'
#' @format A PrepRecipe object.
#' @autoglobal
#' @usage data("test_prep_rec")
"test_prep_rec"
