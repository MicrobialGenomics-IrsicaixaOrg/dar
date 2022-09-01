
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dar <a href="https://microbialgenomics-irsicaixaorg.github.io/dar/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/actions)
[![Codecov test
coverage](https://codecov.io/gh/MicrobialGenomics-IrsicaixaOrg/dar/branch/main/graph/badge.svg)](https://app.codecov.io/gh/MicrobialGenomics-IrsicaixaOrg/dar?branch=main)
[![PRs
Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](https://makeapullrequest.com)
[![GitHub
issues](https://img.shields.io/github/issues/MicrobialGenomics-IrsicaixaOrg/dar)](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/MicrobialGenomics-IrsicaixaOrg/dar)](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/pulls)
<!-- badges: end -->

## Introduction

Differential abundance testing in microbiome data challenges both
parametric and non-parametric statistical methods, due to its sparsity,
high variability and compositional nature. Microbiome-specific
statistical methods often assume classical distribution models or take
into account compositional specifics. These produce results that range
within the specificity vs sensitivity space in such a way that type I
and type II error are difficult to ascertain in real microbiome data
when a single method is used. Recently, a consensus approach based on
multiple differential abundance (DA) methods was recently suggested in
order to increase robustness.

With dar, you can use dplyr-like pipeable sequences of DA methods and
then apply different consensus strategies. In this way we can obtain
more reliable results in a fast, consistent and reproducible way.

## Installation

You can install the development version of dar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MicrobialGenomics-IrsicaixaOrg/dar")
```

## Usage

``` r
library(dar)
#> Registered S3 methods overwritten by 'vegan':
#>   method         from      
#>   reorder.hclust seriation 
#>   rev.hclust     dendextend
data("metaHIV_phy")

## Define recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species") %>%
  step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.03 * length(x))") %>%
  step_metagenomeseq(rm_zeros = 0.01) %>%
  step_maaslin()

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Sad_cake 
#>      ◉ step_filter_taxa() id = filter_taxa__Malsouka 
#> 
#> DA steps:
#> 
#>      ◉ step_metagenomeseq() id = metagenomeseq__Moorkop 
#>      ◉ step_maaslin() id = maaslin__Flaugnarde

## Prep recipe
da_results <- prep(rec, parallel = TRUE)
da_results
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Results:
#> 
#>      ✔ metagenomeseq__Moorkop diff_taxa = 278 
#>      ✔ maaslin__Flaugnarde diff_taxa = 205 
#> 
#>      ℹ 205 taxa are present in all tested methods

## Consensus strategy
n_methods <- 2
da_results <- bake(da_results, count_cutoff = n_methods)
da_results
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Results:
#> 
#>      ✔ metagenomeseq__Moorkop diff_taxa = 278 
#>      ✔ maaslin__Flaugnarde diff_taxa = 205 
#> 
#>      ℹ 205 taxa are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: 2, weights: NULL, exclude: NULL, id: bake__Curry_puff

## Results
cool(da_results)
#> ℹ Bake for count_cutoff = 2
#> # A tibble: 205 × 2
#>    taxa_id taxa                             
#>    <chr>   <chr>                            
#>  1 Otu_1   Methanobrevibacter_smithii       
#>  2 Otu_2   Methanosphaera_stadtmanae        
#>  3 Otu_10  Bifidobacterium_adolescentis     
#>  4 Otu_12  Bifidobacterium_animalis         
#>  5 Otu_13  Bifidobacterium_bifidum          
#>  6 Otu_15  Bifidobacterium_catenulatum      
#>  7 Otu_18  Bifidobacterium_longum           
#>  8 Otu_19  Bifidobacterium_pseudocatenulatum
#>  9 Otu_34  Olsenella_scatoligenes           
#> 10 Otu_35  Collinsella_aerofaciens          
#> # … with 195 more rows
```

## Contributing

-   If you think you have encountered a bug, please [submit an
    issue](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues).

-   Either way, learn how to create and share a
    [reprex](https://reprex.tidyverse.org/articles/articles/learn-reprex.html)
    (a minimal, reproducible example), to clearly communicate about your
    code.

-   Working on your first Pull Request? You can learn how from this
    *free* series [How to Contribute to an Open Source Project on
    GitHub](https://kcd.im/pull-request)

## Code of Conduct

Please note that the dar project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
