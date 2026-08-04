# dar

## Introduction

Differential abundance testing in microbiome data challenges both
parametric and non-parametric statistical methods, due to its sparsity,
high variability and compositional nature. Microbiome-specific
statistical methods often assume classical distribution models or take
into account compositional specifics. These produce results across the
specificity-versus-sensitivity space, while type I and type II errors
are difficult to ascertain in real microbiome data when a single method
is used. Consensus across multiple differential abundance (DA) methods
can increase robustness.

With dar, you can use dplyr-like pipeable sequences of DA methods and
then apply different consensus strategies. In this way we can obtain
more reliable results in a fast, consistent and reproducible way.

## Installation

You can install the development version of dar from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pkg_install("MicrobialGenomics-IrsicaixaOrg/dar")
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
  recipe(metaHIV_phy) |>
  add_model(
    ~ RiskGroup2,
    targets = "RiskGroup2",
    tax_level = "Species"
  ) |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.03 * length(x))") |>
  step_maaslin() |>
  step_aldex()

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#>
#>      ℹ ~RiskGroup2
#>
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Nunt 
#>      ◉ step_filter_taxa() id = filter_taxa__Schaumrolle 
#> 
#> DA steps:
#> 
#>      ◉ step_maaslin() id = maaslin__Crêpes_Suzette 
#>      ◉ step_aldex() id = aldex__Pineapple_bun

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
#>      ✔ maaslin__Crêpes_Suzette diff_taxa = 52 
#>      ✔ aldex__Pineapple_bun diff_taxa = 98 
#> 
#>      ℹ 35 taxa are present in all tested methods

## Uniform, reproducible results across engines
tidy_results(da_results) |>
  dplyr::select(
    taxa_id, comparison, effect_size, effect_metric,
    adj_p_value, significant, method, step_id
  )

## Advanced native options remain explicit and reproducible
rec <- step_deseq(
  rec,
  engine_args = list(fit = list(minReplicatesForReplace = Inf))
)

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
#>      ✔ maaslin__Crêpes_Suzette diff_taxa = 52 
#>      ✔ aldex__Pineapple_bun diff_taxa = 98 
#> 
#>      ℹ 35 taxa are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: 2, weights: NULL, exclude: NULL, id: bake__Ladies_navels

## Results
cool(da_results)
#> ℹ Bake for count_cutoff = 2
#> # A tibble: 35 × 2
#>    taxa_id taxa                        
#>    <chr>   <chr>                       
#>  1 Otu_78  Bacteroides_uniformis       
#>  2 Otu_88  Odoribacter_splanchnicus    
#>  3 Otu_119 Alistipes_putredinis        
#>  4 Otu_129 Parabacteroides_merdae      
#>  5 Otu_125 Parabacteroides_distasonis  
#>  6 Otu_82  Barnesiella_intestinihominis
#>  7 Otu_96  Prevotella_copri            
#>  8 Otu_51  Bacteroides_dorei           
#>  9 Otu_332 Catenibacterium_mitsuokai   
#> 10 Otu_62  Bacteroides_ovatus          
#> # ℹ 25 more rows
```

## Contributing

- If you think you have encountered a bug, please [submit an
  issue](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues).

- Either way, learn how to create and share a
  [reprex](https://reprex.tidyverse.org/articles/articles/learn-reprex.html)
  (a minimal, reproducible example), to clearly communicate about your
  code.

- Working on your first Pull Request? You can learn how from this *free*
  series [How to Contribute to an Open Source Project on
  GitHub](https://kcd.im/pull-request)

## Code of Conduct

Please note that the dar project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
