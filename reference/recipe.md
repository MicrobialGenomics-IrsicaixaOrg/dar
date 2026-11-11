# Create a Recipe for preprocessing data

A Recipe is a description of the steps to be applied to a data set in
order to prepare it for data analysis.

## Usage

``` r
recipe(microbiome_object = NULL, steps = list())
```

## Arguments

- microbiome_object:

  Phyloseq-class object or TreeSummarizedExperiment-class object.

- steps:

  list with steps.

## Value

An object of class `Recipe` containing the microbiome object, an
optional centralized model and the configured processing or DA steps.

## Examples

``` r
data(metaHIV_phy)

## Configure, prepare and obtain a consensus
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
  step_filter_by_prevalence(threshold = 0.3) |>
  step_deseq() |>
  step_maaslin()

da_results <- prep(rec, parallel = FALSE) |>
  bake()
#> Warning: Estimated rdf < 1.0; not estimating variance
#> Warning: Estimated rdf < 1.0; not estimating variance
#> Warning: Estimated rdf < 1.0; not estimating variance
cool(da_results)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 1 × 9
#>   taxa_id taxa    contrast_id comparison contrast_type var   effect method_count
#>   <chr>   <chr>   <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#> 1 Otu_96  Bacter… RiskGroup2… RiskGroup… main          Risk… down              2
#> # ℹ 1 more variable: methods <chr>

## A recipe without a model can be used for preprocessing
rec <- recipe(metaHIV_phy) |>
  step_filter_by_prevalence(threshold = 0.1)
```
