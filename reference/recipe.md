# Create a Recipe for preprocessing data

A Recipe is a description of the steps to be applied to a data set in
order to prepare it for data analysis.

## Usage

``` r
recipe(
  microbiome_object = NULL,
  var_info = NULL,
  tax_info = NULL,
  steps = list()
)
```

## Arguments

- microbiome_object:

  Phyloseq-class object or TreeSummarizedExperiment-class object.

- var_info:

  A character string of column names corresponding to variables that
  will be used in any context. Deprecated; supply `targets` to
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  instead.

- tax_info:

  A character string of taxonomic levels that will be used in any
  context. Deprecated; supply `tax_level` to
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  instead.

- steps:

  list with steps.

## Value

An object of class `Recipe` containing the microbiome object, an
optional centralized model and the configured processing or DA steps.
Legacy selector slots are retained for compatibility during the
deprecation cycle.

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
