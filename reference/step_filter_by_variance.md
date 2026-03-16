# Filter taxa by variance

This is a convenience function around the filter_taxa function. It is
designed to speed up filtering complex experimental objects with one
function call. In the case of run_filter_by_variance, the filtering will
be based on the variance of each taxon. The taxa retained in the dataset
are those where the variance of their abundance is greater than the
provided threshold.

## Usage

``` r
step_filter_by_variance(
  rec,
  threshold = 0.01,
  id = rand_id("filter_by_variance")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- threshold:

  The variance threshold for filtering taxa. The default value is 0.01.

- id:

  A character string that is unique to this step to identify it.

## Value

A Recipe object that has been filtered based on variance.

## Details

The function calculates the variance of all taxa in the phyloseq object.
It then compares this variance to the variance of each individual taxon.
If a taxon's variance is less than or equal to the threshold, that taxon
is removed from the phyloseq object.

## Note

This function modifies `rec` in place, you might want to make a copy of
`rec` before modifying it if you need to preserve the original object.

## See also

[`filter_taxa`](https://rdrr.io/pkg/phyloseq/man/filter_taxa.html)

Other filter phy steps:
[`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md),
[`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- recipe(metaHIV_phy, "RiskGroup2", "Phylum")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> 

## Define step_filter_by_variance step with default parameters
rec <- step_filter_by_variance(rec, threshold = 0.01)
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_filter_by_variance() id = filter_by_variance__Gulab_jamun 
#> 
#> DA steps:
#> 
```
