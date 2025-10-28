# Filter taxa by abundance

This is a convenience wrapper around the filter_taxa function. It is
intended to speed up filtering complex experimental objects with one
function call. In the case of filter_by_abundance, the filtering will be
based on the relative abundance of each taxon. The taxa retained in the
dataset are those where the sum of their abundance is greater than the
product of the total abundance and the provided threshold.

## Usage

``` r
step_filter_by_abundance(
  rec,
  threshold = 0.01,
  id = rand_id("filter_by_abundance")
)

# S4 method for class 'Recipe'
step_filter_by_abundance(
  rec,
  threshold = 0.01,
  id = rand_id("filter_by_abundance")
)

# S4 method for class 'PrepRecipe'
step_filter_by_abundance(
  rec,
  threshold = 0.01,
  id = rand_id("filter_by_abundance")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- threshold:

  The relative abundance threshold for filtering taxa, expressed as a
  proportion of the total abundance. For example, a threshold of 0.01
  means that a taxon must make up at least 1% of the total abundance to
  be retained. The default value is 0.01.

- id:

  A character string that is unique to this step to identify it.

## Value

A Recipe object that has been filtered based on abundance.

## Details

The function calculates the total abundance of all taxa in the phyloseq
object. It then compares this total abundance to the abundance of each
individual taxon. If a taxon's abundance is less than the threshold
times the total abundance, that taxon is removed from the phyloseq
object.

## Note

This function modifies `rec` in place, you might want to make a copy of
`rec` before modifying it if you need to preserve the original object.

## See also

[`filter_taxa`](https://rdrr.io/pkg/phyloseq/man/filter_taxa.html)

Other filter phy steps:
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md),
[`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md),
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

## Define filter_by_abundance step with default parameters
rec <- step_filter_by_abundance(rec, threshold = 0.01)
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
#>      ◉ step_filter_by_abundance() id = filter_by_abundance__Mooncake 
#> 
#> DA steps:
#> 
```
