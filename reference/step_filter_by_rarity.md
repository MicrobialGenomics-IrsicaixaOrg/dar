# Filter taxa by rarity

This is a convenience function around the filter_taxa function. It is
designed to speed up filtering complex experimental objects with one
function call. In the case of run_filter_by_rarity, the filtering will
be based on the rarity of each taxon. The taxa retained in the dataset
are those where the sum of their rarity is less than the provided
threshold.

## Usage

``` r
step_filter_by_rarity(rec, threshold = 0.01, id = rand_id("filter_by_rarity"))

# S4 method for class 'Recipe'
step_filter_by_rarity(rec, threshold = 0.01, id = rand_id("filter_by_rarity"))

# S4 method for class 'PrepRecipe'
step_filter_by_rarity(rec, threshold = 0.01, id = rand_id("filter_by_rarity"))
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- threshold:

  The rarity threshold for filtering taxa, expressed as a proportion of
  the total number of samples. For example, a threshold of 0.01 means
  that a taxon must be present in less than 1% of the samples to be
  retained. The default value is 0.01.

- id:

  A character string that is unique to this step to identify it.

## Value

A Recipe object that has been filtered based on rarity.

## Details

The function calculates the rarity of all taxa in the phyloseq object as
the proportion of samples in which they are present. It then compares
this rarity to the threshold. If a taxon's rarity is greater than the
threshold, that taxon is removed from the phyloseq object.

## Note

This function modifies `rec` in place, you might want to make a copy of
`rec` before modifying it if you need to preserve the original object.

## See also

[`filter_taxa`](https://rdrr.io/pkg/phyloseq/man/filter_taxa.html)

Other filter phy steps:
[`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
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

## Define step_filter_by_rarity step with default parameters
rec <- step_filter_by_rarity(rec, threshold = 0.01)
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
#>      ◉ step_filter_by_rarity() id = filter_by_rarity__Kolache 
#> 
#> DA steps:
#> 
```
