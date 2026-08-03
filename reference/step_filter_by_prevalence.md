# Filter taxa by prevalence

This is a convenience function around the filter_taxa function. It is
designed to speed up filtering complex experimental objects with one
function call. In the case of run_filter_by_prevalence, the filtering
will be based on the prevalence of each taxon. The taxa retained in the
dataset are those where the prevalence is greater than the provided
threshold.

## Usage

``` r
step_filter_by_prevalence(
  rec,
  threshold = 0.01,
  id = rand_id("filter_by_prevalence")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- threshold:

  The prevalence threshold for filtering taxa, expressed as a proportion
  of the total number of samples. For example, a threshold of 0.01 means
  that a taxon must be present in at least 1% of the samples to be
  retained. The default value is 0.01.

- id:

  A character string that is unique to this step to identify it.

## Value

A Recipe object that has been filtered based on prevalence.

## Details

The function calculates the prevalence of all taxa in the phyloseq
object as the proportion of samples in which they are present. It then
compares this prevalence to the threshold. If a taxon's prevalence is
less than the threshold, that taxon is removed from the phyloseq object.

## Note

This function modifies `rec` in place, you might want to make a copy of
`rec` before modifying it if you need to preserve the original object.

## See also

[`filter_taxa`](https://rdrr.io/pkg/phyloseq/man/filter_taxa.html)

Other filter phy steps:
[`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md),
[`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md),
[`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- recipe(metaHIV_phy)
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ✖ undefined taxonomic level. Use add_model() to define the analysis! 
#> 

## Define step_filter_by_prevalence step with default parameters
rec <- step_filter_by_prevalence(rec, threshold = 0.01)
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ✖ undefined taxonomic level. Use add_model() to define the analysis! 
#> Preporcessing steps:
#> 
#>      ◉ step_filter_by_prevalence() id = filter_by_prevalence__Lotus_seed_bun 
#> 
#> DA steps:
#> 
```
