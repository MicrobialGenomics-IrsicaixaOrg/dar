# Subset taxa by taxonomic level

This is a convenience function around the subset_taxa function from the
phyloseq package. It is designed to speed up subsetting complex
experimental objects with one function call. In the case of
run_subset_taxa, the subsetting will be based on the taxonomic level of
each taxon. The taxa retained in the dataset are those where the
taxonomic level matches the provided taxa.

## Usage

``` r
step_subset_taxa(rec, tax_level, taxa, id = rand_id("subset_taxa"))

# S4 method for class 'Recipe'
step_subset_taxa(rec, tax_level, taxa, id = rand_id("subset_taxa"))

# S4 method for class 'PrepRecipe'
step_subset_taxa(rec, tax_level, taxa, id = rand_id("subset_taxa"))
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- tax_level:

  The taxonomic level for subsetting taxa.

- taxa:

  The taxa to be retained in the dataset.

- id:

  A character string that is unique to this step to identify it.

## Value

A Recipe object that has been subsetted based on taxonomic level.

## Details

The function subsets the taxa in the phyloseq object based on the
provided taxonomic level and taxa. Only the taxa that match the provided
taxa at the given taxonomic level are retained in the phyloseq object.

## Note

This function modifies `rec` in place, you might want to make a copy of
`rec` before modifying it if you need to preserve the original object.

## See also

[`subset_taxa`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> 

## Define step_subset_taxa step with default parameters
rec <- step_subset_taxa(
  rec, 
  tax_level = "Kingdom",
  taxa = c("Bacteria", "Archaea")
)
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
#>      ◉ step_subset_taxa() id = subset_taxa__Curry_puff 
#> 
#> DA steps:
#> 
```
