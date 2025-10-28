# Resample an OTU table such that all samples have the same library size.

Please note that the authors of phyloseq do not advocate using this as a
normalization procedure, despite its recent popularity. Our
justifications for using alternative approaches to address disparities
in library sizes have been made available as an article in PLoS
Computational Biology. See phyloseq_to_deseq2 for a recommended
alternative to rarefying directly supported in the phyloseq package, as
well as the supplemental materials for the PLoS-CB article and the
phyloseq extensions repository on GitHub. Nevertheless, for comparison
and demonstration, the rarefying procedure is implemented here in good
faith and with options we hope are useful. This function uses the
standard R sample function to resample from the abundance values in the
otu_table component of the first argument, physeq. Often one of the
major goals of this procedure is to achieve parity in total number of
counts between samples, as an alternative to other formal normalization
procedures, which is why a single value for the sample.size is expected.
This kind of resampling can be performed with and without replacement,
with replacement being the more computationally-efficient, default
setting. See the replace parameter documentation for more details. We
recommended that you explicitly select a random number generator seed
before invoking this function, or, alternatively, that you explicitly
provide a single positive integer argument as rngseed.

## Usage

``` r
step_rarefaction(rec, id = rand_id("rarefaction"))

# S4 method for class 'Recipe'
step_rarefaction(rec, id = rand_id("rarefaction"))

# S4 method for class 'PrepRecipe'
step_rarefaction(rec, id = rand_id("rarefaction"))
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- id:

  A character string that is unique to this step to identify it.

## Value

An object of class `Recipe`

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- 
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.03 * length(x))")

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
#>      ◉ step_subset_taxa() id = subset_taxa__Alexandertorte 
#>      ◉ step_filter_taxa() id = filter_taxa__Baklava 
#> 
#> DA steps:
#> 

## Define step with default parameters and prep
rec <- step_rarefaction(rec) 
  
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
#>      ◉ step_subset_taxa() id = subset_taxa__Alexandertorte 
#>      ◉ step_filter_taxa() id = filter_taxa__Baklava 
#>      ◉ step_rarefaction() id = rarefaction__Cronut 
#> 
#> DA steps:
#> 
```
