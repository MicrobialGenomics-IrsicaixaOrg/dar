# Filter taxa based on across-sample OTU abundance criteria

This function is directly analogous to the genefilter function for
microarray filtering, but is used for filtering OTUs from phyloseq
objects. It applies an arbitrary set of functions — as a function list,
for instance, created by filterfun — as across-sample criteria, one OTU
at a time. It takes as input a phyloseq object, and returns a logical
vector indicating whether or not each OTU passed the criteria.
Alternatively, if the "prune" option is set to FALSE, it returns the
already-trimmed version of the phyloseq object.

## Usage

``` r
step_filter_taxa(rec, .f, id = rand_id("filter_taxa"))
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- .f:

  A function or list of functions that take a vector of abundance values
  and return a logical. Some canned useful function types are included
  in the genefilter-package.

- id:

  A character string that is unique to this step to identify it.

## Value

An object of class `Recipe`

## See also

Other filter phy steps:
[`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md),
[`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md)

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

## Define filter taxa step with default parameters
rec <- 
  step_filter_taxa(rec, .f = "function(x) sum(x > 0) >= (0.03 * length(x))")
  
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ✖ undefined taxonomic level. Use add_model() to define the analysis! 
#> Preporcessing steps:
#> 
#>      ◉ step_filter_taxa() id = filter_taxa__Nazook 
#> 
#> DA steps:
#> 
```
