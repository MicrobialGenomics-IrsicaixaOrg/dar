# Add a legacy taxonomic level of interest to a recipe

`add_tax()` is deprecated. Supply `tax_level` to
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
instead.

## Usage

``` r
add_tax(rec, tax_info)
```

## Arguments

- rec:

  A `Recipe` object (must not be prepped).

- tax_info:

  A character string of taxonomic levels that will be used in any
  context.

## Value

A `Recipe` object.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
rec <- suppressWarnings(add_tax(rec, tax_info = "Species"))
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ℹ taxonomic level Species 
#> 
#> 
```
