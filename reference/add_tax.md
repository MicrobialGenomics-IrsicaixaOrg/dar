# Adds taxonomic level of interest in the Recipe.

Adds taxonomic level of interest in the Recipe.

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
rec <- add_tax(rec, tax_info = "Species")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined variable of interest. Use add_var() to add it to Recipe! 
#>      ℹ taxonomic level Species 
#> 
#> 
```
