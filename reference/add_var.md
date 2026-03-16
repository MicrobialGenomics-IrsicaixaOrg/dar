# Adds variable of interest to the Recipe

Adds variable of interest to the Recipe

## Usage

``` r
add_var(rec, var_info)
```

## Arguments

- rec:

  A `Recipe` object (must not be prepped).

- var_info:

  A character string of column names corresponding to variables that
  will be used in any context.

## Value

A `Recipe` object.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
rec <- add_var(rec, var_info = "RiskGroup2")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ✖ undefined taxonomic level. Use add_tax() to add it to Recipe! 
#> 
```
