# Return legacy taxonomic information from a recipe

`get_tax()` is deprecated. Use `get_model(rec)$tax_level` for modeled
recipes.

## Usage

``` r
get_tax(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

Tibble containing `tax_info`.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
suppressWarnings(get_tax(rec))
#> # A tibble: 1 × 1
#>   tax_lev
#>   <chr>  
#> 1 Species
```
