# Return legacy target information from a recipe

`get_var()` is deprecated. Use `get_model(rec)$targets` for modeled
recipes.

## Usage

``` r
get_var(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

Tibble containing `var_info`.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
suppressWarnings(get_var(rec))
#> # A tibble: 1 × 1
#>   vars      
#>   <chr>     
#> 1 RiskGroup2
```
