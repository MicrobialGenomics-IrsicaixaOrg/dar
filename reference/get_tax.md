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
