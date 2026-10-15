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
