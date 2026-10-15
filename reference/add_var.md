# Add legacy variables of interest to a recipe

`add_var()` is deprecated. Supply `targets` to
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
instead.

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
