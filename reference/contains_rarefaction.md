# Checks if Recipe contains a rarefaction step

Checks if Recipe contains a rarefaction step

## Usage

``` r
contains_rarefaction(rec)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this recipe.

## Value

boolean

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
contains_rarefaction(rec)
#> [1] FALSE
contains_rarefaction(step_rarefaction(rec))
#> [1] TRUE
```
