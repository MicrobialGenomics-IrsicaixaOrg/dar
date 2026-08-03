# Get the centralized statistical model from a recipe

Get the centralized statistical model from a recipe

## Usage

``` r
get_model(rec)
```

## Arguments

- rec:

  A
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  or
  [`PrepRecipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep_recipe.md)
  object.

## Value

`NULL` when no model is defined, otherwise a list containing the
formula, targets, taxonomic level, and its contrast, reference,
missing-value and time controls.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
get_model(rec)
#> $formula
#> ~RiskGroup2
#> <environment: 0x5556cfec3160>
#> 
#> $targets
#> [1] "RiskGroup2"
#> 
#> $tax_level
#> [1] "Species"
#> 
#> $time
#> NULL
#> 
#> $reference
#> NULL
#> 
#> $contrasts
#> [1] "simple"
#> 
#> $na_action
#> [1] "error"
#> 
#> $max_time_levels
#> [1] 12
#> 
```
