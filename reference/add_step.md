# Add a New Operation to the Current Recipe

`add_step` adds a step to the last location in the Recipe. `add_check`
does the same for checks.

## Usage

``` r
add_step(rec, object)
```

## Arguments

- rec:

  A
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  or
  [`PrepRecipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep_recipe.md).

- object:

  A step or check object.

## Value

A updated
[`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
or
[`PrepRecipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep_recipe.md)
with the new operation in the last slot.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")

# Internally, step_ functions use add_step to append themselves
rec <- step_maaslin(rec)
length(rec@steps) # Returns 1
#> [1] 1

# If we try to add the exact same step, it will inform us and skip it
rec <- step_maaslin(rec)
#> ! This step is already defined with the same parameters and will be skipped.
#> ℹ `rec %>% step_maaslin(min_abundance = 0, min_prevalence = 0.1, min_variance =
#>   0, normalization = "TSS", transform = "LOG", max_significance = 0.1,
#>   random_effects = NULL, correction = "BH", standardize = TRUE, reference =
#>   NULL, median_comparison_abundance = TRUE, rarefy = FALSE)`
length(rec@steps) # Still returns 1
#> [1] 1
```
