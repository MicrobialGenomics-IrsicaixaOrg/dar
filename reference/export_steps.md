# Export a recipe analysis configuration as versioned JSON

Writes the centralized model, ordered preprocessing and DA steps, and
any bake configurations using the non-executable `dar-recipe` JSON
schema.

## Usage

``` r
export_steps(rec, file_name)
```

## Arguments

- rec:

  A
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  or
  [`PrepRecipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep_recipe.md)
  object.

- file_name:

  Output path.

## Value

`file_name`, invisibly.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  step_filter_by_prevalence(threshold = 0.1)
path <- tempfile(fileext = ".json")
export_steps(rec, path)
```
