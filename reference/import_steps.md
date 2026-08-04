# Import a recipe configuration from versioned JSON

Recipe files are treated as untrusted input. The complete document is
decoded and validated before supported constructors are applied.
Importing bake configurations into an unprepared recipe never runs
analysis unless `prepare = TRUE` is supplied explicitly.

## Usage

``` r
import_steps(
  rec,
  file,
  prepare = FALSE,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16")
)
```

## Arguments

- rec:

  A
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  or
  [`PrepRecipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep_recipe.md)
  object.

- file:

  Input JSON path.

- prepare:

  Whether an unprepared recipe may run
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  when the file contains bake configurations. Defaults to `FALSE`.

- parallel:

  Passed to
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  only when `prepare = TRUE`.

- workers:

  Passed to
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  only when `prepare = TRUE`.

## Value

The configured `Recipe` or `PrepRecipe`.

## Examples

``` r
data(metaHIV_phy)
configured <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  step_filter_by_prevalence(threshold = 0.1)
path <- tempfile(fileext = ".json")
export_steps(configured, path)
restored <- import_steps(recipe(metaHIV_phy), path)
```
