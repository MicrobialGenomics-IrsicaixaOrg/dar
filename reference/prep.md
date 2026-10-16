# Performs all the steps defined in a Recipe

For a Recipe with at least one preprocessing or DA operation run the
steps in a convenient order.

## Usage

``` r
prep(
  rec,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16"),
  force = FALSE
)
```

## Arguments

- rec:

  A `Recipe` object.

- parallel:

  if FALSE, no palatalization. if TRUE, parallel execution using future
  and furrr packages.

- workers:

  Number of workers for palatalization.

- force:

  Force the reexecution of all steps. This remove previous results.

## Value

A `PrepRecipe` object.

## Details

Recipes containing differential-abundance steps require an explicit
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
before execution. Preprocessing-only recipes remain supported without a
model. A missing model is reported before preprocessing or dependency
checks, so the input data are not modified.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Class") |>
  step_filter_by_prevalence(threshold = 0.03) |>
  step_wilcox()

da_results <- prep(rec, parallel = FALSE)
```
