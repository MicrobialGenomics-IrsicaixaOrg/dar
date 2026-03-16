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

## Examples

``` r
data(metaHIV_phy)
rec <- 
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Class") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = function(x) sum(x > 0) >= (0.03 * length(x))) |>
  step_maaslin()

da_results <- prep(rec)
```
