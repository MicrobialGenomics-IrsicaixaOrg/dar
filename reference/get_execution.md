# Inspect the recorded execution of a prepared recipe

`get_execution()` exposes the manifest recorded by
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
without recomputing contrasts, modifying results or rerunning any
analysis. Samples listed in `dropped_samples` were removed by the
centralized model's missing-value policy, not by method-specific
filtering.

## Usage

``` r
get_execution(rec)
```

## Arguments

- rec:

  A
  [PrepRecipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
  object.

## Value

A list containing the recorded `contrasts`, `dropped_samples`,
`executed_steps` and `skipped_steps`. The contrasts and skipped-step
tables retain their existing columns, including incompatibility reasons.
No centralized model also yields `NULL` contrasts in current recipes.
When supplied by an engine, an additional `engine_details` list records
its effective settings and package version. ALDEx3 records its seed,
scale-model basis, Monte Carlo settings, fitting engine and elapsed
time.

## Details

Older prepared recipes may lack all or part of this manifest.
Unavailable fields are returned as `NULL`: missing provenance does not
imply that no samples were removed or no steps were skipped.

## Examples

``` r
data(metaHIV_phy)
prepared <- recipe(metaHIV_phy) |>
  step_filter_by_prevalence(threshold = 0.1) |>
  prep(parallel = FALSE)
get_execution(prepared)
#> $contrasts
#> NULL
#> 
#> $dropped_samples
#> character(0)
#> 
#> $executed_steps
#> character(0)
#> 
#> $skipped_steps
#> # A tibble: 0 × 4
#> # ℹ 4 variables: step_id <chr>, engine <chr>, compatible <lgl>, reason <chr>
#> 
```
