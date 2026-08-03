# Add a centralized statistical model to a recipe

The model is the single source of truth for statistical design
arguments. `targets` are the variables of inference; other fixed terms
are adjustment variables. `tax_level` defines the taxonomic resolution
used by DA methods. When `time` is supplied, numeric values are treated
as ordered categorical time points while compiling the model.

## Usage

``` r
add_model(
  rec,
  formula,
  targets = NULL,
  tax_level = NULL,
  time = NULL,
  reference = NULL,
  contrasts = c("simple", "all"),
  na_action = c("error", "drop"),
  max_time_levels = 12
)
```

## Arguments

- rec:

  An unprepared
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  object.

- formula:

  A one-sided formula or `NULL` to remove the model.

- targets:

  One or more categorical metadata columns that define the statistical
  targets. Every target must occur in the fixed formula.

- tax_level:

  A single taxonomic rank used by differential-abundance methods.

- time:

  Optional single metadata column identifying time.

- reference:

  Optional named character vector defining factor reference levels.

- contrasts:

  Either `"simple"` for pairwise/simple effects or `"all"` to
  additionally request difference-in-differences contrasts.

- na_action:

  Either `"error"` or `"drop"`. Dropping is performed once for all
  engines during
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md).

- max_time_levels:

  Maximum number of time points that may be created by automatic
  factorization.

## Value

The updated `Recipe`.

## Details

For compatibility with recipes created before version 1.9.3, omitted
`targets` or `tax_level` values are recovered from the deprecated
selector slots when available. New recipes must supply both arguments
explicitly.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(
    ~ RiskGroup2 + Cluster,
    targets = "RiskGroup2",
    tax_level = "Species"
  )
get_model(rec)
#> $formula
#> ~RiskGroup2 + Cluster
#> <environment: 0x5640ee15d418>
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
