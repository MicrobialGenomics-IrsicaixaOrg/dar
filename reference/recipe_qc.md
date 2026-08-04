# Recipe quality-control metrics

`recipe_qc()` summarizes sparsity and sequencing-depth metrics from the
count matrix stored in a
[Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
or
[PrepRecipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
object. The summary is independent of the centralized statistical model
and can therefore be used before model configuration and after recipe
preprocessing.

## Usage

``` r
recipe_qc(rec, group_by = NULL)
```

## Arguments

- rec:

  A
  [Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
  or
  [PrepRecipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
  object.

- group_by:

  `NULL` or the name of one sample metadata column used to stratify the
  QC summary.

## Value

A tibble with grouping identifiers and QC metrics.

## Details

With `group_by = NULL`, the result contains one overall row. When a
sample metadata column is supplied, the overall row is followed by one
row for each observed group. The `group_by` and `group` columns are `NA`
in the overall row, so a real group named `"all"` remains unambiguous.

The returned metrics are:

- `n`: number of values in the count matrix.

- `n_zero`: number of zero values.

- `pct_zero`: percentage of zero values.

- `pct_all_zero`: percentage of taxa absent from every sample.

- `pct_singletons`: percentage of taxa detected in exactly one sample.

- `pct_doubletons`: percentage of taxa detected in exactly two samples.

- `n_samples`: number of samples.

- `lib_size_min`: minimum sample library size.

- `lib_size_max`: maximum sample library size.

- `count_mean`: mean of the per-sample mean counts.

- `count_max`: mean of the per-sample maximum counts.

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy)

## Overall metrics do not require a statistical model
recipe_qc(rec)
#> # A tibble: 1 × 13
#>   group_by group     n n_zero pct_zero pct_all_zero pct_singletons
#>   <chr>    <chr> <dbl>  <dbl>    <dbl>        <dbl>          <dbl>
#> 1 NA       NA    70356  57632     81.9            0           20.6
#> # ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,
#> #   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>

## Add explicit sample-level stratification when useful
recipe_qc(rec, group_by = "RiskGroup2")
#> # A tibble: 4 × 13
#>   group_by   group     n n_zero pct_zero pct_all_zero pct_singletons
#>   <chr>      <chr> <dbl>  <dbl>    <dbl>        <dbl>          <dbl>
#> 1 NA         NA    70356  57632     81.9          0             20.6
#> 2 RiskGroup2 hts   18491  15108     81.7         24.2           22.8
#> 3 RiskGroup2 msm   45100  37019     82.1         16.0           20.2
#> 4 RiskGroup2 pwid   6765   5505     81.4         41.2           16.6
#> # ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,
#> #   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>
```
