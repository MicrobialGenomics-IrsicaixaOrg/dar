# Phyloseq Quality Control Metrics

`phy_qc()` returns a tibble containing important metrics about the
sparsity and sequencing depth of the count matrix. These metrics are
calculated both globally and split by the levels of the categorical
variable of interest. The content of the table is as follows:

## Usage

``` r
phy_qc(rec)
```

## Arguments

- rec:

  A `Recipe` object.

## Value

A tibble with QC metrics.

## Details

- `var_levels`: levels of the categorical variable of interest. "all"
  refers to all rows of the dataset (without splitting by categorical
  levels).

- `n`: total number of values in the count matrix.

- `n_zero`: number of zeros in the count matrix.

- `pct_zero`: percentage of zeros in the count matrix.

- `pct_all_zero`: percentage of taxa with zero counts in all samples.

- `pct_singletons`: percentage of taxa with counts in a single sample.

- `pct_doubletons`: percentage of taxa with counts in exactly two
  samples.

- `n_samples`: total number of samples in the group.

- `lib_size_min`: minimum library size (total counts) across samples in
  the group.

- `lib_size_max`: maximum library size (total counts) across samples in
  the group.

- `count_mean`: average of the mean counts per sample.

- `count_max`: average of the max counts per sample.

## Examples

``` r
data(metaHIV_phy)

## 1. Init Recipe
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")

## 2. Get QC metrics
phy_qc(rec)
#> # A tibble: 4 × 12
#>   var_levels     n n_zero pct_zero pct_all_zero pct_singletons pct_doubletons
#>   <chr>      <int>  <int>    <dbl>        <dbl>          <dbl>          <dbl>
#> 1 all        70356  57632     81.9          0             20.6           8.87
#> 2 hts        18491  15108     81.7         24.2           22.8           8.43
#> 3 msm        45100  37019     82.1         16.0           20.2           9.53
#> 4 pwid        6765   5505     81.4         41.2           16.6           9.31
#> # ℹ 5 more variables: n_samples <int>, lib_size_min <dbl>, lib_size_max <dbl>,
#> #   count_mean <dbl>, count_max <dbl>
```
