# Phyloseq Quality Control Metrics

phy_qc() returns a tibble. It will have information about some important
metrics about the sparsity of the count matrix. The content of the table
is as follows:

- var_levels: levels of the categorical variable of interest. "all"
  refers to all rows of the dataset (without splitting by categorical
  levels).

- n: total number of values in the count matrix.

- n_zero: number of zeros in the count matrix.

- pct_zero: percentage of zeros in the count matrix.

- pct_all_zero: percentage of taxa with zero counts in all samples.

- pct_singletons: percentage of taxa with counts in a single sample.

- pct_doubletons: percentage of taxa with counts in two samples.

- count_mean: average of the mean counts per sample.

- count_min: average of the min counts per sample.

- count_max: average of the max counts per sample.

## Usage

``` r
phy_qc(rec)

# S4 method for class 'Recipe'
phy_qc(rec)
```

## Arguments

- rec:

  A Recipe or Recipe step.

## Value

A tibble

## Examples

``` r
data(metaHIV_phy)

## Define Recipe
rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")

phy_qc(rec)
#> # A tibble: 4 × 10
#>   var_levels     n n_zero pct_zero pct_all_zero pct_singletons pct_doubletons
#>   <chr>      <int>  <int>    <dbl>        <dbl>          <dbl>          <dbl>
#> 1 all        70356  57632     81.9          0             20.6           8.87
#> 2 hts        18491  15108     81.7         24.2           22.8           8.43
#> 3 msm        45100  37019     82.1         16.0           20.2           9.53
#> 4 pwid        6765   5505     81.4         41.2           16.6           9.31
#> # ℹ 3 more variables: count_mean <dbl>, count_min <dbl>, count_max <dbl>
```
