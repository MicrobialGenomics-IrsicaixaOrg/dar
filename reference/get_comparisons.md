# Generate all unique contrasts between levels of a categorical variable.

Generate all unique contrasts between levels of a categorical variable.

## Usage

``` r
get_comparisons(var, phy, as_list = TRUE, n_cut = 1)
```

## Arguments

- var:

  categorical variable

- phy:

  phyloseq object

- as_list:

  boolean indicating if output must be returned as a list.

- n_cut:

  minimum of observations by level.

## Value

tibble or list

## Examples

``` r
data(metaHIV_phy)
dar:::get_comparisons("RiskGroup2", metaHIV_phy)
#> [[1]]
#> [1] "hts" "msm"
#> 
#> [[2]]
#> [1] "hts"  "pwid"
#> 
#> [[3]]
#> [1] "msm"  "pwid"
#> 
```
