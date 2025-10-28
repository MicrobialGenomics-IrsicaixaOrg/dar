# Extract outs with all 0 values in at least on level of the variable

Extract outs with all 0 values in at least on level of the variable

## Usage

``` r
zero_otu(obj, var = NULL, pct_cutoff = 0)

# S4 method for class 'Recipe'
zero_otu(obj, var = NULL, pct_cutoff = 0)

# S4 method for class 'phyloseq'
zero_otu(obj, var = NULL, pct_cutoff = 0)
```

## Arguments

- obj:

  A `Recipe` or `phyloseq` object.

- var:

  Variable of interest. Must be present in the metadata.

- pct_cutoff:

  Minimum of pct counts samples with counts for each taxa.

## Value

character vector

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")

## Extract outs with all 0 values
zero_otu(rec)
#> # A tibble: 1,353 × 5
#>    taxa_id RiskGroup2 no_zero total   pct
#>    <chr>   <chr>        <dbl> <int> <dbl>
#>  1 Otu_100 hts              0    41     0
#>  2 Otu_100 pwid             0    15     0
#>  3 Otu_103 msm              0   100     0
#>  4 Otu_103 pwid             0    15     0
#>  5 Otu_104 pwid             0    15     0
#>  6 Otu_105 hts              0    41     0
#>  7 Otu_105 msm              0   100     0
#>  8 Otu_108 pwid             0    15     0
#>  9 Otu_109 pwid             0    15     0
#> 10 Otu_11  pwid             0    15     0
#> # ℹ 1,343 more rows
```
