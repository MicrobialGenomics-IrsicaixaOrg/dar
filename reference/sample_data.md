# Extracts sample_data from phyloseq inside a Recipe

Extracts sample_data from phyloseq inside a Recipe

## Usage

``` r
sample_data(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

A tibble

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
sample_data(rec)
#> # A tibble: 156 × 2
#>    sample_id  RiskGroup2
#>    <chr>      <chr>     
#>  1 Sample_186 msm       
#>  2 Sample_185 msm       
#>  3 Sample_184 hts       
#>  4 Sample_182 pwid      
#>  5 Sample_181 msm       
#>  6 Sample_170 msm       
#>  7 Sample_169 msm       
#>  8 Sample_168 pwid      
#>  9 Sample_167 hts       
#> 10 Sample_166 msm       
#> # ℹ 146 more rows
```
