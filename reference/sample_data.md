# Extracts sample_data from phyloseq inside a Recipe

Extracts sample_data from phyloseq inside a Recipe

## Usage

``` r
sample_data(rec)

# S4 method for class 'Recipe'
sample_data(rec)
```

## Arguments

- rec:

  A Recipe or Recipe step.

## Value

A tibble

## Examples

``` r
data(metaHIV_phy)

## Define recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")

## Extract sample_data from phyloseq object
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
