# Extract complete sample metadata from a recipe

Extract complete sample metadata from a recipe

## Usage

``` r
sample_data(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

A tibble with `sample_id` followed by every sample metadata column.
Samples are ordered according to
[`phyloseq::sample_names()`](https://rdrr.io/pkg/phyloseq/man/sample_names-methods.html).
If the recipe has no sample metadata, an ID-only tibble is returned.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
sample_data(rec)
#> # A tibble: 156 × 4
#>    sample_id  RiskGroup2 Cluster       HIV_Status
#>    <chr>      <chr>      <chr>         <chr>     
#>  1 Sample_186 msm        "Prevotella"  positive  
#>  2 Sample_185 msm        ""            positive  
#>  3 Sample_184 hts        "Bacteroides" positive  
#>  4 Sample_182 pwid       "Bacteroides" positive  
#>  5 Sample_181 msm        "Prevotella"  positive  
#>  6 Sample_170 msm        "Prevotella"  positive  
#>  7 Sample_169 msm        "Prevotella"  negative  
#>  8 Sample_168 pwid       "Bacteroides" positive  
#>  9 Sample_167 hts        "Bacteroides" negative  
#> 10 Sample_166 msm        "Prevotella"  positive  
#> # ℹ 146 more rows
```
