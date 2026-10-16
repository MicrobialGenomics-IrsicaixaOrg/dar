# Extract results from defined bake

Extract results from defined bake

## Usage

``` r
cool(rec, bake = 1)
```

## Arguments

- rec:

  A `PrepRecipe` object.

- bake:

  Name or index of the bake to extract.

## Value

tbl_df

## Examples

``` r
data(test_prep_rec)
rec <- bake(test_prep_rec)
cool(rec)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 16 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  2 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  3 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  4 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                2
#>  5 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… up                2
#>  6 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… down              2
#>  7 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… up                2
#>  8 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… down              2
#>  9 Otu_35  Colli… RiskGroup2… RiskGroup… main          Risk… up                2
#> 10 Otu_354 Firmi… RiskGroup2… RiskGroup… main          Risk… up                2
#> 11 Otu_78  Bacte… RiskGroup2… RiskGroup… main          Risk… down              2
#> 12 Otu_78  Bacte… RiskGroup2… RiskGroup… main          Risk… up                2
#> 13 Otu_82  Barne… RiskGroup2… RiskGroup… main          Risk… down              2
#> 14 Otu_88  Odori… RiskGroup2… RiskGroup… main          Risk… down              2
#> 15 Otu_88  Odori… RiskGroup2… RiskGroup… main          Risk… up                2
#> 16 Otu_96  Prevo… RiskGroup2… RiskGroup… main          Risk… up                2
#> # ℹ 1 more variable: methods <chr>
```
