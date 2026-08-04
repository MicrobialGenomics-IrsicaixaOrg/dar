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
#> # A tibble: 16 × 2
#>    taxa_id taxa            
#>    <chr>   <chr>           
#>  1 Otu_119 Alistipes       
#>  2 Otu_125 Parabacteroides 
#>  3 Otu_255 Blautia         
#>  4 Otu_259 Coprococcus     
#>  5 Otu_262 Dorea           
#>  6 Otu_307 Faecalibacterium
#>  7 Otu_330 Ruthenibacterium
#>  8 Otu_332 Catenibacterium 
#>  9 Otu_35  Collinsella     
#> 10 Otu_365 Mitsuokella     
#> 11 Otu_433 Haemophilus     
#> 12 Otu_78  Bacteroides     
#> 13 Otu_82  Barnesiella     
#> 14 Otu_88  Odoribacter     
#> 15 Otu_94  Paraprevotella  
#> 16 Otu_96  Prevotella      
```
