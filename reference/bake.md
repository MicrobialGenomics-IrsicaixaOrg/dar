# Define consensus strategies from a Recipe

For a prep Recipe adds a consensus strategy to use for result
extraction.

## Usage

``` r
bake(
  rec,
  count_cutoff = NULL,
  weights = NULL,
  exclude = NULL,
  id = rand_id("bake")
)
```

## Arguments

- rec:

  A PrepRecipe object. The step will be added to the sequence of
  operations for this Recipe.

- count_cutoff:

  Indicates the minimum number of methods in which an OTU must be
  present (Default: NULL). If `count_cutoff` is `NULL`, it is set to
  `length(steps_ids(rec, "da")) - length(exclude)`.

- weights:

  Named numeric vector with the weight for each method.

- exclude:

  Method ids to exclude.

- id:

  A character string that is unique to this step to identify it.

## Value

An object of class `PrepRecipe`

## Examples

``` r
data(test_prep_rec)

## Require agreement among all differential-abundance methods
cool(bake(test_prep_rec))
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

## Retain taxa detected by at least one method
cool(bake(test_prep_rec, count_cutoff = 1))
#> ℹ Baking with count_cutoff = 1
#> # A tibble: 52 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_102 Prevotella     
#>  2 Otu_107 Prevotella     
#>  3 Otu_115 Alistipes      
#>  4 Otu_119 Alistipes      
#>  5 Otu_120 Alistipes      
#>  6 Otu_125 Parabacteroides
#>  7 Otu_129 Parabacteroides
#>  8 Otu_231 Eubacterium    
#>  9 Otu_234 Eubacterium    
#> 10 Otu_255 Blautia        
#> # ℹ 42 more rows

## Exclude a method from the consensus
method_ids <- steps_ids(test_prep_rec, "da")
cool(bake(test_prep_rec, exclude = method_ids[1]))
#> ℹ Baking with count_cutoff = 1
#> Excluding "maaslin__Welsh_cake".
#> # A tibble: 27 × 2
#>    taxa_id taxa            
#>    <chr>   <chr>           
#>  1 Otu_119 Alistipes       
#>  2 Otu_125 Parabacteroides 
#>  3 Otu_231 Eubacterium     
#>  4 Otu_255 Blautia         
#>  5 Otu_259 Coprococcus     
#>  6 Otu_262 Dorea           
#>  7 Otu_285 Roseburia       
#>  8 Otu_307 Faecalibacterium
#>  9 Otu_308 Flavonifractor  
#> 10 Otu_330 Ruthenibacterium
#> # ℹ 17 more rows

## Give the first method twice the weight of the others
weights <- setNames(rep(1, length(method_ids)), method_ids)
weights[1] <- 2
cool(bake(test_prep_rec, weights = weights))
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 41 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_102 Prevotella     
#>  2 Otu_107 Prevotella     
#>  3 Otu_115 Alistipes      
#>  4 Otu_119 Alistipes      
#>  5 Otu_120 Alistipes      
#>  6 Otu_125 Parabacteroides
#>  7 Otu_129 Parabacteroides
#>  8 Otu_234 Eubacterium    
#>  9 Otu_255 Blautia        
#> 10 Otu_258 Coprococcus    
#> # ℹ 31 more rows
```
