# Extracts tax_table from phyloseq inside a Recipe

Extracts tax_table from phyloseq inside a Recipe

## Usage

``` r
tax_table(rec)

# S4 method for class 'Recipe'
tax_table(rec)
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

## Extract tax_table from phyloseq object
tax_table(rec)
#> # A tibble: 451 × 2
#>    taxa_id taxa                         
#>    <chr>   <chr>                        
#>  1 Otu_1   Methanobrevibacter_smithii   
#>  2 Otu_2   Methanosphaera_stadtmanae    
#>  3 Otu_3   Actinomyces_odontolyticus    
#>  4 Otu_4   Actinomyces_radicidentis     
#>  5 Otu_5   Actinomyces_sp_HPA0247       
#>  6 Otu_6   Actinomyces_sp_oral_taxon_181
#>  7 Otu_7   Actinomyces_turicensis       
#>  8 Otu_8   Aeriscardovia_aeriphila      
#>  9 Otu_9   Alloscardovia_omnicolens     
#> 10 Otu_10  Bifidobacterium_adolescentis 
#> # ℹ 441 more rows
```
