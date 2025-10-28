# Returns tax_info from Recipe-class object

Returns tax_info from Recipe-class object

## Usage

``` r
get_tax(rec)

# S4 method for class 'Recipe'
get_tax(rec)
```

## Arguments

- rec:

  A `Recipe` object

## Value

Tibble containing `tax_info`.

## Examples

``` r
data(metaHIV_phy)

## Define recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")

## Extract taxonomic level
get_tax(rec)
#> # A tibble: 1 × 1
#>   tax_lev
#>   <chr>  
#> 1 Species
```
