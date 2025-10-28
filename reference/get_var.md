# Returns var_info from Recipe-class object

Returns var_info from Recipe-class object

## Usage

``` r
get_var(rec)

# S4 method for class 'Recipe'
get_var(rec)
```

## Arguments

- rec:

  A `Recipe` object

## Value

Tibble containing `var_info`.

## Examples

``` r
data(metaHIV_phy)

## Define recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")

## Extract variable of interest
get_var(rec)
#> # A tibble: 1 × 1
#>   vars      
#>   <chr>     
#> 1 RiskGroup2
```
