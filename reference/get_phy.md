# Returns phyloseq from Recipe-class object

Returns phyloseq from Recipe-class object

## Usage

``` r
get_phy(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

Phyloseq class object

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
get_phy(rec)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 451 taxa and 156 samples ]
#> sample_data() Sample Data:       [ 156 samples by 3 sample variables ]
#> tax_table()   Taxonomy Table:    [ 451 taxa by 7 taxonomic ranks ]
```
