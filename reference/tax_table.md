# Extract complete taxonomy from a recipe

Extract complete taxonomy from a recipe

## Usage

``` r
tax_table(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

A tibble with `taxa_id` followed by every taxonomic rank. Taxa are
ordered according to
[`phyloseq::taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html).
If the recipe has no taxonomy, an ID-only tibble is returned.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
tax_table(rec)
#> # A tibble: 451 × 8
#>    taxa_id Kingdom  Phylum         Class           Order    Family Genus Species
#>    <chr>   <chr>    <chr>          <chr>           <chr>    <chr>  <chr> <chr>  
#>  1 Otu_1   Archaea  Euryarchaeota  Methanobacteria Methano… Metha… Meth… Methan…
#>  2 Otu_2   Archaea  Euryarchaeota  Methanobacteria Methano… Metha… Meth… Methan…
#>  3 Otu_3   Bacteria Actinobacteria Actinobacteria  Actinom… Actin… Acti… Actino…
#>  4 Otu_4   Bacteria Actinobacteria Actinobacteria  Actinom… Actin… Acti… Actino…
#>  5 Otu_5   Bacteria Actinobacteria Actinobacteria  Actinom… Actin… Acti… Actino…
#>  6 Otu_6   Bacteria Actinobacteria Actinobacteria  Actinom… Actin… Acti… Actino…
#>  7 Otu_7   Bacteria Actinobacteria Actinobacteria  Actinom… Actin… Acti… Actino…
#>  8 Otu_8   Bacteria Actinobacteria Actinobacteria  Bifidob… Bifid… Aeri… Aerisc…
#>  9 Otu_9   Bacteria Actinobacteria Actinobacteria  Bifidob… Bifid… Allo… Allosc…
#> 10 Otu_10  Bacteria Actinobacteria Actinobacteria  Bifidob… Bifid… Bifi… Bifido…
#> # ℹ 441 more rows
```
