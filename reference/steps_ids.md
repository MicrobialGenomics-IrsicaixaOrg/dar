# Get step_ids from recipe

Get step_ids from recipe

## Usage

``` r
steps_ids(rec, type = "all")
```

## Arguments

- rec:

  A Recipe object.

- type:

  character vector indicating the type class. Options
  `c("all", "da", "prepro")`.

## Value

character vector

## Examples

``` r
data(test_rec)

## We can extract the step identifiers from a Recipe with `step_ids`
ids <- steps_ids(test_rec)
ids
#> [1] "subset_taxa__Boyoz"     "filter_taxa__Charlotte" "maaslin__ChaSan"       
#> [4] "metagenomeseq__Zlebia"  "deseq__Linzer_torte"   

## With the `type` parameter, extract the prepro and da steps separately.
da_ids <- steps_ids(test_rec, type = "da")
da_ids
#> [1] "maaslin__ChaSan"       "metagenomeseq__Zlebia" "deseq__Linzer_torte"  

prepro_ids <- steps_ids(test_rec, type = "prepro")
prepro_ids
#> [1] "subset_taxa__Boyoz"     "filter_taxa__Charlotte"
```
