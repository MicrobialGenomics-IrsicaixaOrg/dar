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
#> [1] "subset_taxa__Bear_claw"   "filter_taxa__Spanakopita"
#> [3] "maaslin__Eccles_cake"     "deseq__Belekoy"          

## With the `type` parameter, extract the prepro and da steps separately.
da_ids <- steps_ids(test_rec, type = "da")
da_ids
#> [1] "maaslin__Eccles_cake" "deseq__Belekoy"      

prepro_ids <- steps_ids(test_rec, type = "prepro")
prepro_ids
#> [1] "subset_taxa__Bear_claw"   "filter_taxa__Spanakopita"
```
