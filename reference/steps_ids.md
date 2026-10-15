# Get step_ids from recipe

Get step_ids from recipe

## Usage

``` r
steps_ids(rec, type = "all", include_skipped = FALSE)
```

## Arguments

- rec:

  A Recipe object.

- type:

  character vector indicating the type class. Options
  `c("all", "da", "prepro")`.

- include_skipped:

  For a prepared modeled recipe, include DA steps that were declared but
  skipped as statistically incompatible.

## Value

character vector

## Examples

``` r
data(test_rec)

steps_ids(test_rec)
#> [1] "subset_taxa__Viennoiserie" "filter_taxa__Karakudamono"
#> [3] "maaslin__Welsh_cake"       "deseq__Coussin_de_Lyon"   

## Select differential-abundance or preprocessing steps
steps_ids(test_rec, type = "da")
#> [1] "maaslin__Welsh_cake"    "deseq__Coussin_de_Lyon"
steps_ids(test_rec, type = "prepro")
#> [1] "subset_taxa__Viennoiserie" "filter_taxa__Karakudamono"
```
