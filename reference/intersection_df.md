# Returns data.frame with OTU intersection between methods

Returns data.frame with OTU intersection between methods

## Usage

``` r
intersection_df(
  rec,
  steps = steps_ids(rec, "da"),
  tidy = FALSE,
  target = NULL,
  contrast_id = NULL
)
```

## Arguments

- rec:

  A `PrepRecipe` object.

- steps:

  character vector with step_ids to take in account.

- tidy:

  Boolean indicating if result must be in tidy format.

- target:

  Optional modeled target to include.

- contrast_id:

  Optional modeled contrast identifier to include.

## Value

data.frame class object

## Examples

``` r
data(test_prep_rec)
df <- intersection_df(test_prep_rec)
head(df)
#>   taxa_id maaslin__Welsh_cake deseq__Coussin_de_Lyon
#> 1   Otu_1                   0                      0
#> 2  Otu_10                   0                      0
#> 3  Otu_18                   0                      0
#> 4  Otu_34                   0                      1
#> 5  Otu_35                   1                      1
#> 6  Otu_37                   1                      0
```
