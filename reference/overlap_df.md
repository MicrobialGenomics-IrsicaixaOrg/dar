# Overlap of significant OTUs between tested methods.

Overlap of significant OTUs between tested methods.

## Usage

``` r
overlap_df(rec, steps = steps_ids(rec, "da"), type = "all")
```

## Arguments

- rec:

  A `PrepRecipe` object.

- steps:

  Character vector with step_ids to take in account.

- type:

  Indicates whether to use all taxa ("all") or only those that are
  differentially abundant in at least one method ("da"). Default as
  "all".

## Value

df

## Examples

``` r
data(test_prep_rec)
df <- overlap_df(test_prep_rec, steps_ids(test_prep_rec, "da"))
head(df)
#>                        maaslin__Welsh_cake deseq__Coussin_de_Lyon
#> maaslin__Welsh_cake              1.0000000              0.6435644
#> deseq__Coussin_de_Lyon           0.6435644              1.0000000
```
