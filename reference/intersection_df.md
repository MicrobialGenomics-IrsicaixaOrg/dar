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
#>   taxa_id          contrast_id           comparison contrast_type        var
#> 1   Otu_1  RiskGroup2[msm-hts]  RiskGroup2[msm-hts]          main RiskGroup2
#> 2   Otu_1 RiskGroup2[pwid-hts] RiskGroup2[pwid-hts]          main RiskGroup2
#> 3   Otu_1 RiskGroup2[pwid-msm] RiskGroup2[pwid-msm]          main RiskGroup2
#> 4  Otu_10  RiskGroup2[msm-hts]  RiskGroup2[msm-hts]          main RiskGroup2
#> 5  Otu_10 RiskGroup2[pwid-hts] RiskGroup2[pwid-hts]          main RiskGroup2
#> 6  Otu_10 RiskGroup2[pwid-msm] RiskGroup2[pwid-msm]          main RiskGroup2
#>   effect maaslin__Welsh_cake deseq__Coussin_de_Lyon
#> 1     up                   0                      1
#> 2     up                   0                      0
#> 3   down                   0                      0
#> 4     up                   0                      0
#> 5     up                   0                      0
#> 6     up                   0                      0
```
