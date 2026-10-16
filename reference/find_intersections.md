# Finds common OTU between method results

Finds common OTU between method results

## Usage

``` r
find_intersections(
  rec,
  steps = steps_ids(rec, "da"),
  target = NULL,
  contrast_id = NULL
)
```

## Arguments

- rec:

  A Recipe object.

- steps:

  character vector with step ids to take in account

- target:

  Optional modeled target to include.

- contrast_id:

  Optional modeled contrast identifier to include.

## Value

tibble

## Examples

``` r
data(test_prep_rec)

## Extract all intersections
find_intersections(test_prep_rec)
#> # A tibble: 46 × 9
#>    taxa_id contrast_id          comparison   contrast_type var   effect step_ids
#>    <chr>   <glue>               <glue>       <chr>         <chr> <chr>  <chr>   
#>  1 Otu_119 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… down   deseq__…
#>  2 Otu_119 RiskGroup2[pwid-msm] RiskGroup2[… main          Risk… up     deseq__…
#>  3 Otu_125 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… down   deseq__…
#>  4 Otu_259 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  5 Otu_262 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  6 Otu_262 RiskGroup2[pwid-msm] RiskGroup2[… main          Risk… down   deseq__…
#>  7 Otu_332 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  8 Otu_332 RiskGroup2[pwid-msm] RiskGroup2[… main          Risk… down   deseq__…
#>  9 Otu_35  RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#> 10 Otu_354 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#> # ℹ 36 more rows
#> # ℹ 2 more variables: sum_methods <int>, taxa <chr>

## Exclude the first method
find_intersections(
  test_prep_rec,
  steps = steps_ids(test_prep_rec, "da")[-1]
)
#> # A tibble: 36 × 9
#>    taxa_id contrast_id          comparison   contrast_type var   effect step_ids
#>    <chr>   <glue>               <glue>       <chr>         <chr> <chr>  <chr>   
#>  1 Otu_1   RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  2 Otu_119 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… down   deseq__…
#>  3 Otu_119 RiskGroup2[pwid-msm] RiskGroup2[… main          Risk… up     deseq__…
#>  4 Otu_125 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… down   deseq__…
#>  5 Otu_231 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  6 Otu_255 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  7 Otu_259 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  8 Otu_262 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#>  9 Otu_262 RiskGroup2[pwid-msm] RiskGroup2[… main          Risk… down   deseq__…
#> 10 Otu_332 RiskGroup2[msm-hts]  RiskGroup2[… main          Risk… up     deseq__…
#> # ℹ 26 more rows
#> # ℹ 2 more variables: sum_methods <int>, taxa <chr>
```
