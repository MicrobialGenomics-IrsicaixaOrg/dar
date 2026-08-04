# Collect differential-abundance results in a tidy contract

`tidy_results()` combines the results of executed differential-abundance
steps without changing the engine-specific tables stored in
`rec@results`. Effect sizes are labelled with their engine-specific
metric and must not be compared directly across incompatible metrics.

## Usage

``` r
tidy_results(rec, steps = steps_ids(rec, "da"), significant_only = FALSE)
```

## Arguments

- rec:

  A
  [PrepRecipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
  object.

- steps:

  Character vector of executed differential-abundance step IDs to
  collect. By default, all executed DA steps are included.

- significant_only:

  If `TRUE`, retain only rows classified as significant by their
  originating step.

## Value

A deterministic tibble with one row per taxon, contrast and configured
step. It contains taxon identifiers (`taxa_id`, `taxa`), contrast
metadata (`contrast_id`, `comparison`, `contrast_type`, `var`,
`numerator`, `denominator`, `at_var`, `at_level`), statistics
(`effect_size`, `effect_metric`, `p_value`, `adj_p_value`,
`significant`), the stable engine in `method`, and the configured
instance in `step_id`.

## Examples

``` r
data(test_prep_rec)
results <- tidy_results(test_prep_rec)
results
#> # A tibble: 448 × 17
#>    taxa_id taxa             contrast_id comparison contrast_type var   numerator
#>    <chr>   <chr>            <chr>       <chr>      <chr>         <chr> <chr>    
#>  1 Otu_1   Methanobrevibac… RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  2 Otu_10  Bifidobacterium  RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  3 Otu_18  Bifidobacterium  RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  4 Otu_34  Olsenella        RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  5 Otu_35  Collinsella      RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  6 Otu_37  Collinsella      RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  7 Otu_39  Enorma           RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  8 Otu_45  Slackia          RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#>  9 Otu_46  Bacteroides      RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#> 10 Otu_47  Bacteroides      RiskGroup2… hts_msm    legacy_pairw… Risk… NA       
#> # ℹ 438 more rows
#> # ℹ 10 more variables: denominator <chr>, at_var <chr>, at_level <chr>,
#> #   effect_size <dbl>, effect_metric <chr>, p_value <dbl>, adj_p_value <dbl>,
#> #   significant <lgl>, method <chr>, step_id <chr>

tidy_results(test_prep_rec, significant_only = TRUE)
#> # A tibble: 101 × 17
#>    taxa_id taxa        contrast_id      comparison contrast_type var   numerator
#>    <chr>   <chr>       <chr>            <chr>      <chr>         <chr> <chr>    
#>  1 Otu_35  Collinsella RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  2 Otu_37  Collinsella RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  3 Otu_46  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  4 Otu_47  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  5 Otu_51  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  6 Otu_52  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  7 Otu_62  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  8 Otu_69  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#>  9 Otu_76  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#> 10 Otu_77  Bacteroides RiskGroup2[hts_… hts_msm    legacy_pairw… Risk… NA       
#> # ℹ 91 more rows
#> # ℹ 10 more variables: denominator <chr>, at_var <chr>, at_level <chr>,
#> #   effect_size <dbl>, effect_metric <chr>, p_value <dbl>, adj_p_value <dbl>,
#> #   significant <lgl>, method <chr>, step_id <chr>
```
