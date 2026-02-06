# Finds common OTU between method results

Finds common OTU between method results

## Usage

``` r
find_intersections(rec, steps = steps_ids(rec, "da"))
```

## Arguments

- rec:

  A Recipe object.

- steps:

  character vector with step ids to take in account

## Value

tibble

## Examples

``` r
data(test_prep_rec)

## From a PrepRecipe we can extract a tibble with all intersections
intersections <- find_intersections(test_prep_rec)
intersections
#> # A tibble: 52 × 4
#>    taxa_id taxa            step_ids                             sum_methods
#>    <chr>   <chr>           <chr>                                      <dbl>
#>  1 Otu_35  Collinsella     maaslin__Eccles_cake, deseq__Belekoy           2
#>  2 Otu_78  Bacteroides     maaslin__Eccles_cake, deseq__Belekoy           2
#>  3 Otu_82  Barnesiella     maaslin__Eccles_cake, deseq__Belekoy           2
#>  4 Otu_88  Odoribacter     maaslin__Eccles_cake, deseq__Belekoy           2
#>  5 Otu_94  Paraprevotella  maaslin__Eccles_cake, deseq__Belekoy           2
#>  6 Otu_96  Prevotella      maaslin__Eccles_cake, deseq__Belekoy           2
#>  7 Otu_119 Alistipes       maaslin__Eccles_cake, deseq__Belekoy           2
#>  8 Otu_125 Parabacteroides maaslin__Eccles_cake, deseq__Belekoy           2
#>  9 Otu_255 Blautia         maaslin__Eccles_cake, deseq__Belekoy           2
#> 10 Otu_259 Coprococcus     maaslin__Eccles_cake, deseq__Belekoy           2
#> # ℹ 42 more rows

## Additionally, we can exclude some methods form the table
intersections <- find_intersections(
  test_prep_rec, 
  steps = steps_ids(test_prep_rec, "da")[-1]
)

intersections
#> # A tibble: 27 × 4
#>    taxa_id taxa            step_ids       sum_methods
#>    <chr>   <chr>           <chr>                <dbl>
#>  1 Otu_34  Olsenella       deseq__Belekoy           1
#>  2 Otu_35  Collinsella     deseq__Belekoy           1
#>  3 Otu_45  Slackia         deseq__Belekoy           1
#>  4 Otu_78  Bacteroides     deseq__Belekoy           1
#>  5 Otu_82  Barnesiella     deseq__Belekoy           1
#>  6 Otu_88  Odoribacter     deseq__Belekoy           1
#>  7 Otu_94  Paraprevotella  deseq__Belekoy           1
#>  8 Otu_96  Prevotella      deseq__Belekoy           1
#>  9 Otu_119 Alistipes       deseq__Belekoy           1
#> 10 Otu_125 Parabacteroides deseq__Belekoy           1
#> # ℹ 17 more rows
```
