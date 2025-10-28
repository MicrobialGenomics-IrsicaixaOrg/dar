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
#> # A tibble: 88 × 4
#>    taxa_id taxa            step_ids                                  sum_methods
#>    <chr>   <chr>           <chr>                                           <dbl>
#>  1 Otu_34  Olsenella       maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  2 Otu_35  Collinsella     maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  3 Otu_45  Slackia         maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  4 Otu_78  Bacteroides     maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  5 Otu_82  Barnesiella     maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  6 Otu_88  Odoribacter     maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  7 Otu_94  Paraprevotella  maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  8 Otu_96  Prevotella      maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#>  9 Otu_119 Alistipes       maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#> 10 Otu_125 Parabacteroides maaslin__ChaSan, metagenomeseq__Zlebia, …           3
#> # ℹ 78 more rows

## Additionally, we can exclude some methods form the table
intersections <- find_intersections(
  test_prep_rec, 
  steps = steps_ids(test_prep_rec, "da")[-1]
)

intersections
#> # A tibble: 36 × 4
#>    taxa_id taxa            step_ids                                  sum_methods
#>    <chr>   <chr>           <chr>                                           <dbl>
#>  1 Otu_34  Olsenella       metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  2 Otu_35  Collinsella     metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  3 Otu_45  Slackia         metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  4 Otu_78  Bacteroides     metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  5 Otu_82  Barnesiella     metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  6 Otu_88  Odoribacter     metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  7 Otu_94  Paraprevotella  metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  8 Otu_96  Prevotella      metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#>  9 Otu_119 Alistipes       metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#> 10 Otu_125 Parabacteroides metagenomeseq__Zlebia, deseq__Linzer_tor…           2
#> # ℹ 26 more rows
```
