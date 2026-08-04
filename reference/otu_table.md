# Extract a canonical count table from a recipe

Extract a canonical count table from a recipe

## Usage

``` r
otu_table(rec)
```

## Arguments

- rec:

  A `Recipe` or `PrepRecipe` object.

## Value

A wide tibble with one row per taxon, `taxa_id` first, and one column
per sample. Taxa and samples follow
[`phyloseq::taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html)
and
[`phyloseq::sample_names()`](https://rdrr.io/pkg/phyloseq/man/sample_names-methods.html)
regardless of the stored OTU-table orientation.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy)
otu_table(rec)
#> # A tibble: 451 × 157
#>    taxa_id Sample_186 Sample_185 Sample_184 Sample_182 Sample_181 Sample_170
#>    <chr>        <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
#>  1 Otu_1         2538          0       9026      12049      22126      15374
#>  2 Otu_2            0          0          0        609       2543        988
#>  3 Otu_3            0          0          0          0          0          0
#>  4 Otu_4            0          0          0          0          0          0
#>  5 Otu_5            0          0          0          0          0          0
#>  6 Otu_6            0          0          0          0          0          0
#>  7 Otu_7            0          0          0          0          0          0
#>  8 Otu_8            0          0          0          0          0          0
#>  9 Otu_9            0          0          0          0          0        351
#> 10 Otu_10           0          0      10575      10764        504          0
#> # ℹ 441 more rows
#> # ℹ 150 more variables: Sample_169 <dbl>, Sample_168 <dbl>, Sample_167 <dbl>,
#> #   Sample_166 <dbl>, Sample_165 <dbl>, Sample_164 <dbl>, Sample_163 <dbl>,
#> #   Sample_162 <dbl>, Sample_160 <dbl>, Sample_159 <dbl>, Sample_158 <dbl>,
#> #   Sample_157 <dbl>, Sample_156 <dbl>, Sample_155 <dbl>, Sample_153 <dbl>,
#> #   Sample_152 <dbl>, Sample_151 <dbl>, Sample_150 <dbl>, Sample_149 <dbl>,
#> #   Sample_148 <dbl>, Sample_147 <dbl>, Sample_146 <dbl>, Sample_145 <dbl>, …
```
