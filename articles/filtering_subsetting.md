# Filtering and Subsetting

The `step_filter_taxa` function is a general function that allows for
flexible filtering of OTUs based on across-sample abundance criteria.
The other functions, `step_filter_by_prevalence`,
`step_filter_by_variance`, `step_filter_by_abundance`, and
`step_filter_by_rarity`, are convenience wrappers around
`step_filter_taxa`, each designed to filter OTUs based on a specific
criterion: prevalence, variance, abundance, and rarity, respectively.

The `step_subset_taxa` function is used to subset taxa based on their
taxonomic level.

The phyloseq or TSE used as input can be pre-filtered using methods that
are most convenient to the user. However, the `dar` package provides
several functions to perform this filtering directly on the recipe
object.

## step_filter_taxa

The `step_filter_taxa` function applies an arbitrary set of functions to
OTUs as across-sample criteria. It takes a phyloseq object as input and
returns a logical vector indicating whether each OTU passed the
criteria. If the “prune” option is set to FALSE, it returns the
already-trimmed version of the phyloseq object.

``` r

library(dar)
data("metaHIV_phy")

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <- 
  step_filter_taxa(rec, .f = "function(x) sum(x > 0) >= (0 * length(x))") |> 
  prep()
```

## Convenience Wrappers

### step_filter_by_abundance

This function filters OTUs based on their abundance. The taxa retained
in the dataset are those where the sum of their abundance is greater
than the product of the total abundance and the provided threshold.

``` r

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <- 
  step_filter_by_abundance(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_prevalence

This function filters OTUs based on their prevalence. The taxa retained
in the dataset are those where the prevalence is greater than the
provided threshold.

``` r

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <- 
  step_filter_by_prevalence(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_rarity

This function filters OTUs based on their rarity. The taxa retained in
the dataset are those where the sum of their rarity is less than the
provided threshold.

``` r

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <- 
  step_filter_by_rarity(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_variance

This function filters OTUs based on their variance. The taxa retained in
the dataset are those where the variance of their abundance is greater
than the provided threshold.

``` r

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <- 
  step_filter_by_variance(rec, threshold = 0.01) |> 
  prep()
```

## subset_taxa

The `subset_taxa` function subsets taxa based on their taxonomic level.
The taxa retained in the dataset are those where the taxonomic level
matches the provided taxa.

``` r

rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec <-
  step_subset_taxa(rec, tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  prep()
```

## Conclusion

These functions provide a powerful and flexible way to filter and subset
OTUs in phyloseq objects contained within a recipe object, making it
easier to work with complex experimental data. By understanding how to
use these functions effectively, you can streamline your data analysis
workflow and focus on the aspects of your data that are most relevant to
your research questions. The `dar` package offers the added convenience
of performing these operations directly on the recipe object.

## Session info

``` r

devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.2 (2025-10-31)
#>  os       Ubuntu 24.04.3 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       UTC
#>  date     2026-02-06
#>  pandoc   3.8.2.1 @ /usr/bin/ (via rmarkdown)
#>  quarto   1.7.32 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package      * version  date (UTC) lib source
#>  ade4           1.7-23   2025-02-14 [1] RSPM (R 4.5.2)
#>  ape            5.8-1    2024-12-16 [1] RSPM (R 4.5.2)
#>  assertthat     0.2.1    2019-03-21 [1] RSPM (R 4.5.0)
#>  Biobase        2.70.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocGenerics   0.56.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  biomformat     1.38.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  Biostrings     2.78.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  bslib          0.10.0   2026-01-26 [2] RSPM (R 4.5.0)
#>  ca             0.71.1   2020-01-24 [1] RSPM (R 4.5.0)
#>  cachem         1.1.0    2024-05-16 [2] RSPM (R 4.5.0)
#>  cli            3.6.5    2025-04-23 [2] RSPM (R 4.5.0)
#>  cluster        2.1.8.1  2025-03-12 [3] CRAN (R 4.5.2)
#>  codetools      0.2-20   2024-03-31 [3] CRAN (R 4.5.2)
#>  crayon         1.5.3    2024-06-20 [2] RSPM (R 4.5.0)
#>  dar          * 1.5.6    2026-02-06 [1] Bioconductor
#>  data.table     1.18.2.1 2026-01-27 [1] RSPM (R 4.5.0)
#>  dendextend     1.19.1   2025-07-15 [1] RSPM (R 4.5.0)
#>  desc           1.4.3    2023-12-10 [2] RSPM (R 4.5.0)
#>  devtools       2.4.6    2025-10-03 [2] RSPM (R 4.5.0)
#>  digest         0.6.39   2025-11-19 [2] RSPM (R 4.5.0)
#>  dplyr          1.2.0    2026-02-03 [1] RSPM (R 4.5.0)
#>  ellipsis       0.3.2    2021-04-29 [2] RSPM (R 4.5.0)
#>  evaluate       1.0.5    2025-08-27 [2] RSPM (R 4.5.0)
#>  farver         2.1.2    2024-05-13 [1] RSPM (R 4.5.0)
#>  fastmap        1.2.0    2024-05-15 [2] RSPM (R 4.5.0)
#>  foreach        1.5.2    2022-02-02 [1] RSPM (R 4.5.0)
#>  fs             1.6.6    2025-04-12 [2] RSPM (R 4.5.0)
#>  furrr          0.3.1    2022-08-15 [1] RSPM (R 4.5.0)
#>  future         1.69.0   2026-01-16 [1] RSPM (R 4.5.0)
#>  generics       0.1.4    2025-05-09 [1] RSPM (R 4.5.0)
#>  ggplot2        4.0.2    2026-02-03 [1] RSPM (R 4.5.0)
#>  globals        0.19.0   2026-02-02 [1] RSPM (R 4.5.0)
#>  glue           1.8.0    2024-09-30 [2] RSPM (R 4.5.0)
#>  gridExtra      2.3      2017-09-09 [1] RSPM (R 4.5.0)
#>  gtable         0.3.6    2024-10-25 [1] RSPM (R 4.5.0)
#>  heatmaply      1.6.0    2025-07-12 [1] RSPM (R 4.5.0)
#>  htmltools      0.5.9    2025-12-04 [2] RSPM (R 4.5.0)
#>  htmlwidgets    1.6.4    2023-12-06 [2] RSPM (R 4.5.0)
#>  httr           1.4.7    2023-08-15 [1] RSPM (R 4.5.0)
#>  igraph         2.2.1    2025-10-27 [1] RSPM (R 4.5.0)
#>  IRanges        2.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  iterators      1.0.14   2022-02-05 [1] RSPM (R 4.5.0)
#>  jquerylib      0.1.4    2021-04-26 [2] RSPM (R 4.5.0)
#>  jsonlite       2.0.0    2025-03-27 [2] RSPM (R 4.5.0)
#>  knitr          1.51     2025-12-20 [2] RSPM (R 4.5.0)
#>  lattice        0.22-7   2025-04-02 [3] CRAN (R 4.5.2)
#>  lazyeval       0.2.2    2019-03-15 [1] RSPM (R 4.5.0)
#>  lifecycle      1.0.5    2026-01-08 [2] RSPM (R 4.5.0)
#>  listenv        0.10.0   2025-11-02 [1] RSPM (R 4.5.0)
#>  magrittr       2.0.4    2025-09-12 [2] RSPM (R 4.5.0)
#>  MASS           7.3-65   2025-02-28 [3] CRAN (R 4.5.2)
#>  Matrix         1.7-4    2025-08-28 [3] CRAN (R 4.5.2)
#>  memoise        2.0.1    2021-11-26 [2] RSPM (R 4.5.0)
#>  mgcv           1.9-4    2025-11-07 [3] RSPM (R 4.5.0)
#>  multtest       2.66.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  nlme           3.1-168  2025-03-31 [3] CRAN (R 4.5.2)
#>  otel           0.2.0    2025-08-29 [2] RSPM (R 4.5.0)
#>  parallelly     1.46.1   2026-01-08 [1] RSPM (R 4.5.0)
#>  permute        0.9-8    2025-06-25 [1] RSPM (R 4.5.0)
#>  phyloseq       1.54.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  pillar         1.11.1   2025-09-17 [2] RSPM (R 4.5.0)
#>  pkgbuild       1.4.8    2025-05-26 [2] RSPM (R 4.5.0)
#>  pkgconfig      2.0.3    2019-09-22 [2] RSPM (R 4.5.0)
#>  pkgdown        2.2.0    2025-11-06 [2] RSPM (R 4.5.0)
#>  pkgload        1.5.0    2026-02-03 [2] RSPM (R 4.5.0)
#>  plotly         4.12.0   2026-01-24 [1] RSPM (R 4.5.0)
#>  plyr           1.8.9    2023-10-02 [1] RSPM (R 4.5.2)
#>  purrr          1.2.1    2026-01-09 [2] RSPM (R 4.5.0)
#>  R6             2.6.1    2025-02-15 [2] RSPM (R 4.5.0)
#>  ragg           1.5.0    2025-09-02 [2] RSPM (R 4.5.0)
#>  RColorBrewer   1.1-3    2022-04-03 [1] RSPM (R 4.5.0)
#>  Rcpp           1.1.1    2026-01-10 [2] RSPM (R 4.5.0)
#>  registry       0.5-1    2019-03-05 [1] RSPM (R 4.5.0)
#>  remotes        2.5.0    2024-03-17 [1] RSPM (R 4.5.0)
#>  reshape2       1.4.5    2025-11-12 [1] RSPM (R 4.5.2)
#>  rhdf5          2.54.1   2025-12-04 [1] Bioconductor 3.22 (R 4.5.2)
#>  rhdf5filters   1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  Rhdf5lib       1.32.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  rlang          1.1.7    2026-01-09 [2] RSPM (R 4.5.0)
#>  rmarkdown      2.30     2025-09-28 [2] RSPM (R 4.5.0)
#>  S4Vectors      0.48.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  S7             0.2.1    2025-11-14 [1] RSPM (R 4.5.0)
#>  sass           0.4.10   2025-04-11 [2] RSPM (R 4.5.0)
#>  scales         1.4.0    2025-04-24 [1] RSPM (R 4.5.0)
#>  Seqinfo        1.0.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  seriation      1.5.8    2025-08-20 [1] RSPM (R 4.5.0)
#>  sessioninfo    1.2.3    2025-02-05 [2] RSPM (R 4.5.0)
#>  stringi        1.8.7    2025-03-27 [2] RSPM (R 4.5.0)
#>  stringr        1.6.0    2025-11-04 [2] RSPM (R 4.5.0)
#>  survival       3.8-6    2026-01-16 [3] RSPM (R 4.5.0)
#>  systemfonts    1.3.1    2025-10-01 [2] RSPM (R 4.5.0)
#>  textshaping    1.0.4    2025-10-10 [2] RSPM (R 4.5.0)
#>  tibble         3.3.1    2026-01-11 [2] RSPM (R 4.5.0)
#>  tidyr          1.3.2    2025-12-19 [1] RSPM (R 4.5.0)
#>  tidyselect     1.2.1    2024-03-11 [1] RSPM (R 4.5.0)
#>  TSP            1.2.6    2025-11-27 [1] RSPM (R 4.5.0)
#>  UpSetR         1.4.0    2019-05-22 [1] RSPM (R 4.5.0)
#>  usethis        3.2.1    2025-09-06 [2] RSPM (R 4.5.0)
#>  vctrs          0.7.1    2026-01-23 [2] RSPM (R 4.5.0)
#>  vegan          2.7-2    2025-10-08 [1] RSPM (R 4.5.0)
#>  viridis        0.6.5    2024-01-29 [1] RSPM (R 4.5.0)
#>  viridisLite    0.4.3    2026-02-04 [1] RSPM (R 4.5.0)
#>  webshot        0.5.5    2023-06-26 [1] RSPM (R 4.5.0)
#>  xfun           0.56     2026-01-18 [2] RSPM (R 4.5.0)
#>  XVector        0.50.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  yaml           2.3.12   2025-12-10 [2] RSPM (R 4.5.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
