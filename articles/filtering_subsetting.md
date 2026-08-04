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

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
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

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec <- 
  step_filter_by_abundance(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_prevalence

This function filters OTUs based on their prevalence. The taxa retained
in the dataset are those where the prevalence is greater than the
provided threshold.

``` r

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec <- 
  step_filter_by_prevalence(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_rarity

This function filters OTUs based on their rarity. The taxa retained in
the dataset are those where the sum of their rarity is less than the
provided threshold.

``` r

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec <- 
  step_filter_by_rarity(rec, threshold = 0.01) |> 
  prep()
```

### step_filter_by_variance

This function filters OTUs based on their variance. The taxa retained in
the dataset are those where the variance of their abundance is greater
than the provided threshold.

``` r

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec <- 
  step_filter_by_variance(rec, threshold = 0.01) |> 
  prep()
```

## subset_taxa

The `subset_taxa` function subsets taxa based on their taxonomic level.
The taxa retained in the dataset are those where the taxonomic level
matches the provided taxa.

``` r

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
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
#>  version  R version 4.6.1 (2026-06-24)
#>  os       Ubuntu 24.04.4 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       UTC
#>  date     2026-08-04
#>  pandoc   3.10 @ /usr/bin/ (via rmarkdown)
#>  quarto   1.9.38 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package                  * version date (UTC) lib source
#>  abind                      1.4-8   2024-09-12 [1] RSPM (R 4.6.0)
#>  ade4                       1.7-24  2026-03-21 [1] RSPM (R 4.6.0)
#>  ape                        5.8-1   2024-12-16 [1] RSPM (R 4.6.0)
#>  backports                  1.5.1   2026-04-03 [1] RSPM (R 4.6.0)
#>  beachmat                   2.29.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  beeswarm                   0.4.0   2021-06-01 [1] RSPM (R 4.6.0)
#>  Biobase                    2.73.2  2026-07-29 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocGenerics               0.59.10 2026-07-07 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocNeighbors              2.7.2   2026-05-10 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocParallel               1.47.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocSingular               1.29.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  biomformat                 1.41.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  Biostrings                 2.81.6  2026-07-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  bluster                    1.23.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  bslib                      0.12.0  2026-08-04 [2] RSPM (R 4.6.0)
#>  cachem                     1.1.0   2024-05-16 [2] RSPM (R 4.6.0)
#>  checkmate                  2.3.4   2026-02-03 [1] RSPM (R 4.6.0)
#>  cli                        3.6.6   2026-04-09 [2] RSPM (R 4.6.0)
#>  cluster                    2.1.8.3 2026-07-30 [2] RSPM (R 4.6.0)
#>  codetools                  0.2-20  2024-03-31 [3] CRAN (R 4.6.1)
#>  crayon                     1.5.3   2024-06-20 [2] RSPM (R 4.6.0)
#>  dar                      * 1.9.12  2026-08-04 [1] Bioconductor
#>  data.table                 1.18.4  2026-05-06 [1] RSPM (R 4.6.0)
#>  DBI                        1.3.0   2026-02-25 [1] RSPM (R 4.6.0)
#>  DECIPHER                   3.9.2   2026-07-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  decontam                   1.33.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  DelayedArray               0.39.3  2026-06-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  DelayedMatrixStats         1.35.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  desc                       1.4.3   2023-12-10 [2] RSPM (R 4.6.0)
#>  devtools                   2.5.2   2026-04-30 [2] RSPM (R 4.6.0)
#>  digest                     0.6.39  2025-11-19 [2] RSPM (R 4.6.0)
#>  DirichletMultinomial       1.55.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  dplyr                      1.2.1   2026-04-03 [1] RSPM (R 4.6.0)
#>  ellipsis                   0.3.3   2026-04-04 [2] RSPM (R 4.6.0)
#>  evaluate                   1.0.5   2025-08-27 [2] RSPM (R 4.6.0)
#>  farver                     2.1.2   2024-05-13 [1] RSPM (R 4.6.0)
#>  fastmap                    1.2.0   2024-05-15 [2] RSPM (R 4.6.0)
#>  foreach                    1.5.2   2022-02-02 [1] RSPM (R 4.6.0)
#>  fs                         2.1.0   2026-04-18 [2] RSPM (R 4.6.0)
#>  furrr                      0.4.0   2026-03-31 [1] RSPM (R 4.6.0)
#>  future                     1.75.0  2026-07-20 [1] RSPM (R 4.6.0)
#>  generics                   0.1.4   2025-05-09 [1] RSPM (R 4.6.0)
#>  GenomicRanges              1.65.1  2026-07-08 [1] Bioconductor 3.24 (R 4.6.1)
#>  ggbeeswarm                 0.7.3   2025-11-29 [1] RSPM (R 4.6.0)
#>  ggplot2                    4.0.3   2026-04-22 [1] RSPM (R 4.6.0)
#>  ggrepel                    0.9.8   2026-03-17 [1] RSPM (R 4.6.0)
#>  globals                    0.19.1  2026-03-13 [1] RSPM (R 4.6.0)
#>  glue                       1.8.1   2026-04-17 [2] RSPM (R 4.6.0)
#>  gridExtra                  2.3.1   2026-06-25 [1] RSPM (R 4.6.0)
#>  gtable                     0.3.6   2024-10-25 [1] RSPM (R 4.6.0)
#>  htmltools                  0.5.9   2025-12-04 [2] RSPM (R 4.6.0)
#>  htmlwidgets                1.6.4   2023-12-06 [2] RSPM (R 4.6.0)
#>  igraph                     2.3.3   2026-06-26 [1] RSPM (R 4.6.0)
#>  IRanges                    2.47.2  2026-06-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  irlba                      2.3.7   2026-01-30 [1] RSPM (R 4.6.0)
#>  iterators                  1.0.14  2022-02-05 [1] RSPM (R 4.6.0)
#>  jquerylib                  0.1.4   2021-04-26 [2] RSPM (R 4.6.0)
#>  jsonlite                   2.0.0   2025-03-27 [2] RSPM (R 4.6.0)
#>  knitr                      1.51    2025-12-20 [2] RSPM (R 4.6.0)
#>  lattice                    0.22-9  2026-02-09 [3] CRAN (R 4.6.1)
#>  lazyeval                   0.2.3   2026-04-04 [1] RSPM (R 4.6.0)
#>  lifecycle                  1.0.5   2026-01-08 [2] RSPM (R 4.6.0)
#>  listenv                    1.0.0   2026-06-22 [1] RSPM (R 4.6.0)
#>  magrittr                   2.0.5   2026-04-04 [2] RSPM (R 4.6.0)
#>  MASS                       7.3-66  2026-07-15 [2] RSPM (R 4.6.0)
#>  Matrix                     1.7-6   2026-07-25 [2] RSPM (R 4.6.0)
#>  MatrixGenerics             1.25.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  matrixStats                1.5.0   2025-01-07 [1] RSPM (R 4.6.0)
#>  memoise                    2.0.1   2021-11-26 [2] RSPM (R 4.6.0)
#>  mgcv                       1.9-4   2025-11-07 [3] CRAN (R 4.6.1)
#>  mia                        1.21.6  2026-07-31 [1] Bioconductor 3.24 (R 4.6.1)
#>  MultiAssayExperiment       1.39.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  multtest                   2.69.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  nlme                       3.1-170 2026-07-15 [2] RSPM (R 4.6.0)
#>  otel                       0.2.0   2025-08-29 [2] RSPM (R 4.6.0)
#>  parallelly                 1.48.0  2026-06-29 [1] RSPM (R 4.6.0)
#>  permute                    0.9-10  2026-02-06 [1] RSPM (R 4.6.0)
#>  phyloseq                   1.57.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  pillar                     1.11.1  2025-09-17 [2] RSPM (R 4.6.0)
#>  pkgbuild                   1.4.8   2025-05-26 [2] RSPM (R 4.6.0)
#>  pkgconfig                  2.0.3   2019-09-22 [2] RSPM (R 4.6.0)
#>  pkgdown                    2.2.1   2026-07-07 [2] RSPM (R 4.6.0)
#>  pkgload                    1.5.3   2026-06-15 [2] RSPM (R 4.6.0)
#>  plyr                       1.8.9   2023-10-02 [1] RSPM (R 4.6.0)
#>  purrr                      1.2.2   2026-04-10 [2] RSPM (R 4.6.0)
#>  R6                         2.6.1   2025-02-15 [2] RSPM (R 4.6.0)
#>  ragg                       1.5.2   2026-03-23 [2] RSPM (R 4.6.0)
#>  rappdirs                   0.3.4   2026-01-17 [2] RSPM (R 4.6.0)
#>  rbibutils                  2.4.1   2026-01-21 [1] RSPM (R 4.6.0)
#>  RColorBrewer               1.1-3   2022-04-03 [1] RSPM (R 4.6.0)
#>  Rcpp                       1.1.2   2026-07-05 [2] RSPM (R 4.6.0)
#>  Rdpack                     2.6.6   2026-02-08 [1] RSPM (R 4.6.0)
#>  reformulas                 0.4.4   2026-02-02 [1] RSPM (R 4.6.0)
#>  reshape2                   1.4.5   2025-11-12 [1] RSPM (R 4.6.0)
#>  rlang                      1.3.0   2026-07-05 [2] RSPM (R 4.6.0)
#>  rmarkdown                  2.31    2026-03-26 [2] RSPM (R 4.6.0)
#>  rsvd                       1.0.5   2021-04-16 [1] RSPM (R 4.6.0)
#>  S4Arrays                   1.13.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  S4Vectors                  0.51.6  2026-07-29 [1] Bioconductor 3.24 (R 4.6.1)
#>  S7                         0.2.2   2026-04-22 [1] CRAN (R 4.6.1)
#>  sass                       0.4.10  2025-04-11 [2] RSPM (R 4.6.0)
#>  ScaledMatrix               1.21.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  scales                     1.4.0   2025-04-24 [1] RSPM (R 4.6.0)
#>  scater                     1.41.2  2026-06-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  scuttle                    1.23.1  2026-05-18 [1] Bioconductor 3.24 (R 4.6.1)
#>  Seqinfo                    1.3.0   2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  sessioninfo                1.2.4   2026-06-04 [2] RSPM (R 4.6.0)
#>  SingleCellExperiment       1.35.2  2026-07-16 [1] Bioconductor 3.24 (R 4.6.1)
#>  SparseArray                1.13.2  2026-05-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  sparseMatrixStats          1.25.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  stringi                    1.8.7   2025-03-27 [2] RSPM (R 4.6.0)
#>  stringr                    1.6.0   2025-11-04 [1] RSPM (R 4.6.0)
#>  SummarizedExperiment       1.43.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  survival                   3.8-9   2026-07-08 [2] RSPM (R 4.6.0)
#>  systemfonts                1.3.2   2026-03-05 [2] RSPM (R 4.6.0)
#>  textshaping                1.0.5   2026-03-06 [2] RSPM (R 4.6.0)
#>  tibble                     3.3.1   2026-01-11 [2] RSPM (R 4.6.0)
#>  tidyr                      1.3.2   2025-12-19 [1] RSPM (R 4.6.0)
#>  tidyselect                 1.2.1   2024-03-11 [1] RSPM (R 4.6.0)
#>  tidytree                   0.4.8   2026-07-02 [1] RSPM (R 4.6.0)
#>  treeio                     1.37.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  TreeSummarizedExperiment   2.21.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  usethis                    3.2.1   2025-09-06 [2] RSPM (R 4.6.0)
#>  vctrs                      0.7.3   2026-04-11 [2] RSPM (R 4.6.0)
#>  vegan                      2.7-5   2026-05-25 [1] RSPM (R 4.6.0)
#>  vipor                      0.4.7   2023-12-18 [1] RSPM (R 4.6.0)
#>  viridis                    0.6.5   2024-01-29 [1] RSPM (R 4.6.0)
#>  viridisLite                0.4.3   2026-02-04 [1] RSPM (R 4.6.0)
#>  xfun                       0.60    2026-07-09 [2] RSPM (R 4.6.0)
#>  XVector                    0.53.0  2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  yaml                       2.3.12  2025-12-10 [2] RSPM (R 4.6.0)
#>  yulab.utils                0.2.4   2026-02-02 [1] RSPM (R 4.6.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
