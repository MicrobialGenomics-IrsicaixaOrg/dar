# Introduction to dar

## An Example

The package includes a dataset from a study by [Noguera-Julian, M., et
al. 2016](https://doi.org/10.1016/j.ebiom.2016.01.032), which
investigated the differential abundance of microbial species between men
who have sex with men (MSM) and non-MSM (hts). This data is stored as an
object of the `phyloseq` class, a standard input format for creating
recipes with dar in conjunction with `TreeSummarizedExperiment`. To
begin the analysis, we first load and inspect the data:

``` r

library(dar)

data("metaHIV_phy", package = "dar")

metaHIV_phy
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 451 taxa and 156 samples ]
#> sample_data() Sample Data:       [ 156 samples by 3 sample variables ]
#> tax_table()   Taxonomy Table:    [ 451 taxa by 7 taxonomic ranks ]
```

## An Initial Recipe

First, we will create a recipe object from the original data and then
specify the processing and differential analysis steps.

Recipes can be created manually by sequentially adding roles to
variables in a data set.

The easiest way to create the initial recipe is:

``` r

rec_obj <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec_obj
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2
```

The `var_info` argument corresponds to the variable to be considered in
the modeling process and `tax_info` indicates the taxonomic level that
will be used for the analyses.

## Preprocessing Steps

From here, preprocessing steps for some step X can be added sequentially
in one of two ways:

``` r
rec_obj <- step_{X}(rec_obj, arguments)   
## or
rec_obj <- rec_obj |> step_{X}(arguments)
```

Note that all `step_ancom` and the other functions will always return
updated recipes.

We have two types of steps, those in charge of processing (prepro) and
those destined to define the methods of differential analysis (da).

The prepro steps are used to modify the data loaded into the recipe
which will then be used for the da steps. The `dar` package include 3
main preprocessing functionalities.

- `step_subset_taxa`: Is used for subsetting the columns and values
  within the tax_table.

- `step_filter_taxa`: Is used for filtering OTUs from recipe objects.

- `step_rarefaction`: Is used to resample an OTU table such that all
  samples have the same library size.

Additionally, the `dar` package provides convenient wrappers for the
`step_filter_taxa` function, designed to filter Operational Taxonomic
Units (OTUs) based on specific criteria: prevalence, variance,
abundance, and rarity.

- `step_filter_by_prevalence`: Filters OTUs according to the number of
  samples in which the OTU appears.
- `step_filter_by_variance`: Filters OTUs based on the variance of the
  OTU’s presence across samples.
- `step_filter_by_abundance`: Filters OTUs according to the OTU’s
  abundance across samples.
- `step_filter_by_rarity`: Filters OTUs based on the rarity of the OTU
  across samples.

For our data, we can add an operation to preprocessing the data stored
in the initial recpie. First, we will use `step_subset_taxa` to retain
only Bacteria and Archaea OTUs from the Kingdom taxonomic level. We will
then filter out OTUs where at least 3% of the samples have counts
greater than 0.

``` r

rec_obj <- rec_obj |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_by_prevalence(0.03)
  
rec_obj
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preprocessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Puff_pastry 
#>      ◉ step_filter_by_prevalence() id = filter_by_prevalence__Paris–Brest 
#> 
#> DA steps:
```

## Differential Analysis

Now that we have defined the preprocessing of the input data for all the
da methods that will be used, we need to define them. For this
introduction we will use **maaslin2** method with default parameters
(those defined by the authors of each method).

``` r

rec_obj <- rec_obj |>
  step_deseq() |>
  step_maaslin()

rec_obj
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preprocessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Puff_pastry 
#>      ◉ step_filter_by_prevalence() id = filter_by_prevalence__Paris–Brest 
#> 
#> DA steps:
#> 
#>      ◉ step_deseq() id = deseq__Spina_santa 
#>      ◉ step_maaslin() id = maaslin__Dutch_Baby_Pancake
```

The `dar` package includes more da steps than those defined above. Below
is the full list:

``` r

grep(
  "_new|_to_expr|filter|subset|rarefaction",
  grep("^step_", ls("package:dar"), value = TRUE),
  value = TRUE,
  invert = TRUE
)
#> [1] "step_aldex"   "step_ancom"   "step_corncob" "step_deseq"   "step_lefse"  
#> [6] "step_linda"   "step_maaslin" "step_wilcox"
```

## Prep

To ensure the reproducibility and consistency of the generated results,
all the steps defined in the recipe are executed at the same time using
the `prep` function.

``` r

da_results <- prep(rec_obj, parallel = TRUE)
da_results
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Results:
#> 
#>      ✔ deseq__Spina_santa diff_taxa = 155 
#>      ✔ maaslin__Dutch_Baby_Pancake diff_taxa = 41 
#> 
#>      ℹ 34 taxon-contrast effects are present in all tested methods
```

Note that the resulting object print shows information about the amount
of differentially abundant OTUs in each of the methods, as well as the
number of OTUs that have been detected by all methods (consensus).

## Bake and cool

Now that we have the results we need to extract them, however for this
we first need to define a consensus strategy with the `bake`. For this
example we are only interested in those OTUs detected as differentially
abundant in the three methods used.

``` r

## Number of used methods
count <- steps_ids(da_results, type = "da") |> length()

## Define the bake 
da_results <- bake(da_results, count_cutoff = count)
```

Finally we can extract the table with the results using the `cool`
function.

``` r

cool(da_results)
#> # A tibble: 34 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_102 Prevo… RiskGroup2… RiskGroup… main          Risk… up                2
#>  2 Otu_115 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  3 Otu_115 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  4 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  5 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  6 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  7 Otu_129 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  8 Otu_130 Parab… RiskGroup2… RiskGroup… main          Risk… up                2
#>  9 Otu_255 Rumin… RiskGroup2… RiskGroup… main          Risk… up                2
#> 10 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                2
#> # ℹ 24 more rows
#> # ℹ 1 more variable: methods <chr>
```

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
#>  package                  * version    date (UTC) lib source
#>  abind                      1.4-8      2024-09-12 [1] RSPM (R 4.6.0)
#>  ade4                       1.7-24     2026-03-21 [1] RSPM (R 4.6.0)
#>  ape                        5.8-1      2024-12-16 [1] RSPM (R 4.6.0)
#>  apeglm                     1.35.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  ashr                       2.2-63     2023-08-21 [1] RSPM (R 4.6.0)
#>  backports                  1.5.1      2026-04-03 [1] RSPM (R 4.6.0)
#>  bbmle                      1.0.25.1   2023-12-09 [1] RSPM (R 4.6.0)
#>  bdsmatrix                  1.3-7      2024-03-02 [1] RSPM (R 4.6.0)
#>  beachmat                   2.29.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  beeswarm                   0.4.0      2021-06-01 [1] RSPM (R 4.6.0)
#>  Biobase                    2.73.2     2026-07-29 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocGenerics               0.59.10    2026-07-07 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocNeighbors              2.7.2      2026-05-10 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocParallel               1.47.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  BiocSingular               1.29.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  biomformat                 1.41.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  Biostrings                 2.81.6     2026-07-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  bluster                    1.23.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  bslib                      0.12.0     2026-08-04 [2] RSPM (R 4.6.0)
#>  cachem                     1.1.0      2024-05-16 [2] RSPM (R 4.6.0)
#>  checkmate                  2.3.4      2026-02-03 [1] RSPM (R 4.6.0)
#>  cli                        3.6.6      2026-04-09 [2] RSPM (R 4.6.0)
#>  cluster                    2.1.8.3    2026-07-30 [2] RSPM (R 4.6.0)
#>  coda                       0.19-4.1   2024-01-31 [1] RSPM (R 4.6.0)
#>  codetools                  0.2-20     2024-03-31 [3] CRAN (R 4.6.1)
#>  crayon                     1.5.3      2024-06-20 [2] RSPM (R 4.6.0)
#>  dar                      * 1.9.12     2026-08-04 [1] Bioconductor
#>  data.table                 1.18.4     2026-05-06 [1] RSPM (R 4.6.0)
#>  DBI                        1.3.0      2026-02-25 [1] RSPM (R 4.6.0)
#>  DECIPHER                   3.9.2      2026-07-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  decontam                   1.33.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  DelayedArray               0.39.3     2026-06-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  DelayedMatrixStats         1.35.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  desc                       1.4.3      2023-12-10 [2] RSPM (R 4.6.0)
#>  DESeq2                     1.53.2     2026-07-20 [1] Bioconductor 3.24 (R 4.6.1)
#>  devtools                   2.5.2      2026-04-30 [2] RSPM (R 4.6.0)
#>  digest                     0.6.39     2025-11-19 [2] RSPM (R 4.6.0)
#>  DirichletMultinomial       1.55.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  dplyr                      1.2.1      2026-04-03 [1] RSPM (R 4.6.0)
#>  ellipsis                   0.3.3      2026-04-04 [2] RSPM (R 4.6.0)
#>  emdbook                    1.3.14     2025-07-23 [1] RSPM (R 4.6.0)
#>  evaluate                   1.0.5      2025-08-27 [2] RSPM (R 4.6.0)
#>  farver                     2.1.2      2024-05-13 [1] RSPM (R 4.6.0)
#>  fastmap                    1.2.0      2024-05-15 [2] RSPM (R 4.6.0)
#>  foreach                    1.5.2      2022-02-02 [1] RSPM (R 4.6.0)
#>  fs                         2.1.0      2026-04-18 [2] RSPM (R 4.6.0)
#>  furrr                      0.4.0      2026-03-31 [1] RSPM (R 4.6.0)
#>  future                     1.75.0     2026-07-20 [1] RSPM (R 4.6.0)
#>  generics                   0.1.4      2025-05-09 [1] RSPM (R 4.6.0)
#>  GenomicRanges              1.65.1     2026-07-08 [1] Bioconductor 3.24 (R 4.6.1)
#>  ggbeeswarm                 0.7.3      2025-11-29 [1] RSPM (R 4.6.0)
#>  ggplot2                    4.0.3      2026-04-22 [1] RSPM (R 4.6.0)
#>  ggrepel                    0.9.8      2026-03-17 [1] RSPM (R 4.6.0)
#>  globals                    0.19.1     2026-03-13 [1] RSPM (R 4.6.0)
#>  glue                       1.8.1      2026-04-17 [2] RSPM (R 4.6.0)
#>  gridExtra                  2.3.1      2026-06-25 [1] RSPM (R 4.6.0)
#>  gtable                     0.3.6      2024-10-25 [1] RSPM (R 4.6.0)
#>  htmltools                  0.5.9      2025-12-04 [2] RSPM (R 4.6.0)
#>  htmlwidgets                1.6.4      2023-12-06 [2] RSPM (R 4.6.0)
#>  igraph                     2.3.3      2026-06-26 [1] RSPM (R 4.6.0)
#>  invgamma                   1.2        2025-07-02 [1] RSPM (R 4.6.0)
#>  IRanges                    2.47.2     2026-06-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  irlba                      2.3.7      2026-01-30 [1] RSPM (R 4.6.0)
#>  iterators                  1.0.14     2022-02-05 [1] RSPM (R 4.6.0)
#>  jquerylib                  0.1.4      2021-04-26 [2] RSPM (R 4.6.0)
#>  jsonlite                   2.0.0      2025-03-27 [2] RSPM (R 4.6.0)
#>  knitr                      1.51       2025-12-20 [2] RSPM (R 4.6.0)
#>  lattice                    0.22-9     2026-02-09 [3] CRAN (R 4.6.1)
#>  lazyeval                   0.2.3      2026-04-04 [1] RSPM (R 4.6.0)
#>  lifecycle                  1.0.5      2026-01-08 [2] RSPM (R 4.6.0)
#>  listenv                    1.0.0      2026-06-22 [1] RSPM (R 4.6.0)
#>  locfit                     1.5-9.12   2025-03-05 [1] RSPM (R 4.6.0)
#>  maaslin3                   1.5.3      2026-06-23 [1] Bioconductor 3.24 (R 4.6.1)
#>  magrittr                   2.0.5      2026-04-04 [2] RSPM (R 4.6.0)
#>  MASS                       7.3-66     2026-07-15 [2] RSPM (R 4.6.0)
#>  Matrix                     1.7-6      2026-07-25 [2] RSPM (R 4.6.0)
#>  MatrixGenerics             1.25.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  matrixStats                1.5.0      2025-01-07 [1] RSPM (R 4.6.0)
#>  memoise                    2.0.1      2021-11-26 [2] RSPM (R 4.6.0)
#>  mgcv                       1.9-4      2025-11-07 [3] CRAN (R 4.6.1)
#>  mia                        1.21.6     2026-07-31 [1] Bioconductor 3.24 (R 4.6.1)
#>  mixsqp                     0.3-54     2023-12-20 [1] RSPM (R 4.6.0)
#>  MultiAssayExperiment       1.39.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  multtest                   2.69.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  mvtnorm                    1.4-2      2026-07-12 [1] RSPM (R 4.6.0)
#>  nlme                       3.1-170    2026-07-15 [2] RSPM (R 4.6.0)
#>  numDeriv                   2016.8-1.1 2019-06-06 [1] RSPM (R 4.6.0)
#>  optparse                   1.8.2      2026-04-17 [1] RSPM (R 4.6.0)
#>  otel                       0.2.0      2025-08-29 [2] RSPM (R 4.6.0)
#>  parallelly                 1.48.0     2026-06-29 [1] RSPM (R 4.6.0)
#>  permute                    0.9-10     2026-02-06 [1] RSPM (R 4.6.0)
#>  phyloseq                   1.57.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  pillar                     1.11.1     2025-09-17 [2] RSPM (R 4.6.0)
#>  pkgbuild                   1.4.8      2025-05-26 [2] RSPM (R 4.6.0)
#>  pkgconfig                  2.0.3      2019-09-22 [2] RSPM (R 4.6.0)
#>  pkgdown                    2.2.1      2026-07-07 [2] RSPM (R 4.6.0)
#>  pkgload                    1.5.3      2026-06-15 [2] RSPM (R 4.6.0)
#>  plyr                       1.8.9      2023-10-02 [1] RSPM (R 4.6.0)
#>  purrr                      1.2.2      2026-04-10 [2] RSPM (R 4.6.0)
#>  R6                         2.6.1      2025-02-15 [2] RSPM (R 4.6.0)
#>  ragg                       1.5.2      2026-03-23 [2] RSPM (R 4.6.0)
#>  rappdirs                   0.3.4      2026-01-17 [2] RSPM (R 4.6.0)
#>  rbibutils                  2.4.1      2026-01-21 [1] RSPM (R 4.6.0)
#>  RColorBrewer               1.1-3      2022-04-03 [1] RSPM (R 4.6.0)
#>  Rcpp                       1.1.2      2026-07-05 [2] RSPM (R 4.6.0)
#>  Rdpack                     2.6.6      2026-02-08 [1] RSPM (R 4.6.0)
#>  reformulas                 0.4.4      2026-02-02 [1] RSPM (R 4.6.0)
#>  reshape2                   1.4.5      2025-11-12 [1] RSPM (R 4.6.0)
#>  rlang                      1.3.0      2026-07-05 [2] RSPM (R 4.6.0)
#>  rmarkdown                  2.31       2026-03-26 [2] RSPM (R 4.6.0)
#>  rsvd                       1.0.5      2021-04-16 [1] RSPM (R 4.6.0)
#>  S4Arrays                   1.13.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  S4Vectors                  0.51.6     2026-07-29 [1] Bioconductor 3.24 (R 4.6.1)
#>  S7                         0.2.2      2026-04-22 [1] CRAN (R 4.6.1)
#>  sass                       0.4.10     2025-04-11 [2] RSPM (R 4.6.0)
#>  ScaledMatrix               1.21.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  scales                     1.4.0      2025-04-24 [1] RSPM (R 4.6.0)
#>  scater                     1.41.2     2026-06-30 [1] Bioconductor 3.24 (R 4.6.1)
#>  scuttle                    1.23.1     2026-05-18 [1] Bioconductor 3.24 (R 4.6.1)
#>  Seqinfo                    1.3.0      2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  sessioninfo                1.2.4      2026-06-04 [2] RSPM (R 4.6.0)
#>  SingleCellExperiment       1.35.2     2026-07-16 [1] Bioconductor 3.24 (R 4.6.1)
#>  SparseArray                1.13.2     2026-05-01 [1] Bioconductor 3.24 (R 4.6.1)
#>  sparseMatrixStats          1.25.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  SQUAREM                    2026.1     2026-03-12 [1] RSPM (R 4.6.0)
#>  stringi                    1.8.7      2025-03-27 [2] RSPM (R 4.6.0)
#>  stringr                    1.6.0      2025-11-04 [1] RSPM (R 4.6.0)
#>  SummarizedExperiment       1.43.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  survival                   3.8-9      2026-07-08 [2] RSPM (R 4.6.0)
#>  systemfonts                1.3.2      2026-03-05 [2] RSPM (R 4.6.0)
#>  textshaping                1.0.5      2026-03-06 [2] RSPM (R 4.6.0)
#>  tibble                     3.3.1      2026-01-11 [2] RSPM (R 4.6.0)
#>  tidyr                      1.3.2      2025-12-19 [1] RSPM (R 4.6.0)
#>  tidyselect                 1.2.1      2024-03-11 [1] RSPM (R 4.6.0)
#>  tidytree                   0.4.8      2026-07-02 [1] RSPM (R 4.6.0)
#>  treeio                     1.37.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  TreeSummarizedExperiment   2.21.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  truncnorm                  1.0-9      2023-03-20 [1] RSPM (R 4.6.0)
#>  usethis                    3.2.1      2025-09-06 [2] RSPM (R 4.6.0)
#>  utf8                       1.2.6      2025-06-08 [2] RSPM (R 4.6.0)
#>  vctrs                      0.7.3      2026-04-11 [2] RSPM (R 4.6.0)
#>  vegan                      2.7-5      2026-05-25 [1] RSPM (R 4.6.0)
#>  vipor                      0.4.7      2023-12-18 [1] RSPM (R 4.6.0)
#>  viridis                    0.6.5      2024-01-29 [1] RSPM (R 4.6.0)
#>  viridisLite                0.4.3      2026-02-04 [1] RSPM (R 4.6.0)
#>  withr                      3.0.3      2026-06-19 [2] RSPM (R 4.6.0)
#>  xfun                       0.60       2026-07-09 [2] RSPM (R 4.6.0)
#>  XVector                    0.53.0     2026-04-28 [1] Bioconductor 3.24 (R 4.6.1)
#>  yaml                       2.3.12     2025-12-10 [2] RSPM (R 4.6.0)
#>  yulab.utils                0.2.4      2026-02-02 [1] RSPM (R 4.6.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
