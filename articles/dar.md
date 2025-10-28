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

rec_obj <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species") 
rec_obj
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species
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
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Puff_pastry 
#>      ◉ step_filter_by_prevalence() id = filter_by_prevalence__Paris–Brest 
#> 
#> DA steps:
```

## Differential Analysis

Now that we have defined the preprocessing of the input data for all the
da methods that will be used, we need to define them. For this
introduction we will use **metagenomeSeq** and **maaslin2** methods with
default parameters (those defined by the authors of each method).

``` r

rec_obj <- rec_obj |>
  step_deseq() |>
  step_metagenomeseq(rm_zeros = 0.01) |>
  step_maaslin()

rec_obj
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Puff_pastry 
#>      ◉ step_filter_by_prevalence() id = filter_by_prevalence__Paris–Brest 
#> 
#> DA steps:
#> 
#>      ◉ step_deseq() id = deseq__Spina_santa 
#>      ◉ step_metagenomeseq() id = metagenomeseq__Dutch_Baby_Pancake 
#>      ◉ step_maaslin() id = maaslin__Knieküchle
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
#> [1] "step_aldex"         "step_ancom"         "step_corncob"      
#> [4] "step_deseq"         "step_lefse"         "step_maaslin"      
#> [7] "step_metagenomeseq" "step_wilcox"
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
#>      ✔ deseq__Spina_santa diff_taxa = 166 
#>      ✔ metagenomeseq__Dutch_Baby_Pancake diff_taxa = 236 
#>      ✔ maaslin__Knieküchle diff_taxa = 146 
#> 
#>      ℹ 65 taxa are present in all tested methods
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
#> # A tibble: 65 × 2
#>    taxa_id taxa                        
#>    <chr>   <chr>                       
#>  1 Otu_15  Bifidobacterium_catenulatum 
#>  2 Otu_34  Olsenella_scatoligenes      
#>  3 Otu_35  Collinsella_aerofaciens     
#>  4 Otu_37  Collinsella_stercoris       
#>  5 Otu_38  Enorma_massiliensis         
#>  6 Otu_45  Slackia_isoflavoniconvertens
#>  7 Otu_47  Bacteroides_cellulosilyticus
#>  8 Otu_48  Bacteroides_clarus          
#>  9 Otu_63  Bacteroides_plebeius        
#> 10 Otu_69  Bacteroides_sp_CAG_530      
#> # ℹ 55 more rows
```

## Session info

``` r

devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.1 (2025-06-13)
#>  os       Ubuntu 24.04.3 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       UTC
#>  date     2025-10-28
#>  pandoc   3.8.1 @ /usr/bin/ (via rmarkdown)
#>  quarto   1.7.32 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package          * version    date (UTC) lib source
#>  ade4               1.7-23     2025-02-14 [1] RSPM (R 4.5.0)
#>  ape                5.8-1      2024-12-16 [1] RSPM (R 4.5.0)
#>  assertthat         0.2.1      2019-03-21 [1] RSPM (R 4.5.0)
#>  Biobase            2.68.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  BiocGenerics       0.54.1     2025-10-12 [1] Bioconductor 3.21 (R 4.5.1)
#>  biomformat         1.36.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  Biostrings         2.76.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  brio               1.1.5      2024-04-24 [2] RSPM (R 4.5.0)
#>  bslib              0.9.0      2025-01-30 [2] RSPM (R 4.5.0)
#>  ca                 0.71.1     2020-01-24 [1] RSPM (R 4.5.0)
#>  cachem             1.1.0      2024-05-16 [2] RSPM (R 4.5.0)
#>  cli                3.6.5      2025-04-23 [2] RSPM (R 4.5.0)
#>  cluster            2.1.8.1    2025-03-12 [3] CRAN (R 4.5.1)
#>  codetools          0.2-20     2024-03-31 [3] CRAN (R 4.5.1)
#>  crayon             1.5.3      2024-06-20 [2] RSPM (R 4.5.0)
#>  dar              * 1.5.3      2025-10-28 [1] Bioconductor
#>  data.table         1.17.8     2025-07-10 [1] RSPM (R 4.5.0)
#>  dendextend         1.19.1     2025-07-15 [1] RSPM (R 4.5.0)
#>  desc               1.4.3      2023-12-10 [2] RSPM (R 4.5.0)
#>  devtools           2.4.6      2025-10-03 [2] RSPM (R 4.5.0)
#>  digest             0.6.37     2024-08-19 [2] RSPM (R 4.5.0)
#>  dplyr              1.1.4      2023-11-17 [1] RSPM (R 4.5.0)
#>  ellipsis           0.3.2      2021-04-29 [2] RSPM (R 4.5.0)
#>  evaluate           1.0.5      2025-08-27 [2] RSPM (R 4.5.0)
#>  farver             2.1.2      2024-05-13 [1] RSPM (R 4.5.0)
#>  fastmap            1.2.0      2024-05-15 [2] RSPM (R 4.5.0)
#>  foreach            1.5.2      2022-02-02 [1] RSPM (R 4.5.0)
#>  fs                 1.6.6      2025-04-12 [2] RSPM (R 4.5.0)
#>  furrr              0.3.1      2022-08-15 [1] RSPM (R 4.5.0)
#>  future             1.67.0     2025-07-29 [1] RSPM (R 4.5.0)
#>  generics           0.1.4      2025-05-09 [1] RSPM (R 4.5.0)
#>  GenomeInfoDb       1.44.3     2025-09-21 [1] Bioconductor 3.21 (R 4.5.1)
#>  GenomeInfoDbData   1.2.14     2025-10-28 [1] Bioconductor
#>  ggplot2            4.0.0      2025-09-11 [1] RSPM (R 4.5.0)
#>  globals            0.18.0     2025-05-08 [1] RSPM (R 4.5.0)
#>  glue               1.8.0      2024-09-30 [2] RSPM (R 4.5.0)
#>  gridExtra          2.3        2017-09-09 [1] RSPM (R 4.5.0)
#>  gtable             0.3.6      2024-10-25 [1] RSPM (R 4.5.0)
#>  heatmaply          1.6.0      2025-07-12 [1] RSPM (R 4.5.0)
#>  htmltools          0.5.8.1    2024-04-04 [2] RSPM (R 4.5.0)
#>  htmlwidgets        1.6.4      2023-12-06 [2] RSPM (R 4.5.0)
#>  httr               1.4.7      2023-08-15 [1] RSPM (R 4.5.0)
#>  igraph             2.2.0      2025-10-13 [1] RSPM (R 4.5.0)
#>  IRanges            2.42.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  iterators          1.0.14     2022-02-05 [1] RSPM (R 4.5.0)
#>  jquerylib          0.1.4      2021-04-26 [2] RSPM (R 4.5.0)
#>  jsonlite           2.0.0      2025-03-27 [2] RSPM (R 4.5.0)
#>  knitr              1.50       2025-03-16 [2] RSPM (R 4.5.0)
#>  lattice            0.22-7     2025-04-02 [3] CRAN (R 4.5.1)
#>  lazyeval           0.2.2      2019-03-15 [1] RSPM (R 4.5.0)
#>  lifecycle          1.0.4      2023-11-07 [2] RSPM (R 4.5.0)
#>  listenv            0.9.1      2024-01-29 [1] RSPM (R 4.5.0)
#>  magrittr           2.0.4      2025-09-12 [2] RSPM (R 4.5.0)
#>  MASS               7.3-65     2025-02-28 [3] CRAN (R 4.5.1)
#>  Matrix             1.7-4      2025-08-28 [2] RSPM (R 4.5.0)
#>  memoise            2.0.1      2021-11-26 [2] RSPM (R 4.5.0)
#>  mgcv               1.9-3      2025-04-04 [3] CRAN (R 4.5.1)
#>  multtest           2.64.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  nlme               3.1-168    2025-03-31 [3] CRAN (R 4.5.1)
#>  parallelly         1.45.1     2025-07-24 [1] RSPM (R 4.5.0)
#>  permute            0.9-8      2025-06-25 [1] RSPM (R 4.5.0)
#>  phyloseq           1.52.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  pillar             1.11.1     2025-09-17 [2] RSPM (R 4.5.0)
#>  pkgbuild           1.4.8      2025-05-26 [2] RSPM (R 4.5.0)
#>  pkgconfig          2.0.3      2019-09-22 [2] RSPM (R 4.5.0)
#>  pkgdown            2.1.3.9000 2025-10-28 [1] Github (r-lib/pkgdown@f1778e3)
#>  pkgload            1.4.1      2025-09-23 [2] RSPM (R 4.5.0)
#>  plotly             4.11.0     2025-06-19 [1] RSPM (R 4.5.0)
#>  plyr               1.8.9      2023-10-02 [1] RSPM (R 4.5.0)
#>  purrr              1.1.0      2025-07-10 [2] RSPM (R 4.5.0)
#>  R6                 2.6.1      2025-02-15 [2] RSPM (R 4.5.0)
#>  ragg               1.5.0      2025-09-02 [2] RSPM (R 4.5.0)
#>  RColorBrewer       1.1-3      2022-04-03 [1] RSPM (R 4.5.0)
#>  Rcpp               1.1.0      2025-07-02 [2] RSPM (R 4.5.0)
#>  registry           0.5-1      2019-03-05 [1] RSPM (R 4.5.0)
#>  remotes            2.5.0      2024-03-17 [1] RSPM (R 4.5.0)
#>  reshape2           1.4.4      2020-04-09 [1] RSPM (R 4.5.0)
#>  rhdf5              2.52.1     2025-06-08 [1] Bioconductor 3.21 (R 4.5.1)
#>  rhdf5filters       1.20.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  Rhdf5lib           1.30.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  rlang              1.1.6      2025-04-11 [2] RSPM (R 4.5.0)
#>  rmarkdown          2.30       2025-09-28 [2] RSPM (R 4.5.0)
#>  S4Vectors          0.46.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  S7                 0.2.0      2024-11-07 [1] RSPM (R 4.5.0)
#>  sass               0.4.10     2025-04-11 [2] RSPM (R 4.5.0)
#>  scales             1.4.0      2025-04-24 [1] RSPM (R 4.5.0)
#>  seriation          1.5.8      2025-08-20 [1] RSPM (R 4.5.0)
#>  sessioninfo        1.2.3      2025-02-05 [2] RSPM (R 4.5.0)
#>  stringi            1.8.7      2025-03-27 [2] RSPM (R 4.5.0)
#>  stringr            1.5.2      2025-09-08 [2] RSPM (R 4.5.0)
#>  survival           3.8-3      2024-12-17 [3] CRAN (R 4.5.1)
#>  systemfonts        1.3.1      2025-10-01 [2] RSPM (R 4.5.0)
#>  testthat           3.2.3      2025-01-13 [2] RSPM (R 4.5.0)
#>  textshaping        1.0.4      2025-10-10 [2] RSPM (R 4.5.0)
#>  tibble             3.3.0      2025-06-08 [2] RSPM (R 4.5.0)
#>  tidyr              1.3.1      2024-01-24 [1] RSPM (R 4.5.0)
#>  tidyselect         1.2.1      2024-03-11 [1] RSPM (R 4.5.0)
#>  TSP                1.2-5      2025-05-27 [1] RSPM (R 4.5.0)
#>  UCSC.utils         1.4.0      2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  UpSetR             1.4.0      2019-05-22 [1] RSPM (R 4.5.0)
#>  usethis            3.2.1      2025-09-06 [2] RSPM (R 4.5.0)
#>  utf8               1.2.6      2025-06-08 [2] RSPM (R 4.5.0)
#>  vctrs              0.6.5      2023-12-01 [2] RSPM (R 4.5.0)
#>  vegan              2.7-2      2025-10-08 [1] RSPM (R 4.5.0)
#>  viridis            0.6.5      2024-01-29 [1] RSPM (R 4.5.0)
#>  viridisLite        0.4.2      2023-05-02 [1] RSPM (R 4.5.0)
#>  webshot            0.5.5      2023-06-26 [1] RSPM (R 4.5.0)
#>  withr              3.0.2      2024-10-28 [2] RSPM (R 4.5.0)
#>  xfun               0.53       2025-08-19 [2] RSPM (R 4.5.0)
#>  XVector            0.48.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
#>  yaml               2.3.10     2024-07-26 [2] RSPM (R 4.5.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
