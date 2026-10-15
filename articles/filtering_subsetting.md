# Quality control and preprocessing

Preprocessing steps describe changes to the data stored in a recipe.
Adding a step does not execute it:
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
applies preprocessing in the configured order, before running
differential-abundance methods. A recipe without a statistical model is
valid for preprocessing alone.

The examples below use the included data and run only preprocessing, not
statistical models. The thresholds illustrate the API; they are not
universal recommendations for microbiome studies.

## Inspect overall and grouped QC

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dar`](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar)`)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"metaHIV_phy"``, package ``=`` ``"dar"``)`\
\
`rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`\
[`recipe_qc`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)`(``rec``)`\
`#> ``# A tibble: 1 × 13`\
`#>   group_by group     n n_zero pct_zero pct_all_zero pct_singletons`\
`#>   ``<chr>``    ``<chr>`` ``<dbl>``  ``<dbl>``    ``<dbl>``        ``<dbl>``          ``<dbl>`\
`#> ``1`` ``NA``       ``NA``    ``70``356  ``57``632     81.9            0           20.6`\
`#> ``# ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,`\
`#> ``#   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>`\
[`recipe_qc`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)`(``rec``, group_by ``=`` ``"RiskGroup2"``)`\
`#> ``# A tibble: 4 × 13`\
`#>   group_by   group     n n_zero pct_zero pct_all_zero pct_singletons`\
`#>   ``<chr>``      ``<chr>`` ``<dbl>``  ``<dbl>``    ``<dbl>``        ``<dbl>``          ``<dbl>`\
`#> ``1`` ``NA``         ``NA``    ``70``356  ``57``632     81.9          0             20.6`\
`#> ``2`` RiskGroup2 hts   ``18``491  ``15``108     81.7         24.2           22.8`\
`#> ``3`` RiskGroup2 msm   ``45``100  ``37``019     82.1         16.0           20.2`\
`#> ``4`` RiskGroup2 pwid   ``6``765   ``5``505     81.4         41.2           16.6`\
`#> ``# ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,`\
`#> ``#   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>`

[`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
summarizes the stored abundance matrix independently of any model. It
reports sample counts, sparsity, taxa detected in zero, one or two
samples, and library-size and count summaries. Singleton and doubleton
metrics refer to the number of samples with a detection, not to a
taxon’s total read count. The `n` column counts matrix entries, not
biological replicates.

With `group_by`, the overall row is followed by one row per observed
metadata group. The overall row has missing `group_by` and `group`
identifiers; it is not a group called `"all"`. Grouping is explicit and
is not inferred from the statistical target.

Use these summaries to inspect depth or sparsity differences between
groups. They do not establish the absence of confounding or that a
particular filter will improve statistical accuracy.

## Define and apply a shared preprocessing sequence

Here we first retain Bacteria and Archaea, then keep taxa observed in
more than 3% of samples. These operations act on the stored features;
selecting a rank later with
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
determines the taxonomic analysis view used by the DA methods.

\
`filtered_rec`` ``<-`` ``rec`` ``|>`\
`  `[`step_subset_taxa`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_subset_taxa.md)`(`\
`    tax_level ``=`` ``"Kingdom"``, taxa ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bacteria"``, ``"Archaea"``)`\
`  ``)`` ``|>`\
`  `[`step_filter_by_prevalence`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md)`(``threshold ``=`` ``0.03``)`\
\
`filtered`` ``<-`` `[`prep`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)`(``filtered_rec``, parallel ``=`` ``FALSE``)`

[`step_subset_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_subset_taxa.md)
matches annotations at the specified rank. The prevalence filter then
uses all samples in the recipe, not one threshold calculated separately
within each QC group. Preprocessing preserves the original input object;
the prepared recipe contains the filtered data.

Compare the stored data before and after preprocessing:

\
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  stage ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"before"``, ``"after"``)``,`\
`  n_taxa ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``dar``::`[`otu_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)`(``rec``)``)``, `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``dar``::`[`otu_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)`(``filtered``)``)``)`\
`)`\
`#>    stage n_taxa`\
`#> 1 before    451`\
`#> 2  after    278`\
[`recipe_qc`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)`(``filtered``)`\
`#> ``# A tibble: 1 × 13`\
`#>   group_by group     n n_zero pct_zero pct_all_zero pct_singletons`\
`#>   ``<chr>``    ``<chr>`` ``<dbl>``  ``<dbl>``    ``<dbl>``        ``<dbl>``          ``<dbl>`\
`#> ``1`` ``NA``       ``NA``    ``43``368  ``30``986     71.4            0              0`\
`#> ``# ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,`\
`#> ``#   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>`\
[`recipe_qc`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)`(``filtered``, group_by ``=`` ``"RiskGroup2"``)`\
`#> ``# A tibble: 4 × 13`\
`#>   group_by   group     n n_zero pct_zero pct_all_zero pct_singletons`\
`#>   ``<chr>``      ``<chr>`` ``<dbl>``  ``<dbl>``    ``<dbl>``        ``<dbl>``          ``<dbl>`\
`#> ``1`` ``NA``         ``NA``    ``43``368  ``30``986     71.4         0              0   `\
`#> ``2`` RiskGroup2 hts   ``11``398   ``8``117     71.2         7.19          12.6 `\
`#> ``3`` RiskGroup2 msm   ``27``800  ``19``912     71.6         2.52           6.12`\
`#> ``4`` RiskGroup2 pwid   ``4``170   ``2``957     70.9        18.7           15.8 `\
`#> ``# ℹ 6 more variables: pct_doubletons <dbl>, n_samples <dbl>,`\
`#> ``#   lib_size_min <dbl>, lib_size_max <dbl>, count_mean <dbl>, count_max <dbl>`

Order matters. An abundance threshold expressed relative to all reads
uses the total remaining at that point in the sequence. Resampling
before a prevalence filter can change which taxa pass it. Specify the
intended order explicitly, and remember that engines may also apply
their own internal filters. A shared input does not make every engine’s
tested-feature set or assumptions identical.

## Choose a filtering predicate

The convenience wrappers apply across-sample criteria to each stored
taxon:

| Step | Criterion for retaining a taxon |
|:---|:---|
| `step_filter_by_prevalence(threshold)` | Proportion of samples with nonzero abundance is strictly greater than the threshold. |
| `step_filter_by_abundance(threshold)` | Its total abundance is strictly greater than the threshold times the total abundance of the current table. |
| `step_filter_by_variance(threshold)` | Its across-sample variance is strictly greater than the threshold, on the stored abundance scale. |
| `step_filter_by_rarity(threshold)` | Proportion of samples with nonzero abundance is strictly below the threshold. |

[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md)
selects rare taxa; it does not remove rare taxa. Use
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md)
for that purpose. The variance threshold is not a proportion and depends
on the abundance units.

For a custom rule,
[`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)
accepts a function with exactly one argument named `x` and no default.
This example includes taxa exactly on the 3% boundary, unlike the strict
prevalence wrapper above.

\
`custom_rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`` ``|>`\
`  `[`step_filter_taxa`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)`(``.f ``=`` ``function``(``x``)`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``x`` ``>`` ``0``)`` ``>=`` ``0.03`` ``*`` `[`length`](https://rdrr.io/r/base/length.html)`(``x``)``)`

The same predicate may be supplied as a character string:

\
`custom_string_rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`` ``|>`\
`  `[`step_filter_taxa`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)`(``.f ``=`` ``"function(x) sum(x > 0) >= 0.03 * length(x)"``)`

Both forms return a recipe, not filtered data or a logical vector. The
predicate must return one logical value when applied to each abundance
vector during
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md);
the wrapper has no `prune` argument.

To keep recipes safely serializable, the predicate body must be a single
expression. Braces, assignments, external variables, namespace-qualified
calls, user-defined helpers and lists of functions are not supported,
even in local R use. Write thresholds as literal numbers. Accepted
expressions include arithmetic, comparisons, logical operators and
approved summaries such as `sum`, `mean`, `median`, `var` and `sd`; see
[`?step_filter_taxa`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)
for the complete syntax and supported `na.rm` arguments. The
restrictions are checked when the step is added, not only on export.

## Rarefaction is optional

The workflow above does not rarefy.
[`step_rarefaction()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rarefaction.md)
is available when an explicitly chosen analysis calls for shared
equal-depth resampling:

\
`rarefied_rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`` ``|>`\
`  `[`step_rarefaction`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rarefaction.md)`(``)`\
`rarefied`` ``<-`` `[`prep`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)`(``rarefied_rec``, parallel ``=`` ``FALSE``)`

This is not a general recommendation to normalize DA inputs by
rarefaction. Subsampling discards observations and changes the data
supplied to downstream methods. Select preprocessing in the context of
the scientific question and the methods’ assumptions; equal library
sizes do not guarantee valid inference.
[`rarefaction_help()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/rarefaction_help.md)
describes the cautions in more detail.

Some DA steps also accept `rarefy = TRUE` to resample only for that
step, rather than altering the shared preprocessing data. For those
configured DA steps, the informational note about rarefaction’s
controversy appears during
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md),
not when the step is added. A shared
[`step_rarefaction()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rarefaction.md)
and a method-specific `rarefy` option are different configurations; do
not combine them unless repeat resampling is intentional.

For selecting a central statistical model after defining preprocessing,
see [Statistical models and
contrasts](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/statistical_models.md).

## Session info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.1 (2026-06-24)`\
`#> Platform: x86_64-pc-linux-gnu`\
`#> Running under: Ubuntu 24.04.4 LTS`\
`#> `\
`#> Matrix products: default`\
`#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 `\
`#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0`\
`#> `\
`#> locale:`\
`#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              `\
`#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    `\
`#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   `\
`#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 `\
`#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            `\
`#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       `\
`#> `\
`#> time zone: UTC`\
`#> tzcode source: system (glibc)`\
`#> `\
`#> attached base packages:`\
`#> [1] stats     graphics  grDevices utils     datasets  methods   base     `\
`#> `\
`#> other attached packages:`\
`#> [1] dar_1.9.12`\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>   [1] DBI_1.3.0                       gridExtra_2.3.1                `\
`#>   [3] permute_0.9-10                  rlang_1.3.0                    `\
`#>   [5] magrittr_2.0.5                  ade4_1.7-24                    `\
`#>   [7] scater_1.41.2                   otel_0.2.0                     `\
`#>   [9] matrixStats_1.5.0               compiler_4.6.1                 `\
`#>  [11] DelayedMatrixStats_1.35.0       mgcv_1.9-4                     `\
`#>  [13] systemfonts_1.3.2               vctrs_0.7.3                    `\
`#>  [15] reshape2_1.4.5                  stringr_1.6.0                  `\
`#>  [17] pkgconfig_2.0.3                 crayon_1.5.3                   `\
`#>  [19] fastmap_1.2.0                   backports_1.5.1                `\
`#>  [21] XVector_0.53.0                  scuttle_1.23.2                 `\
`#>  [23] utf8_1.2.6                      rmarkdown_2.32                 `\
`#>  [25] ggbeeswarm_0.7.3                DirichletMultinomial_1.55.1    `\
`#>  [27] ragg_1.5.2                      purrr_1.2.2                    `\
`#>  [29] bluster_1.23.1                  MultiAssayExperiment_1.39.1    `\
`#>  [31] xfun_0.61                       cachem_1.1.0                   `\
`#>  [33] beachmat_2.29.3                 jsonlite_2.0.0                 `\
`#>  [35] biomformat_1.41.0               DelayedArray_0.39.8            `\
`#>  [37] BiocParallel_1.47.0             irlba_2.4.1                    `\
`#>  [39] parallel_4.6.1                  cluster_2.1.8.3                `\
`#>  [41] R6_2.6.1                        bslib_0.12.0                   `\
`#>  [43] stringi_1.8.9                   RColorBrewer_1.1-3             `\
`#>  [45] GenomicRanges_1.65.4            jquerylib_0.1.4                `\
`#>  [47] Rcpp_1.1.2                      Seqinfo_1.3.2                  `\
`#>  [49] SummarizedExperiment_1.43.0     iterators_1.0.14               `\
`#>  [51] knitr_1.52                      DECIPHER_3.9.3                 `\
`#>  [53] IRanges_2.47.5                  mia_1.21.11                    `\
`#>  [55] Matrix_1.7-6                    splines_4.6.1                  `\
`#>  [57] igraph_2.3.4                    tidyselect_1.2.1               `\
`#>  [59] viridis_0.6.5                   abind_1.4-8                    `\
`#>  [61] yaml_2.3.12                     vegan_2.7-6                    `\
`#>  [63] TreeSummarizedExperiment_2.21.0 codetools_0.2-20               `\
`#>  [65] lattice_0.23-1                  tibble_3.3.1                   `\
`#>  [67] plyr_1.8.9                      withr_3.0.3                    `\
`#>  [69] Biobase_2.73.2                  treeio_1.37.1                  `\
`#>  [71] S7_0.2.2                        evaluate_1.0.5                 `\
`#>  [73] desc_1.4.3                      survival_3.8-12                `\
`#>  [75] Biostrings_2.81.9               pillar_1.11.1                  `\
`#>  [77] phyloseq_1.57.0                 MatrixGenerics_1.25.0          `\
`#>  [79] checkmate_2.3.4                 foreach_1.5.2                  `\
`#>  [81] stats4_4.6.1                    generics_0.1.4                 `\
`#>  [83] S4Vectors_0.51.10               ggplot2_4.0.3                  `\
`#>  [85] sparseMatrixStats_1.25.0        scales_1.4.0                   `\
`#>  [87] tidytree_0.4.8                  glue_1.8.1                     `\
`#>  [89] lazyeval_0.2.3                  tools_4.6.1                    `\
`#>  [91] BiocNeighbors_2.7.3             data.table_1.18.6.1            `\
`#>  [93] ScaledMatrix_1.21.0             fs_2.1.0                       `\
`#>  [95] grid_4.6.1                      tidyr_1.3.2                    `\
`#>  [97] ape_5.8-1                       SingleCellExperiment_1.35.2    `\
`#>  [99] nlme_3.1-171                    beeswarm_0.4.0                 `\
`#> [101] BiocSingular_1.29.2             vipor_0.4.7                    `\
`#> [103] cli_3.6.6                       rsvd_1.0.5                     `\
`#> [105] rappdirs_0.3.4                  textshaping_1.0.5              `\
`#> [107] viridisLite_0.4.3               S4Arrays_1.13.2                `\
`#> [109] dplyr_1.2.1                     gtable_0.3.6                   `\
`#> [111] yulab.utils_0.2.5               sass_0.4.10                    `\
`#> [113] digest_0.6.39                   BiocGenerics_0.59.12           `\
`#> [115] ggrepel_0.9.8                   SparseArray_1.13.4             `\
`#> [117] decontam_1.33.0                 htmlwidgets_1.6.4              `\
`#> [119] farver_2.1.2                    htmltools_0.5.9                `\
`#> [121] pkgdown_2.2.1                   multtest_2.69.0                `\
`#> [123] lifecycle_1.0.5                 MASS_7.3-66`
