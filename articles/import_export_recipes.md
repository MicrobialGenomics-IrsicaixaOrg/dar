# Saving and sharing analysis recipes

[`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
saves an analysis configuration;
[`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
restores it onto data supplied by the user. The versioned,
non-executable `dar-recipe` JSON format stores the central model,
ordered preprocessing and differential abundance steps, their
identifiers and native options. Optional consensus configurations are
stored separately when exporting a prepared recipe.

The file does not contain the input data, fitted results, active
random-number state or software environment. Sharing the JSON alone
therefore does not guarantee identical results on another computer.

## Save a configuration

This example configures an analysis without fitting a model. The
temporary path keeps the example self-contained; use a persistent path
when saving a recipe for later work.

\
[`data`](https://rdrr.io/r/utils/data.html)`(``metaHIV_phy``)`\
\
`configured`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`` ``|>`\
`  `[`add_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)`(`\
`    ``~`` ``RiskGroup2``,`\
`    targets ``=`` ``"RiskGroup2"``,`\
`    tax_level ``=`` ``"Species"`\
`  ``)`` ``|>`\
`  `[`step_subset_taxa`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_subset_taxa.md)`(``tax_level ``=`` ``"Kingdom"``, taxa ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Bacteria"``, ``"Archaea"``)``)`` ``|>`\
`  `[`step_filter_by_prevalence`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md)`(``threshold ``=`` ``0.03``)`` ``|>`\
`  `[`step_maaslin`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md)`(``)`\
\
`recipe_file`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``fileext ``=`` ``".json"``)`\
[`export_steps`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)`(``configured``, ``recipe_file``)`

The model record includes references, contrasts, the missing-value
policy and any time settings. See [Statistical models and
contrasts](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/statistical_models.md)
for configuring those fields.

## Restore onto data

Start with a new recipe containing the intended data. The JSON supplies
the model and steps; importing this file does not run the analysis.

\
`restored`` ``<-`` `[`import_steps`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)`(`[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)``, ``recipe_file``)`\
[`get_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)`(``restored``)`\
`#> $formula`\
`#> ~RiskGroup2`\
`#> <environment: base>`\
`#> `\
`#> $targets`\
`#> [1] "RiskGroup2"`\
`#> `\
`#> $tax_level`\
`#> [1] "Species"`\
`#> `\
`#> $time`\
`#> NULL`\
`#> `\
`#> $reference`\
`#> NULL`\
`#> `\
`#> $contrasts`\
`#> [1] "simple"`\
`#> `\
`#> $na_action`\
`#> [1] "error"`\
`#> `\
`#> $max_time_levels`\
`#> [1] 12`\
[`steps_ids`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/steps_ids.md)`(``restored``)`\
`#> [1] "subset_taxa__Ensaïmada"              `\
`#> [2] "filter_by_prevalence__Puff_pastry"   `\
`#> [3] "maaslin__Fazuelos_Fijuelas_or_Deblas"`

Using other data reuses the configuration, not the original analysis.
Check that the required metadata, factor levels, reference categories
and taxonomic ranks are present and appropriate for the new study.

Run
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
explicitly when ready. Set the execution seed immediately before the
run; it is not restored from the JSON:

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``130``)`\
`results`` ``<-`` `[`prep`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)`(``restored``, parallel ``=`` ``FALSE``)`\
[`tidy_results`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_results.md)`(``results``)`

[`bake()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/bake.md)
is not required to execute a recipe or inspect its per-method results.
If an exported file does contain consensus configurations, importing it
into an unprepared recipe requires explicit `prepare = TRUE`; otherwise
import stops without fitting models. The `parallel` and `workers`
arguments apply only to that explicit preparation path. Importing into
an already prepared recipe requires its model and analysis steps to
match the file.

## What to archive alongside the recipe

For a reproducible handoff, keep the input data and preprocessing
provenance, the execution script with its seed and execution settings,
fitted results, and the R, dar and engine-package versions.
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) records
loaded versions; an environment lockfile or container can additionally
preserve the dependencies. Stochastic methods and execution settings can
affect results even when the recipe parameters match.

Recipe files are treated as untrusted data: constructors, fields,
formulas and filter predicates are validated against explicit allowlists
before the recipe is changed. R functions in `engine_args` are not
serializable. Legacy executable exports are rejected; recreate their
configuration using the current API and export it again rather than
executing an old recipe file.

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
`#>   [1] Rdpack_2.6.6                    DBI_1.3.0                      `\
`#>   [3] gridExtra_2.3.1                 permute_0.9-10                 `\
`#>   [5] rlang_1.3.0                     magrittr_2.0.5                 `\
`#>   [7] ade4_1.7-24                     scater_1.41.2                  `\
`#>   [9] otel_0.2.0                      matrixStats_1.5.0              `\
`#>  [11] compiler_4.6.1                  DelayedMatrixStats_1.35.0      `\
`#>  [13] mgcv_1.9-4                      systemfonts_1.3.2              `\
`#>  [15] vctrs_0.7.3                     reshape2_1.4.5                 `\
`#>  [17] stringr_1.6.0                   pkgconfig_2.0.3                `\
`#>  [19] crayon_1.5.3                    fastmap_1.2.0                  `\
`#>  [21] backports_1.5.1                 XVector_0.53.0                 `\
`#>  [23] scuttle_1.23.2                  rmarkdown_2.32                 `\
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
`#>  [67] plyr_1.8.9                      Biobase_2.73.2                 `\
`#>  [69] treeio_1.37.1                   S7_0.2.2                       `\
`#>  [71] evaluate_1.0.5                  desc_1.4.3                     `\
`#>  [73] survival_3.8-12                 Biostrings_2.81.9              `\
`#>  [75] pillar_1.11.1                   phyloseq_1.57.0                `\
`#>  [77] MatrixGenerics_1.25.0           checkmate_2.3.4                `\
`#>  [79] foreach_1.5.2                   stats4_4.6.1                   `\
`#>  [81] reformulas_0.4.4                generics_0.1.4                 `\
`#>  [83] S4Vectors_0.51.10               ggplot2_4.0.3                  `\
`#>  [85] sparseMatrixStats_1.25.0        scales_1.4.0                   `\
`#>  [87] tidytree_0.4.8                  glue_1.8.1                     `\
`#>  [89] lazyeval_0.2.3                  tools_4.6.1                    `\
`#>  [91] maaslin3_1.5.7                  BiocNeighbors_2.7.3            `\
`#>  [93] data.table_1.18.6.1             ScaledMatrix_1.21.0            `\
`#>  [95] fs_2.1.0                        grid_4.6.1                     `\
`#>  [97] optparse_1.8.2                  tidyr_1.3.2                    `\
`#>  [99] ape_5.8-1                       rbibutils_2.4.1                `\
`#> [101] SingleCellExperiment_1.35.2     nlme_3.1-171                   `\
`#> [103] beeswarm_0.4.0                  BiocSingular_1.29.2            `\
`#> [105] vipor_0.4.7                     cli_3.6.6                      `\
`#> [107] rsvd_1.0.5                      rappdirs_0.3.4                 `\
`#> [109] textshaping_1.0.5               viridisLite_0.4.3              `\
`#> [111] S4Arrays_1.13.2                 dplyr_1.2.1                    `\
`#> [113] gtable_0.3.6                    yulab.utils_0.2.5              `\
`#> [115] sass_0.4.10                     digest_0.6.39                  `\
`#> [117] BiocGenerics_0.59.12            ggrepel_0.9.8                  `\
`#> [119] SparseArray_1.13.4              decontam_1.33.0                `\
`#> [121] htmlwidgets_1.6.4               farver_2.1.2                   `\
`#> [123] htmltools_0.5.9                 pkgdown_2.2.1                  `\
`#> [125] multtest_2.69.0                 lifecycle_1.0.5                `\
`#> [127] MASS_7.3-66`
