# Statistical models and contrasts

[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
defines one statistical design shared by every compatible
differential-abundance step. Its `targets` and `tax_level` arguments are
the single source of truth for the analysis configuration; other fixed
terms provide covariate adjustment. Method-specific thresholds,
normalization and rarefaction settings remain active. A model-free
recipe remains valid for preprocessing, but differential-abundance
execution without a model is deprecated.

## Small example data

These simulated counts illustrate model configuration, not biological
findings, statistical performance or power. There are 12 subjects, each
measured at three visits. The cross-sectional example uses only their
baseline samples. Messages and warnings remain visible; this small
example does not establish estimation quality or validate a method’s
assumptions.

\
`metadata`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  visit ``=`` ``0``:``2``, subject_index ``=`` ``1``:``6``,`\
`  condition ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"control"``, ``"treated"``)``,`\
`  KEEP.OUT.ATTRS ``=`` ``FALSE``, stringsAsFactors ``=`` ``FALSE`\
`)`\
`metadata``$``subject`` ``<-`` `[`paste`](https://rdrr.io/r/base/paste.html)`(``metadata``$``condition``, ``metadata``$``subject_index``, sep ``=`` ``"_"``)`\
`metadata``$``batch`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``metadata``$``subject_index`` `[`%%`](https://rdrr.io/r/base/Arithmetic.html)` ``2`` ``==`` ``0``, ``"A"``, ``"B"``)``)`\
`metadata``$``age`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``25``, ``42``, ``31``, ``56``, ``38``, ``49``)``, each ``=`` ``3``, times ``=`` ``2``)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``metadata``)`` ``<-`` `[`paste`](https://rdrr.io/r/base/paste.html)`(``metadata``$``subject``, ``metadata``$``visit``, sep ``=`` ``"_"``)`\
`means`` ``<-`` `[`outer`](https://rdrr.io/r/base/outer.html)`(`[`seq`](https://rdrr.io/r/base/seq.html)`(``40``, ``150``, length.out ``=`` ``12``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``1``, `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``metadata``)``)``)`\
`means``[``1``, ``]`` ``<-`` ``means``[``1``, ``]`` ``*`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``metadata``$``condition`` ``==`` ``"treated"``, ``2``, ``1``)`\
`means``[``2``, ``]`` ``<-`` ``means``[``2``, ``]`` ``*`` `[`exp`](https://rdrr.io/r/base/Log.html)`(`\
`  ``0.3`` ``*`` ``metadata``$``visit`` ``*`` ``(``metadata``$``condition`` ``==`` ``"treated"``)`\
`)`\
`counts`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  `[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``means``)``, mu ``=`` `[`as.vector`](https://rdrr.io/r/base/vector.html)`(``means``)``, size ``=`` ``20``)``,`\
`  nrow ``=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``means``)``,`\
`  dimnames ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(``"taxon_"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``means``)``)``)``, `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``metadata``)``)`\
`)`\
`taxonomy`` ``<-`` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(`\
`  Genus ``=`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"Genus_"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``counts``)``)``)``,`\
`  Species ``=`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"Species_"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``counts``)``)``)`\
`)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``taxonomy``)`` ``<-`` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`\
`longitudinal_phy`` ``<-`` ``phyloseq``::`[`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)`(`\
`  ``phyloseq``::`[`otu_table`](https://rdrr.io/pkg/phyloseq/man/otu_table-methods.html)`(``counts``, taxa_are_rows ``=`` ``TRUE``)``,`\
`  ``phyloseq``::`[`tax_table`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)`(``taxonomy``)``,`\
`  ``phyloseq``::`[`sample_data`](https://rdrr.io/pkg/phyloseq/man/sample_data-methods.html)`(``metadata``)`\
`)`\
`phy`` ``<-`` ``phyloseq``::`[`prune_samples`](https://rdrr.io/pkg/phyloseq/man/prune_samples-methods.html)`(``metadata``$``visit`` ``==`` ``0``, ``longitudinal_phy``)`

## Cross-sectional adjustment

\
`rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``phy``)`` ``|>`\
`  `[`add_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)`(`\
`    ``~`` ``condition`` ``+`` ``batch`` ``+`` ``age``,`\
`    targets ``=`` ``"condition"``,`\
`    tax_level ``=`` ``"Species"``,`\
`    reference ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``condition ``=`` ``"control"``)``,`\
`    na_action ``=`` ``"error"`\
`  ``)`` ``|>`\
`  `[`step_deseq`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md)`(``)`\
\
[`get_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)`(``rec``)`\
`#> $formula`\
`#> ~condition + batch + age`\
`#> `\
`#> $targets`\
`#> [1] "condition"`\
`#> `\
`#> $tax_level`\
`#> [1] "Species"`\
`#> `\
`#> $time`\
`#> NULL`\
`#> `\
`#> $reference`\
`#> condition `\
`#> "control" `\
`#> `\
`#> $contrasts`\
`#> [1] "simple"`\
`#> `\
`#> $na_action`\
`#> [1] "error"`\
`#> `\
`#> $max_time_levels`\
`#> [1] 12`\
`cross_sectional_fit`` ``<-`` `[`prep`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)`(``rec``, parallel ``=`` ``FALSE``)`

The reference is the denominator, so a positive effect for
`condition[treated-control]` means higher abundance in `treated`. With
`na_action = "drop"`,
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
removes one common complete-case cohort before any engine runs and
records the sample identifiers in the execution manifest. Use
`na_action = "error"` when missing model variables should stop the
analysis instead of changing the cohort.

## Longitudinal adjustment

Numeric time values are converted internally to a factor whose levels
are sorted by numeric value. The original sample metadata are not
modified.

\
`longitudinal`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``longitudinal_phy``)`` ``|>`\
`  `[`add_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)`(`\
`    ``~`` ``condition`` ``*`` ``visit`` ``+`` ``batch`` ``+`` ``(``1`` ``|`` ``subject``)``,`\
`    targets ``=`` ``"condition"``,`\
`    tax_level ``=`` ``"Species"``,`\
`    time ``=`` ``"visit"``,`\
`    reference ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``condition ``=`` ``"control"``, visit ``=`` ``"0"``)``,`\
`    contrasts ``=`` ``"all"``,`\
`    na_action ``=`` ``"drop"`\
`  ``)`` ``|>`\
`  `[`step_linda`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md)`(``)`\
\
`fit`` ``<-`` `[`prep`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)`(``longitudinal``, parallel ``=`` ``FALSE``)`\
[`tidy_results`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_results.md)`(``fit``)`` ``|>`\
`  ``dplyr``::`[`distinct`](https://dplyr.tidyverse.org/reference/distinct.html)`(``contrast_id``, ``comparison``, ``contrast_type``)`\
`#> ``# A tibble: 12 × 3`\
`#>    contrast_id                           comparison                contrast_type`\
`#>    ``<chr>``                                 ``<chr>``                     ``<chr>``        `\
`#> `` 1`` condition[treated-control]@visit[0]   condition[treated-contro… condition_wi…`\
`#> `` 2`` condition[treated-control]@visit[1]   condition[treated-contro… condition_wi…`\
`#> `` 3`` condition[treated-control]@visit[2]   condition[treated-contro… condition_wi…`\
`#> `` 4`` visit[1-0]@condition[control]         visit[1-0]@condition[con… time_within_…`\
`#> `` 5`` visit[2-0]@condition[control]         visit[2-0]@condition[con… time_within_…`\
`#> `` 6`` visit[2-1]@condition[control]         visit[2-1]@condition[con… time_within_…`\
`#> `` 7`` visit[1-0]@condition[treated]         visit[1-0]@condition[tre… time_within_…`\
`#> `` 8`` visit[2-0]@condition[treated]         visit[2-0]@condition[tre… time_within_…`\
`#> `` 9`` visit[2-1]@condition[treated]         visit[2-1]@condition[tre… time_within_…`\
`#> ``10`` condition[treated-control]:visit[1-0] condition[treated-contro… difference_i…`\
`#> ``11`` condition[treated-control]:visit[2-0] condition[treated-contro… difference_i…`\
`#> ``12`` condition[treated-control]:visit[2-1] condition[treated-contro… difference_i…`

For `condition * visit`, the simple plan contains condition comparisons
within every visit and all pairwise visit comparisons within every
condition. Later time values are numerators. `contrasts = "all"`
additionally requests difference-in-differences. More than 12 time
points requires an explicit larger `max_time_levels`, preventing
accidental contrast explosions.

`(1 | subject)` specifies a subject-level random intercept for repeated
measurements. Declaring `time` does not itself add a random effect or
establish that observations are independent. Include the subject
structure required by the study design and choose a compatible method.

## Engine capabilities

| Engine    | Fixed adjustment | Time interaction | Random effects |
|:----------|:----------------:|:----------------:|:--------------:|
| DESeq2    |       yes        |       yes        |       no       |
| ALDEx2    |       yes        |       yes        |       no       |
| ANCOM-BC2 |       yes        |       yes        |      yes       |
| corncob   |       yes        |       yes        |       no       |
| LinDA     |       yes        |       yes        |      yes       |
| MaAsLin3  |       yes        |       yes        |      yes       |
| Wilcoxon  |        no        |        no        |       no       |
| LEfSe     |        no        |        no        |       no       |

Fixed adjustment includes an additive time term, even when that column
is declared with `time`. Wilcoxon and LEfSe support one categorical
target without adjustment; they cannot execute `~ condition + time` with
`targets = "condition"`.

Incompatible steps are warned about when the chain is built and
rechecked after preprocessing.
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
records them rather than treating them as failed fits. If no method is
compatible, preparation stops. Real engine failures and missing planned
contrasts remain errors.

Use `contrast_id`, `comparison` and `contrast_type` to identify the
hypothesis represented by each result. See [Exploring results and
optional
consensus](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/article.md)
for the result columns, effect scales, selectors and plotting examples.

For centralized ALDEx2 analyses, `effect` is the contrast-oriented
coefficient from the same GLM as the adjusted P value, averaged over
Monte Carlo instances and labeled `aldex_glm_coefficient`. With default
settings it is on the base-2 CLR response scale, not an
absolute-abundance fold change or ALDEx2’s standardized two-group
effect. Advanced native settings can change its scale. Recompute
previously prepared centralized ALDEx2 results to obtain the GLM
coefficients; deprecated model-free analyses retain their standardized
effects.

## Advanced native-engine arguments

The eight DA steps accept `engine_args`, a named list grouped by
native-call stage. This keeps advanced settings attached to the
configured step and avoids an ambiguous `...`. For example:

\
`advanced`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``phy``)`` ``|>`\
`  `[`add_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)`(`\
`    ``~`` ``condition`` ``+`` ``batch`` ``+`` ``age``, targets ``=`` ``"condition"``, tax_level ``=`` ``"Species"``,`\
`    reference ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``condition ``=`` ``"control"``)`\
`  ``)`` ``|>`\
`  `[`step_deseq`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md)`(`\
`    engine_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`      fit ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``minReplicatesForReplace ``=`` ``Inf``)``,`\
`      shrink ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``lfcThreshold ``=`` ``1``)`\
`    ``)`\
`  ``)`` ``|>`\
`  `[`step_ancom`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md)`(`\
`    engine_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`      fit ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``iter_control ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``tol ``=`` ``1e-3``, max_iter ``=`` ``50``, verbose ``=`` ``FALSE``)``)`\
`    ``)`\
`  ``)`\
`#> Registered S3 method overwritten by 'httr':`\
`#>   method         from  `\
`#>   print.response rmutil`\
`advanced`\
`#> ── ``DAR Recipe`` ──────────────────────────────────────────────────────────────────`\
`#> Inputs:`\
`#> `\
`#>      ``ℹ`` phyloseq object with ``12`` taxa and ``12`` samples `\
`#>      ``ℹ`` variable of interes ``condition`` (class: character, levels: control, treated) `\
`#>      ``ℹ`` taxonomic level ``Species`` `\
`#> `\
`#> Statistical model:`\
`#> `\
`#>      ``ℹ`` ``~condition + batch + age`` `\
`#> `\
`#> Preprocessing steps:`\
`#> `\
`#> `\
`#> DA steps:`\
`#> `\
`#>      ``◉`` step_deseq() ``id = deseq__Burek`` `\
`#>      ``◉`` step_ancom() ``id = ancom__Marillenknödel`

Nested native control lists replace the corresponding native argument;
they are not merged with its defaults. For ANCOM-BC2, supply the
complete `iter_control` list, including `verbose`, rather than only the
values changed.

Available stages are `size_factors`, `fit`, `results`, and `shrink` for
DESeq2; `clr` and `fit` for ALDEx2; `fit` for ANCOM-BC2, corncob, and
LinDA; `fit` and `contrast` for MaAsLin3; `test` for Wilcoxon; and
`transform` and `fit` for LEfSe. Data, formulas, contrasts, output paths
and arguments already exposed by each dar wrapper are protected.
Non-empty `engine_args` require a centralized model. A configured stage
that is not used by the selected execution route is reported as an
error. ALDEx2’s former `effect` stage is rejected because centralized
effects now come directly from the GLM fit.

See [Saving and sharing analysis
recipes](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/import_export_recipes.md)
for exporting model settings and native options without executing an
analysis.

## Migrating legacy selectors

The former `recipe(..., var_info=, tax_info=)`,
[`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
[`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md),
[`get_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_var.md)
and
[`get_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_tax.md)
interfaces, and model-free DA execution, are deprecated in Bioconductor
3.24, become defunct in 3.25 and are scheduled for removal in 3.26.
Preprocessing without a model remains supported. Move the selectors into
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
and inspect them with `get_model(rec)`:

\
`# Deprecated`\
`old`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``phy``, var_info ``=`` ``"condition"``, tax_info ``=`` ``"Species"``)`\
\
`# Current API`\
`current`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``phy``)`` ``|>`\
`  `[`add_model`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)`(`\
`    ``~`` ``condition``,`\
`    targets ``=`` ``"condition"``,`\
`    tax_level ``=`` ``"Species"`\
`  ``)`

[`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
is defunct in 3.24; use
[`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md).
The retired
[`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
is no longer exported. Legacy executable recipe files must be recreated
in the current API before export; they are not accepted by the current
importer.

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
`#>   [1] fs_2.1.0                        matrixStats_1.5.0              `\
`#>   [3] DirichletMultinomial_1.55.1     doParallel_1.0.17              `\
`#>   [5] httr_1.4.9                      RColorBrewer_1.1-3             `\
`#>   [7] numDeriv_2016.8-1.1             doRNG_1.8.6.3                  `\
`#>   [9] tools_4.6.1                     backports_1.5.1                `\
`#>  [11] utf8_1.2.6                      R6_2.6.1                       `\
`#>  [13] vegan_2.7-6                     lazyeval_0.2.3                 `\
`#>  [15] mgcv_1.9-4                      apeglm_1.35.0                  `\
`#>  [17] permute_0.9-10                  withr_3.0.3                    `\
`#>  [19] gridExtra_2.3.1                 cli_3.6.6                      `\
`#>  [21] Biobase_2.73.2                  textshaping_1.0.5              `\
`#>  [23] sandwich_3.1-3                  sass_0.4.10                    `\
`#>  [25] SQUAREM_2026.1                  mvtnorm_1.4-2                  `\
`#>  [27] S7_0.2.2                        readr_2.2.0                    `\
`#>  [29] proxy_0.4-29                    pkgdown_2.2.1                  `\
`#>  [31] mixsqp_0.3-54                   mlr3tuning_1.7.0               `\
`#>  [33] systemfonts_1.3.2               yulab.utils_0.2.5              `\
`#>  [35] foreign_0.8-91                  paradox_1.0.1                  `\
`#>  [37] scater_1.41.2                   decontam_1.33.0                `\
`#>  [39] parallelly_1.48.0               invgamma_1.2                   `\
`#>  [41] bbmle_1.0.26                    readxl_1.5.0.1                 `\
`#>  [43] rstudioapi_0.19.0               generics_0.1.4                 `\
`#>  [45] gtools_3.9.5                    dplyr_1.2.1                    `\
`#>  [47] Matrix_1.7-6                    biomformat_1.41.0              `\
`#>  [49] ggbeeswarm_0.7.3                DescTools_0.99.60              `\
`#>  [51] S4Vectors_0.51.10               DECIPHER_3.9.3                 `\
`#>  [53] abind_1.4-8                     lifecycle_1.0.5                `\
`#>  [55] multcomp_1.4-32                 yaml_2.3.12                    `\
`#>  [57] SummarizedExperiment_1.43.0     SparseArray_1.13.4             `\
`#>  [59] grid_4.6.1                      crayon_1.5.3                   `\
`#>  [61] bdsmatrix_1.3-7                 lattice_0.23-1                 `\
`#>  [63] haven_2.5.5                     beachmat_2.29.3                `\
`#>  [65] pillar_1.11.1                   knitr_1.52                     `\
`#>  [67] GenomicRanges_1.65.4            statip_0.2.3                   `\
`#>  [69] boot_1.3-32                     gld_2.6.8                      `\
`#>  [71] codetools_0.2-20                glue_1.8.1                     `\
`#>  [73] spacefillr_0.4.0                data.table_1.18.6.1            `\
`#>  [75] MultiAssayExperiment_1.39.1     vctrs_0.7.3                    `\
`#>  [77] treeio_1.37.1                   Rdpack_2.6.6                   `\
`#>  [79] cellranger_1.1.0                gtable_0.3.6                   `\
`#>  [81] emdbook_1.3.14                  cachem_1.1.0                   `\
`#>  [83] xfun_0.61                       rbibutils_2.4.1                `\
`#>  [85] S4Arrays_1.13.2                 Seqinfo_1.3.2                  `\
`#>  [87] coda_0.19-4.1                   reformulas_0.4.4               `\
`#>  [89] modeest_2.5.0                   survival_3.8-12                `\
`#>  [91] timeDate_4052.112               SingleCellExperiment_1.35.2    `\
`#>  [93] iterators_1.0.14                statmod_1.5.2                  `\
`#>  [95] bluster_1.23.1                  TH.data_1.1-5                  `\
`#>  [97] nlme_3.1-171                    ANCOMBC_2.15.1                 `\
`#>  [99] phyloseq_1.57.0                 bbotk_1.13.1                   `\
`#> [101] fBasics_4052.98                 bslib_0.12.0                   `\
`#> [103] MicrobiomeStat_1.4              mlr3_1.8.0                     `\
`#> [105] irlba_2.4.1                     vipor_0.4.7                    `\
`#> [107] otel_0.2.0                      rpart_4.1.27                   `\
`#> [109] colorspace_2.1-4                Hmisc_5.3-0                    `\
`#> [111] mlr3mbo_1.2.1                   moocore_0.3.2                  `\
`#> [113] BiocGenerics_0.59.12            DBI_1.3.0                      `\
`#> [115] nnet_7.3-21                     ade4_1.7-24                    `\
`#> [117] mlr3misc_0.23.0                 Exact_3.3                      `\
`#> [119] DESeq2_1.53.6                   tidyselect_1.2.1               `\
`#> [121] timeSeries_4052.112             compiler_4.6.1                 `\
`#> [123] htmlTable_2.5.0                 BiocNeighbors_2.7.3            `\
`#> [125] lgr_0.5.2                       expm_1.0-1                     `\
`#> [127] desc_1.4.3                      DelayedArray_0.39.8            `\
`#> [129] checkmate_2.3.4                 scales_1.4.0                   `\
`#> [131] quadprog_1.5-8                  spatial_7.3-19                 `\
`#> [133] rappdirs_0.3.4                  palmerpenguins_0.1.1           `\
`#> [135] stringr_1.6.0                   digest_0.6.39                  `\
`#> [137] minqa_1.2.8                     rmarkdown_2.32                 `\
`#> [139] XVector_0.53.0                  base64enc_0.1-6                `\
`#> [141] htmltools_0.5.9                 pkgconfig_2.0.3                `\
`#> [143] lme4_2.0-6                      sparseMatrixStats_1.25.0       `\
`#> [145] MatrixGenerics_1.25.0           stabledist_0.7-2               `\
`#> [147] fastmap_1.2.0                   rlang_1.3.0                    `\
`#> [149] htmlwidgets_1.6.4               DelayedMatrixStats_1.35.0      `\
`#> [151] farver_2.1.2                    jquerylib_0.1.4                `\
`#> [153] zoo_1.9-1                       energy_1.7-12                  `\
`#> [155] jsonlite_2.0.0                  BiocParallel_1.47.0            `\
`#> [157] BiocSingular_1.29.2             magrittr_2.0.5                 `\
`#> [159] Formula_1.2-6                   scuttle_1.23.2                 `\
`#> [161] Rcpp_1.1.2                      ape_5.8-1                      `\
`#> [163] viridis_0.6.5                   stringi_1.8.9                  `\
`#> [165] rootSolve_1.8.2.4               stable_1.1.7                   `\
`#> [167] MASS_7.3-66                     plyr_1.8.9                     `\
`#> [169] parallel_4.6.1                  listenv_1.1.0                  `\
`#> [171] ggrepel_0.9.8                   forcats_1.0.1                  `\
`#> [173] lmom_3.3                        Biostrings_2.81.9              `\
`#> [175] splines_4.6.1                   multtest_2.69.0                `\
`#> [177] hms_1.1.4                       locfit_1.5-9.12                `\
`#> [179] igraph_2.3.4                    uuid_1.2-2                     `\
`#> [181] rngtools_1.5.2                  reshape2_1.4.5                 `\
`#> [183] stats4_4.6.1                    ScaledMatrix_1.21.0            `\
`#> [185] rmutil_1.1.10                   evaluate_1.0.5                 `\
`#> [187] nloptr_2.2.1                    tzdb_0.5.0                     `\
`#> [189] foreach_1.5.2                   tidyr_1.3.2                    `\
`#> [191] purrr_1.2.2                     future_1.76.0                  `\
`#> [193] clue_0.3-68                     ggplot2_4.0.3                  `\
`#> [195] ashr_2.2-63                     rsvd_1.0.5                     `\
`#> [197] e1071_1.7-17                    tidytree_0.4.8                 `\
`#> [199] viridisLite_0.4.3               class_7.3-24                   `\
`#> [201] ragg_1.5.2                      gsl_2.1-9                      `\
`#> [203] truncnorm_1.0-9                 tibble_3.3.1                   `\
`#> [205] lmerTest_3.2-1                  beeswarm_0.4.0                 `\
`#> [207] IRanges_2.47.5                  cluster_2.1.8.3                `\
`#> [209] TreeSummarizedExperiment_2.21.0 globals_0.20.0                 `\
`#> [211] mia_1.21.11`
