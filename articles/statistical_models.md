# Central statistical models and longitudinal contrasts

[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
defines one statistical design shared by every compatible
differential-abundance step. Its `targets` and `tax_level` arguments are
the single source of truth for the analysis configuration; other fixed
terms adjust for confounding. Method-specific thresholds, normalization
and rarefaction settings remain active. A model-free recipe remains
valid for preprocessing, but differential-abundance execution without a
model is deprecated.

## Cross-sectional adjustment

``` r

rec <- recipe(phy) |>
  add_model(
    ~ condition + batch + age,
    targets = "condition",
    tax_level = "Species",
    reference = c(condition = "control"),
    na_action = "error"
  ) |>
  step_deseq() |>
  step_aldex() |>
  step_ancom() |>
  step_linda()

get_model(rec)
fit <- prep(rec)
```

The reference is the denominator, so a positive effect for
`condition[treated-control]` means higher abundance in `treated`. With
`na_action = "drop"`,
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
removes one common complete-case cohort before any engine runs and
records the sample identifiers in the execution manifest.

## Longitudinal adjustment

Numeric time values are converted internally to a factor whose levels
are sorted by numeric value. The original sample metadata are not
modified.

``` r

longitudinal <- recipe(longitudinal_phy) |>
  add_model(
    ~ condition * visit + batch + (1 | subject),
    targets = "condition",
    tax_level = "Species",
    time = "visit",
    reference = c(condition = "control", visit = "0"),
    contrasts = "all",
    na_action = "drop"
  ) |>
  step_ancom() |>
  step_linda() |>
  step_maaslin()

fit <- prep(longitudinal)
fit@execution$contrasts
```

For `condition * visit`, the simple plan contains condition comparisons
within every visit and all pairwise visit comparisons within every
condition. Later time values are numerators. `contrasts = "all"`
additionally requests difference-in-differences. More than 12 time
points requires an explicit larger `max_time_levels`, preventing
accidental contrast explosions.

## Engine capabilities

| Engine    | Fixed confounders | Time interaction | Random effects |
|:----------|:-----------------:|:----------------:|:--------------:|
| DESeq2    |        yes        |       yes        |       no       |
| ALDEx2    |        yes        |       yes        |       no       |
| ANCOM-BC2 |        yes        |       yes        |      yes       |
| corncob   |        yes        |       yes        |       no       |
| LinDA     |        yes        |       yes        |      yes       |
| MaAsLin3  |        yes        |       yes        |      yes       |
| Wilcoxon  |        no         |        no        |       no       |
| LEfSe     |        no         |        no        |       no       |

Incompatible steps are warned about when the chain is built and
rechecked after preprocessing.
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
records them rather than treating them as failed fits. If no method is
compatible, preparation stops. Real engine failures and missing planned
contrasts remain errors.

Prepared modeled results contain `contrast_id`, `comparison`,
`contrast_type` and `var`. Consensus operations group by `contrast_id`
and effect direction, so different longitudinal hypotheses cannot
contribute to the same vote.

## A uniform tidy result contract

[`tidy_results()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_results.md)
exposes the same columns for every DA engine while preserving the
engine-specific tables in `fit@results`. The `method` column identifies
the stable engine (`deseq`, `aldex`, and so on), whereas `step_id`
identifies the configured instance and therefore distinguishes repeated
configurations of one engine.

``` r

results <- tidy_results(
  fit,
  steps = c("ancom__primary", "maaslin__sensitivity"),
  significant_only = TRUE
)

results |>
  dplyr::select(
    taxa_id, contrast_id, comparison, effect_size, effect_metric,
    adj_p_value, significant, method, step_id
  )
```

Effect values retain their native scientific scale. Use `effect_metric`
when interpreting them: for example, a DESeq2 log2 fold change and an
ALDEx2 standardized CLR effect are labelled explicitly but are not
converted to a common scale and should not be compared numerically.

## Advanced native-engine arguments

The eight DA steps accept `engine_args`, a named list grouped by
native-call stage. This keeps advanced settings attached to the
configured step and avoids an ambiguous `...`. For example:

``` r

rec |>
  step_deseq(
    engine_args = list(
      fit = list(minReplicatesForReplace = Inf),
      shrink = list(lfcThreshold = 1)
    )
  ) |>
  step_ancom(
    engine_args = list(
      fit = list(iter_control = list(tol = 1e-3, max_iter = 50))
    )
  )
```

Available stages are `size_factors`, `fit`, `results`, and `shrink` for
DESeq2; `clr`, `fit`, and `effect` for ALDEx2; `fit` for ANCOM-BC2,
corncob, and LinDA; `fit` and `contrast` for MaAsLin3; `test` for
Wilcoxon; and `transform` and `fit` for LEfSe. Data, formulas,
contrasts, output paths and arguments already exposed by each dar
wrapper are protected. Non-empty `engine_args` require a centralized
model. A configured stage that is not used by the selected execution
route is reported as an error.

[`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
and
[`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
preserve the formula, targets, taxonomic level, time column, references,
contrast mode, missing-value policy and time-level limit in the
versioned `dar-recipe` JSON schema. Executable files created by the
legacy exporter are rejected with migration guidance.

## Migrating legacy selectors

The former `recipe(..., var_info=, tax_info=)`,
[`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
[`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md),
[`get_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_var.md)
and
[`get_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_tax.md)
interfaces remain available during the first Bioconductor deprecation
cycle. They emit classified warnings and synchronize with the
centralized model where applicable. New code should inspect
`get_model(rec)$targets` and `get_model(rec)$tax_level` instead. Legacy
recipes can be migrated by moving both selectors into
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md):

``` r

# Deprecated
old <- recipe(phy, var_info = "condition", tax_info = "Species")

# Current API
current <- recipe(phy) |>
  add_model(
    ~ condition,
    targets = "condition",
    tax_level = "Species"
  )
```

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dar_1.9.9
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-24         tidyselect_1.2.1    viridisLite_0.4.3  
#>  [4] dplyr_1.2.1         farver_2.1.2        viridis_0.6.5      
#>  [7] Biostrings_2.80.1   S7_0.2.2            fastmap_1.2.0      
#> [10] TSP_1.2.7           phyloseq_1.56.0     digest_0.6.39      
#> [13] lifecycle_1.0.5     cluster_2.1.8.3     survival_3.8-9     
#> [16] magrittr_2.0.5      compiler_4.6.1      rlang_1.3.0        
#> [19] sass_0.4.10         tools_4.6.1         igraph_2.3.3       
#> [22] yaml_2.3.12         data.table_1.18.4   knitr_1.51         
#> [25] htmlwidgets_1.6.4   plyr_1.8.9          RColorBrewer_1.1-3 
#> [28] registry_0.5-1      ca_0.71.1           purrr_1.2.2        
#> [31] BiocGenerics_0.58.1 desc_1.4.3          grid_4.6.1         
#> [34] stats4_4.6.1        multtest_2.68.0     biomformat_1.40.0  
#> [37] ggplot2_4.0.3       scales_1.4.0        iterators_1.0.14   
#> [40] MASS_7.3-66         cli_3.6.6           vegan_2.7-5        
#> [43] UpSetR_1.4.1        rmarkdown_2.31      crayon_1.5.3       
#> [46] ragg_1.5.2          generics_0.1.4      otel_0.2.0         
#> [49] heatmaply_1.6.0     httr_1.4.8          reshape2_1.4.5     
#> [52] ape_5.8-1           cachem_1.1.0        stringr_1.6.0      
#> [55] splines_4.6.1       assertthat_0.2.1    parallel_4.6.1     
#> [58] XVector_0.52.0      vctrs_0.7.3         webshot_0.5.5      
#> [61] Matrix_1.7-6        jsonlite_2.0.0      seriation_1.5.8    
#> [64] IRanges_2.46.0      S4Vectors_0.50.1    systemfonts_1.3.2  
#> [67] dendextend_1.19.1   foreach_1.5.2       plotly_4.12.1      
#> [70] tidyr_1.3.2         jquerylib_0.1.4     glue_1.8.1         
#> [73] pkgdown_2.2.1       codetools_0.2-20    stringi_1.8.7      
#> [76] gtable_0.3.6        tibble_3.3.1        pillar_1.11.1      
#> [79] htmltools_0.5.9     Seqinfo_1.2.0       R6_2.6.1           
#> [82] textshaping_1.0.5   evaluate_1.0.5      Biobase_2.72.0     
#> [85] lattice_0.22-9      bslib_0.12.0        Rcpp_1.1.2         
#> [88] permute_0.9-10      gridExtra_2.3.1     nlme_3.1-170       
#> [91] mgcv_1.9-4          xfun_0.60           fs_2.1.0           
#> [94] pkgconfig_2.0.3
```
