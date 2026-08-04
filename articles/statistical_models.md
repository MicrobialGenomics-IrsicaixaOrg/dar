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
tidy_results(fit) |>
  dplyr::distinct(contrast_id, comparison, contrast_type)
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
engine-specific tables stored internally. The `method` column identifies
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
interfaces are deprecated in Bioconductor 3.24, become defunct in 3.25
and are scheduled for removal in 3.26. Model-free DA execution follows
the same schedule; preprocessing without a model remains supported. New
code should inspect `get_model(rec)$targets` and
`get_model(rec)$tax_level`. Legacy recipe objects loaded from RDS can be
migrated by moving both selectors into
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

The retired executable pseudo-JSON format cannot be imported safely.
Recreate those steps in a current recipe and export them with
[`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
to obtain the versioned, non-executable JSON schema.
[`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
is defunct in 3.24 and is replaced by
[`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md);
[`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
has completed its lifecycle and is no longer exported.

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
#> [1] dar_1.9.12
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6         xfun_0.60            bslib_0.12.0        
#>  [4] ggplot2_4.0.3        htmlwidgets_1.6.4    Biobase_2.73.2      
#>  [7] lattice_0.22-9       vctrs_0.7.3          tools_4.6.1         
#> [10] generics_0.1.4       biomformat_1.41.0    stats4_4.6.1        
#> [13] parallel_4.6.1       tibble_3.3.1         cluster_2.1.8.3     
#> [16] pkgconfig_2.0.3      Matrix_1.7-6         data.table_1.18.4   
#> [19] RColorBrewer_1.1-3   S7_0.2.2             desc_1.4.3          
#> [22] S4Vectors_0.51.6     lifecycle_1.0.5      compiler_4.6.1      
#> [25] farver_2.1.2         stringr_1.6.0        textshaping_1.0.5   
#> [28] Biostrings_2.81.6    Seqinfo_1.3.0        codetools_0.2-20    
#> [31] permute_0.9-10       htmltools_0.5.9      sass_0.4.10         
#> [34] yaml_2.3.12          pillar_1.11.1        pkgdown_2.2.1       
#> [37] crayon_1.5.3         jquerylib_0.1.4      MASS_7.3-66         
#> [40] cachem_1.1.0         vegan_2.7-5          iterators_1.0.14    
#> [43] foreach_1.5.2        nlme_3.1-170         tidyselect_1.2.1    
#> [46] digest_0.6.39        stringi_1.8.7        dplyr_1.2.1         
#> [49] reshape2_1.4.5       splines_4.6.1        ade4_1.7-24         
#> [52] fastmap_1.2.0        grid_4.6.1           cli_3.6.6           
#> [55] magrittr_2.0.5       survival_3.8-9       ape_5.8-1           
#> [58] scales_1.4.0         rmarkdown_2.31       XVector_0.53.0      
#> [61] igraph_2.3.3         multtest_2.69.0      otel_0.2.0          
#> [64] ragg_1.5.2           phyloseq_1.57.0      evaluate_1.0.5      
#> [67] knitr_1.51           IRanges_2.47.2       mgcv_1.9-4          
#> [70] rlang_1.3.0          Rcpp_1.1.2           glue_1.8.1          
#> [73] BiocGenerics_0.59.10 jsonlite_2.0.0       R6_2.6.1            
#> [76] plyr_1.8.9           systemfonts_1.3.2    fs_2.1.0
```
