# dar 2.0.0

Changes since dar 1.8.0.

## Breaking changes

* R 4.6 or later is required. Differential-abundance recipes configure their
  targets and taxonomic rank with `add_model()`. The old selectors and recipe
  slots are removed. Rebuild saved recipes from their inputs.
* Use `step_aldex2()` for ALDEx2 and `recipe_qc()` for quality control.
  `step_aldex()`, `phy_qc()` and `step_metagenomeseq()` are removed.
  ALDEx2 effects now use the fitted GLM coefficient.
* Recipes use a versioned JSON exchange format. Importing bakes into an
  unprepared recipe requires `prepare = TRUE`.
* `otu_table()` returns counts with taxa in rows. Sample and taxonomy
  accessors return the full tables. ComplexHeatmap, gplots, heatmaply and
  UpSetR are optional.

## Models and results

* Shared models support covariate adjustment, longitudinal contrasts and
  compatible mixed designs. Added LinDA, ALDEx3 and radEmu, with per-stage
  `engine_args` for advanced options.
* `tidy_results()`, `tidy_abundance()`, `get_execution()` and `effect_plt()`
  expose per-method results, aggregated abundances, execution details and
  effect estimates.
* Abundance plots support shared or union discoveries and explicit taxa.
  Automatic selection requires a contrast when ambiguous. Consensus outputs
  retain the contrast and direction of each discovery.

## Fixes and maintenance

* Random-intercept models accept singleton participants alongside participants
  with repeated observations. Checks use observed groups and still reject
  designs with one group or no repeated observations (#189).
* Recipes validate steps before execution, run preprocessing in declared
  order, reserve result column names and assign unique automatic IDs (#187).
* DESeq2 uses native `poscounts` normalization and applies the significance
  threshold and result options before shrinkage (#188).
* Corrected contrast handling in MaAsLin3 and LEfSe, robust inference in
  ANCOM-BC2 and corncob, and depth-dependent zero handling in LinDA.
* Five vignettes document the current workflow. Docker Hub publication runs
  separately from GitHub package checks (#190).

# dar 1.5.6

## Breaking changes

* `step_metagenomeseq()` is defunct; its model implementation and related
  references are removed.

# dar 1.5.5

## Breaking changes

* `step_maaslin()` uses MaAsLin3 in place of MaAsLin2. Review method-specific
  arguments when migrating existing recipes.

# dar 1.5.4

* Parallel execution limits connection setup to support servers with more
  than 128 cores.
* Parallel execution uses a temporary future plan, preserving the user's
  global future configuration.

# dar 1.5.3

* ANCOM result tables join taxonomy by `taxa_id`; corncob, LEfSe and
  Wilcoxon tables place `taxa_id` first, restoring consensus extraction
  through `bake()` and `cool()` (#117).

# dar 1.1.3

* `step_ancom()` is temporarily disabled.

# dar 1.1.2

* The data import vignette is corrected.

# dar 0.99.13

* The HUMAnN example is removed from the data import vignette.

# dar 0.99.10

* `step_corncob()` is reintroduced after corncob returns to CRAN.

# dar 0.99.9

* `prep()` uses four workers to avoid issues with Bioconductor build
  infrastructure.

# dar 0.99.8

* Examples take less time to execute.

# dar 0.99.7

* Examples take less time to execute.

# dar 0.99.6

* Dependency installation guidance recommends `BiocManager::install()`.

# dar 0.99.5

* Vignettes take less time to execute.

# dar 0.99.4

* Tests and examples take less time to execute.
* The Linux GitHub Actions workflow correctly installs rlang.

# dar 0.99.3

## Breaking changes

* `Recipe` and `PrepRecipe` replace the former `recipe` and
  `prep_recipe` classes. Update code that inspects the old class names.
* `required_deps()` is internal and is no longer exported. Follow the
  installation guidance emitted by recipe steps for missing dependencies.
* R 4.4.0 or later is required.
* The package stops re-exporting `%>%` and `:=`. Use `|>` in
  recipes and qualify functions from their originating packages when needed.

## New features

* `recipe()` accepts both phyloseq and TreeSummarizedExperiment objects.
* Filtering and subsetting steps support chainable recipe construction.
* `step_filter_by_abundance()`, `step_filter_by_prevalence()`,
  `step_filter_by_rarity()` and `step_filter_by_variance()` provide
  dedicated filtering criteria.

## Documentation and improvements

* The import tutorial covers BIOM, QIIME, mothur, MetaPhlAn and HUMAnN
  conversions to phyloseq and TreeSummarizedExperiment.
* Installation guidance uses `BiocManager::install()` and
  `install.packages()` instead of unconventional installation methods.
* The data.table dependency is removed.
* Rarefaction warns about its fixed seed and the implications for results.

## Bug fixes

* `run_aldex()` corrects its internal seed handling.
* The data documentation uses the correct dataset name.

# dar 0.99.0

* Initial Bioconductor submission.
