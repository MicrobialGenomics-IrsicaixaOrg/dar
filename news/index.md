# Changelog

## dar 1.9.12

### Breaking changes

- R 4.6 or later is required.
- [`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  is removed after its defunct release cycle. Replace it with a
  supported DA step when rebuilding legacy recipes
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).
- DA step IDs must not use reserved consensus column names such as
  `effect` or `sum`. Rename conflicting IDs, for example to `da_effect`,
  before constructing or importing a recipe. Other IDs, including those
  beginning with a dot, remain unchanged.
- ComplexHeatmap, gplots, heatmaply and UpSetR are optional
  dependencies. Install the relevant package explicitly to use its
  plotting functions; missing dependencies produce installation guidance
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).

### Deprecations

- Legacy selectors and model-free DA execution remain available in
  Bioconductor 3.24, become defunct in 3.25 and are scheduled for
  removal in 3.26. Use
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  to configure DA analyses and
  [`get_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)
  to inspect their targets and taxonomic rank
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).

### Bug fixes

- Modeled engines, consensus tables and abundance plots share one
  aggregated taxonomic view with stable IDs. Abundance transformations
  use aggregated counts, and public accessors retain the complete
  underlying data.
- [`abundance_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/abundance_plt.md),
  [`intersection_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/intersection_plt.md)
  and
  [`mutual_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/mutual_plt.md)
  check selectors against planned model hypotheses, including those with
  no significant taxa. For empty selections, discovery plots report a
  message and return invisible `NULL`. Abundance plots with explicit
  `taxa_ids` remain available regardless of significance.
- [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  override warnings are emitted only for explicitly supplied
  engine-specific designs; default ANCOM-BC2 and corncob steps remain
  quiet.
- [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  emits the DA rarefaction caution when rarefying steps execute, rather
  than when they are constructed.
- [`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md)
  preserves stable taxon IDs after ANCOM-BC2 aggregation.
- [`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md)
  uses ANCOM-BC2’s robust significance decision when pseudocount
  sensitivity analysis is enabled, in both centralized and deprecated
  model-free execution. Explicitly disabling sensitivity analysis
  retains the basic native decision.
- [`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md)
  honors the requested absolute effect threshold in deprecated
  model-free execution.
- [`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md)
  uses native sandwich covariance for centralized Wald contrasts when
  `robust = TRUE` and forwards the requested `B` during deprecated
  model-free bootstrap execution.
- [`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md)
  locally corrects MicrobiomeStat 1.4’s adaptive zero-handling
  capitalization error, restoring depth-dependent imputation without
  modifying the installed dependency or the remaining native algorithm.
- [`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)
  executes allowed [`median()`](https://rdrr.io/r/stats/median.html),
  [`var()`](https://rdrr.io/r/stats/cor.html) and
  [`sd()`](https://rdrr.io/r/stats/sd.html) predicates locally and after
  JSON import, using fixed trusted bindings without expanding the
  permitted syntax or capturing external functions.
- [`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md)
  joins modeled results by taxon ID, keeping identically named taxa from
  different lineages distinct.
- [`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)
  and
  [`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md)
  skip unsupported additive time adjustment, including when the time
  column is explicitly declared.

### Documentation

- Function examples illustrate supported workflows without redundant
  setup, repeated defaults or intentional error demonstrations.
- [`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)
  documents its single, restricted serializable predicate and filtering
  during
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md),
  rather than arbitrary function lists or a `prune` argument.
- Six focused vignettes cover getting started, data inputs,
  preprocessing, models and contrasts, results and optional consensus,
  and recipe sharing. The duplicate case study is merged into the
  results guide. Examples use the central model API, execute adjusted
  and longitudinal designs, and import a bundled BIOM file into both
  supported input classes.
- Results guides start with per-method findings and distinguish
  agreement from calibrated error control and comparable effect scales.
  Input guidance distinguishes relative abundances from raw counts; the
  advanced ANCOM-BC2 example supplies a complete native control list.

### Development infrastructure

- The development-only Monte Carlo validation suite and its dedicated
  workflow are archived in the frozen
  `codex/scientific-validation-frozen-20261015` branch and removed from
  `devel`. The historical baseline is preserved without recalibration or
  relabeling. Previously generated reports are retained separately.
- Statistical correctness and native-engine regression tests remain in
  `tests/testthat`, alongside the ordinary package checks and vignette
  builds.

## dar 1.9.11

### Bug fixes

- Consensus tables, plots and bake results keep modeled targets,
  contrasts and effect directions as separate hypotheses
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Modeled intersections use the complete taxon-by-contrast universe,
  including non-significant hypotheses, and classify zero effects as
  neutral
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Method exclusions are applied before consensus support is counted.
  Invalid cutoffs, weights and exclusions produce classified errors
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Modeled exploration functions accept explicit `target` and
  `contrast_id` selectors; abundance plots require them when automatic
  selection is ambiguous
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).

## dar 1.9.10

### Breaking changes

- [`otu_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)
  always returns counts with taxa in rows and deterministic taxon and
  sample ordering. Adapt code that assumes the stored phyloseq
  orientation
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).
- [`sample_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)
  and
  [`tax_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)
  expose complete recipe metadata and taxonomy instead of
  model-restricted views. Use
  [`get_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)
  to inspect the analysis configuration, and select columns explicitly
  when needed
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).

### Improvements

- DA engines and plots use explicit internal target-only and
  selected-rank views, keeping analysis inputs independent of public
  accessor semantics
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).

## dar 1.9.9

### Breaking changes

- [`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
  and
  [`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
  use a versioned, non-executable JSON schema. Legacy pseudo-JSON files
  are rejected. Recreate trusted legacy recipes with the current API and
  export them again
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).
- [`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
  requires explicit `prepare = TRUE` to execute an analysis when
  importing bake configurations into an unprepared recipe
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).

### Improvements

- Imported constructors, arguments, formulas, nested values and filter
  predicates are allowlisted and validated transactionally before the
  recipe is changed
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).

## dar 1.9.8

- Recipe steps are classified by their concrete classes, so custom IDs
  do not change preprocessing, DA or bake semantics
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).
- Preprocessing operations execute in their configured order, including
  recipes that interleave filtering, subsetting and rarefaction
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).
- Unsupported, misplaced and duplicate-ID steps fail with classified
  validation errors before execution
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).

## dar 1.9.7

### Breaking changes

- [`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  is defunct. Use
  [`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
  with an explicit `group_by` when grouped quality-control metrics are
  required
  ([\#112](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/112)).

### New features

- [`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
  provides model-independent overall or explicitly grouped
  quality-control metrics for `Recipe` and `PrepRecipe` objects.
  Summaries expose unambiguous `group_by` and `group` columns and handle
  either count-table orientation
  ([\#112](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/112)).

## dar 1.9.6

- [`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md)
  integrates bias-corrected LinDA linear and mixed-effects models with
  every compatible centralized contrast. LinDA participates in the tidy
  result contract and Monte Carlo validation with explicit log2
  fold-change semantics
  ([\#128](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/128)).

## dar 1.9.5

- The seven centralized DA engines accept structured `engine_args` for
  advanced arguments at explicit native-call stages
  ([\#126](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/126)).
- Step export and import preserve nested advanced arguments and reject
  values that cannot be reconstructed faithfully
  ([\#126](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/126)).

## dar 1.9.4

- [`tidy_results()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_results.md)
  provides a stable, engine-labelled result contract while preserving
  the native tables in `PrepRecipe@results`
  ([\#127](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/127)).
- Result consumers validate taxon-contrast-step keys and use explicit
  effect metrics and p-value mappings
  ([\#127](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/127)).

## dar 1.9.3

### Breaking changes

- [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  becomes the single source of truth for DA targets and taxonomic
  resolution. Supply `targets` and `tax_level` there
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).
- [`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
  [`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md),
  [`get_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_var.md),
  [`get_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_tax.md),
  the `var_info` and `tax_info` arguments to
  [`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md),
  and model-free DA execution are deprecated. They remain functional
  during the first Bioconductor deprecation cycle; migrate to
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  and
  [`get_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).

### Improvements

- Legacy centralized models, serialized recipes and imported step files
  are normalized to the new model representation without changing
  original microbiome metadata
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).

### Bug fixes

- Modeled exclusion and mutual plots preserve contrast and
  effect-direction keys, and normalize adjusted p-value aliases without
  duplicate columns
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).
- LEfSe executes every planned pairwise contrast for multilevel targets
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).

## dar 1.9.2

### New features

- [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  and
  [`get_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md)
  define and inspect a shared statistical design with fixed confounders,
  target-time interactions, supported random effects, explicit
  references and a common missing-value policy
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).
- Centralized models generate deterministic condition-within-time,
  time-within-condition and optional difference-in-differences
  contrasts, with a common result contract across DESeq2, ALDEx2,
  ANCOM-BC2, corncob, MaAsLin3, Wilcoxon and LEfSe
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).

### Improvements

- Prepared recipes record executed and skipped methods, validate engine
  capabilities before execution and separate consensus results by
  contrast and effect direction
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).
- Step export and import preserve centralized model specifications
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).

## dar 1.9.1

- `Recipe` and `PrepRecipe` enforce structural invariants through S4
  validity checks
  ([\#124](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/124)).
- [`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
  [`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md)
  and preprocessing results are validated early, with contextual errors
  identifying the failing step
  ([\#124](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/124)).

## dar 1.5.6

### Breaking changes

- [`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  is defunct; its model implementation and related references are
  removed.

## dar 1.5.5

### Breaking changes

- [`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md)
  uses MaAsLin3 in place of MaAsLin2. Review method-specific arguments
  when migrating existing recipes.

## dar 1.5.4

- Parallel execution limits connection setup to support servers with
  more than 128 cores.
- Parallel execution uses a temporary future plan, preserving the user’s
  global future configuration.

## dar 1.5.3

- ANCOM result tables join taxonomy by `taxa_id`; corncob, LEfSe and
  Wilcoxon tables place `taxa_id` first, restoring consensus extraction
  through
  [`bake()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/bake.md)
  and
  [`cool()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/cool.md)
  ([\#117](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/117)).

## dar 1.1.3

- [`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md)
  is temporarily disabled.

## dar 1.1.2

- The data import vignette is corrected.

## dar 0.99.13

- The HUMAnN example is removed from the data import vignette.

## dar 0.99.10

- [`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md)
  is reintroduced after corncob returns to CRAN.

## dar 0.99.9

- [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  uses four workers to avoid issues with Bioconductor build
  infrastructure.

## dar 0.99.8

- Examples take less time to execute.

## dar 0.99.7

- Examples take less time to execute.

## dar 0.99.6

- Dependency installation guidance recommends
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).

## dar 0.99.5

- Vignettes take less time to execute.

## dar 0.99.4

- Tests and examples take less time to execute.
- The Linux GitHub Actions workflow correctly installs rlang.

## dar 0.99.3

### Breaking changes

- `Recipe` and `PrepRecipe` replace the former `recipe` and
  `prep_recipe` classes. Update code that inspects the old class names.
- [`required_deps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/required_deps.md)
  is internal and is no longer exported. Follow the installation
  guidance emitted by recipe steps for missing dependencies.
- R 4.4.0 or later is required.
- The package stops re-exporting `%>%` and `:=`. Use `|>` in recipes and
  qualify functions from their originating packages when needed.

### New features

- [`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  accepts both phyloseq and TreeSummarizedExperiment objects.
- Filtering and subsetting steps support chainable recipe construction.
- [`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
  [`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
  [`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md)
  and
  [`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md)
  provide dedicated filtering criteria.

### Documentation and improvements

- The import tutorial covers BIOM, QIIME, mothur, MetaPhlAn and HUMAnN
  conversions to phyloseq and TreeSummarizedExperiment.
- Installation guidance uses
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  and
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  instead of unconventional installation methods.
- The data.table dependency is removed.
- Rarefaction warns about its fixed seed and the implications for
  results.

### Bug fixes

- `run_aldex()` corrects its internal seed handling.
- The data documentation uses the correct dataset name.

## dar 0.99.0

- Initial Bioconductor submission.
