# dar 1.9.12

## Breaking changes

* R 4.6 or later is required.
* `step_metagenomeseq()` is removed after its defunct release cycle. Replace
  it with a supported DA step when rebuilding legacy recipes (#157).
* DA step IDs must not use reserved consensus column names such as `effect`
  or `sum`. Rename conflicting IDs, for example to `da_effect`, before
  constructing or importing a recipe. Other IDs, including those beginning
  with a dot, remain unchanged.
* ComplexHeatmap, gplots, heatmaply and UpSetR are optional dependencies.
  Install the relevant package explicitly to use its plotting functions;
  missing dependencies produce installation guidance (#157).
* Centralized ALDEx2 results report `aldex_glm_coefficient` rather than a
  standardized two-group effect. Recompute previously prepared results to
  obtain the corrected effects; their former effect values are not relabeled.
  The obsolete `engine_args$effect` stage is rejected. Deprecated model-free
  ALDEx2 execution retains its standardized effect.

## Deprecations

* Legacy selectors and model-free DA execution remain available in
  Bioconductor 3.24, become defunct in 3.25 and are scheduled for removal in
  3.26. Use `add_model()` to configure DA analyses and `get_model()` to
  inspect their targets and taxonomic rank (#157).

## Bug fixes

* `step_aldex()` extracts centralized effects from the same native GLM and
  contrast as its adjusted P values, preserving adjustment and interaction
  semantics instead of combining them with a post-hoc binary-group effect.
* Modeled engines, consensus tables and abundance plots share one aggregated
  taxonomic view with stable IDs. Abundance transformations use aggregated
  counts, and public accessors retain the complete underlying data.
* `abundance_plt()`, `intersection_plt()` and `mutual_plt()` check selectors
  against planned model hypotheses, including those with no significant taxa.
  For empty selections, discovery plots report a message and return invisible
  `NULL`. Abundance plots with explicit `taxa_ids` remain available regardless
  of significance.
* `add_model()` override warnings are emitted only for explicitly supplied
  engine-specific designs; default ANCOM-BC2 and corncob steps remain quiet.
* `prep()` emits the DA rarefaction caution when rarefying steps execute,
  rather than when they are constructed.
* `step_ancom()` preserves stable taxon IDs after ANCOM-BC2 aggregation.
* `step_ancom()` uses ANCOM-BC2's robust significance decision when
  pseudocount sensitivity analysis is enabled, in both centralized and
  deprecated model-free execution. Explicitly disabling sensitivity analysis
  retains the basic native decision.
* `step_corncob()` honors the requested absolute effect threshold in
  deprecated model-free execution.
* `step_corncob()` uses native sandwich covariance for centralized Wald
  contrasts when `robust = TRUE` and forwards the requested `B` during
  deprecated model-free bootstrap execution.
* `step_linda()` locally corrects MicrobiomeStat 1.4's adaptive zero-handling
  capitalization error, restoring depth-dependent imputation without
  modifying the installed dependency or the remaining native algorithm.
* `step_filter_taxa()` executes allowed `median()`, `var()` and `sd()`
  predicates locally and after JSON import, using fixed trusted bindings
  without expanding the permitted syntax or capturing external functions.
* `step_lefse()` joins modeled results by taxon ID, keeping identically named
  taxa from different lineages distinct.
* `step_wilcox()` and `step_lefse()` skip unsupported additive time
  adjustment, including when the time column is explicitly declared.

## Documentation

* Function examples illustrate supported workflows without redundant setup,
  repeated defaults or intentional error demonstrations.
* `step_filter_taxa()` documents its single, restricted serializable
  predicate and filtering during `prep()`, rather than arbitrary function
  lists or a `prune` argument.
* Six focused vignettes cover getting started, data inputs, preprocessing,
  models and contrasts, results and optional consensus, and recipe sharing.
  The duplicate case study is merged into the results guide. Examples use
  the central model API, execute adjusted and longitudinal designs, and
  import a bundled BIOM file into both supported input classes.
* Results guides start with per-method findings and distinguish agreement
  from calibrated error control and comparable effect scales. Input guidance
  distinguishes relative abundances from raw counts; the advanced ANCOM-BC2
  example supplies a complete native control list.

## Development infrastructure

* The development-only Monte Carlo validation suite and its dedicated
  workflow are archived in the frozen
  `codex/scientific-validation-frozen-20261015` branch and removed from
  `devel`. The historical baseline is preserved without recalibration or
  relabeling. Previously generated reports are retained separately.
* Statistical correctness and native-engine regression tests remain in
  `tests/testthat`, alongside the ordinary package checks and vignette builds.

# dar 1.9.11

## Bug fixes

* Consensus tables, plots and bake results keep modeled targets, contrasts
  and effect directions as separate hypotheses (#154).
* Modeled intersections use the complete taxon-by-contrast universe,
  including non-significant hypotheses, and classify zero effects as
  neutral (#154).
* Method exclusions are applied before consensus support is counted.
  Invalid cutoffs, weights and exclusions produce classified errors (#154).
* Modeled exploration functions accept explicit `target` and `contrast_id`
  selectors; abundance plots require them when automatic selection is
  ambiguous (#154).

# dar 1.9.10

## Breaking changes

* `otu_table()` always returns counts with taxa in rows and deterministic
  taxon and sample ordering. Adapt code that assumes the stored phyloseq
  orientation (#153).
* `sample_data()` and `tax_table()` expose complete recipe metadata and
  taxonomy instead of model-restricted views. Use `get_model()` to inspect
  the analysis configuration, and select columns explicitly when needed
  (#153).

## Improvements

* DA engines and plots use explicit internal target-only and selected-rank
  views, keeping analysis inputs independent of public accessor semantics
  (#153).

# dar 1.9.9

## Breaking changes

* `export_steps()` and `import_steps()` use a versioned, non-executable
  JSON schema. Legacy pseudo-JSON files are rejected. Recreate trusted legacy
  recipes with the current API and export them again (#152).
* `import_steps()` requires explicit `prepare = TRUE` to execute an
  analysis when importing bake configurations into an unprepared recipe
  (#152).

## Improvements

* Imported constructors, arguments, formulas, nested values and filter
  predicates are allowlisted and validated transactionally before the
  recipe is changed (#152).

# dar 1.9.8

* Recipe steps are classified by their concrete classes, so custom IDs do not
  change preprocessing, DA or bake semantics (#151).
* Preprocessing operations execute in their configured order, including
  recipes that interleave filtering, subsetting and rarefaction (#151).
* Unsupported, misplaced and duplicate-ID steps fail with classified
  validation errors before execution (#151).

# dar 1.9.7

## Breaking changes

* `phy_qc()` is defunct. Use `recipe_qc()` with an explicit `group_by`
  when grouped quality-control metrics are required (#112).

## New features

* `recipe_qc()` provides model-independent overall or explicitly grouped
  quality-control metrics for `Recipe` and `PrepRecipe` objects.
  Summaries expose unambiguous `group_by` and `group` columns and handle
  either count-table orientation (#112).

# dar 1.9.6

* `step_linda()` integrates bias-corrected LinDA linear and mixed-effects
  models with every compatible centralized contrast. LinDA participates in
  the tidy result contract and Monte Carlo validation with explicit log2
  fold-change semantics (#128).

# dar 1.9.5

* The seven centralized DA engines accept structured `engine_args` for
  advanced arguments at explicit native-call stages (#126).
* Step export and import preserve nested advanced arguments and reject
  values that cannot be reconstructed faithfully (#126).

# dar 1.9.4

* `tidy_results()` provides a stable, engine-labelled result contract while
  preserving the native tables in `PrepRecipe@results` (#127).
* Result consumers validate taxon-contrast-step keys and use explicit effect
  metrics and p-value mappings (#127).

# dar 1.9.3

## Breaking changes

* `add_model()` becomes the single source of truth for DA targets and
  taxonomic resolution. Supply `targets` and `tax_level` there (#142).
* `add_var()`, `add_tax()`, `get_var()`, `get_tax()`, the
  `var_info` and `tax_info` arguments to `recipe()`, and model-free
  DA execution are deprecated. They remain functional during the first
  Bioconductor deprecation cycle; migrate to `add_model()` and
  `get_model()` (#142).

## Improvements

* Legacy centralized models, serialized recipes and imported step files are
  normalized to the new model representation without changing original
  microbiome metadata (#142).

## Bug fixes

* Modeled exclusion and mutual plots preserve contrast and effect-direction
  keys, and normalize adjusted p-value aliases without duplicate columns
  (#142).
* LEfSe executes every planned pairwise contrast for multilevel targets
  (#142).

# dar 1.9.2

## New features

* `add_model()` and `get_model()` define and inspect a shared statistical
  design with fixed confounders, target-time interactions, supported random
  effects, explicit references and a common missing-value policy (#130).
* Centralized models generate deterministic condition-within-time,
  time-within-condition and optional difference-in-differences contrasts,
  with a common result contract across DESeq2, ALDEx2, ANCOM-BC2, corncob,
  MaAsLin3, Wilcoxon and LEfSe (#130).

## Improvements

* Prepared recipes record executed and skipped methods, validate engine
  capabilities before execution and separate consensus results by contrast
  and effect direction (#130).
* Step export and import preserve centralized model specifications (#130).

# dar 1.9.1

* `Recipe` and `PrepRecipe` enforce structural invariants through S4
  validity checks (#124).
* `add_var()`, `add_tax()` and preprocessing results are validated
  early, with contextual errors identifying the failing step (#124).

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
