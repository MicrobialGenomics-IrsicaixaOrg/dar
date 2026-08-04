# Changelog

## dar 1.9.12

### SCIENTIFIC VALIDATION

- IMPROVED: Full Monte Carlo validation now enforces calibrated LinDA
  bias and coverage gates, compares results with a reviewed stochastic
  baseline, records dependency versions and publishes concise execution
  summaries
  ([\#155](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/155)).
- IMPROVED: Full validation reports are retained as compact
  machine-readable history on the `validation-results` branch, while
  weak longitudinal contrasts remain descriptive rather than receiving
  inappropriate power gates
  ([\#155](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/155)).

### LIFECYCLE

- REMOVED:
  [`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  after a complete defunct release cycle;
  [`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  remains defunct with migration guidance
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).
- DEPRECATED: Legacy selectors and model-free differential-abundance
  execution remain available in Bioconductor 3.24, become defunct in
  3.25 and are scheduled for removal in 3.26
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).
- CHANGED: Require R 4.6 and make ComplexHeatmap, gplots, heatmaply and
  UpSetR optional plotting dependencies with classified installation
  guidance
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).
- IMPROVED: Remove generated test artifacts and address actionable
  BiocCheck findings in documentation, vignettes and package metadata
  ([\#157](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/157)).

## dar 1.9.11

### FIXED

- Consensus tables, plots and bake results now keep modeled targets,
  contrasts and effect directions as separate hypotheses instead of
  mixing their support counts
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Modeled intersections use the complete taxon-by-contrast universe,
  including non-significant hypotheses, and classify zero effects as
  neutral
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Method exclusions are applied before consensus support is counted, and
  invalid cutoffs, weights and exclusions now produce classified errors
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).
- Modeled exploration functions accept explicit `target` and
  `contrast_id` selectors; abundance plots require them whenever
  automatic selection would otherwise be ambiguous
  ([\#154](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/154)).

## dar 1.9.10

### BREAKING CHANGES

- CHANGED:
  [`sample_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)
  and
  [`tax_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)
  now return the complete recipe metadata and taxonomy instead of
  model-restricted analysis views
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).
- CHANGED:
  [`otu_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)
  now always returns a wide taxa-by-sample tibble with deterministic
  taxon and sample ordering, regardless of the stored phyloseq
  orientation
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).
- IMPROVED: Differential-abundance engines and plots use explicit
  internal target-only and selected-rank views, preventing public
  accessor semantics from changing analysis inputs
  ([\#153](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/153)).

## dar 1.9.9

### BREAKING CHANGES

- CHANGED:
  [`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
  and
  [`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
  now use a real, versioned, non-executable JSON schema. Files produced
  by the legacy pseudo-JSON format are rejected with migration guidance
  instead of being evaluated
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).
- CHANGED: Imported bake configurations require an explicit
  `prepare = TRUE` before
  [`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
  may execute an analysis on an unprepared recipe
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).
- IMPROVED: Imported constructors, arguments, formulas, nested values
  and filter predicates are allowlisted and validated transactionally
  before a recipe is changed
  ([\#152](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/152)).

## dar 1.9.8

- FIXED: Recipe steps are classified by their concrete classes, so
  arbitrary user-defined IDs no longer change preprocessing,
  differential-abundance or bake semantics
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).
- FIXED: Preprocessing operations execute in their configured order,
  including recipes that interleave filtering, subsetting and
  rarefaction
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).
- IMPROVED: Unsupported, misplaced and duplicate-ID recipe steps now
  fail with classified validation errors before execution
  ([\#151](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/151)).

## dar 1.9.7

### BREAKING CHANGES

- DEFUNCT:
  [`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-defunct.md)
  now directs users to
  [`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
  instead of inferring an analysis target implicitly
  ([\#112](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/112)).
- NEW:
  [`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
  provides model-independent overall or explicitly grouped
  quality-control metrics for both `Recipe` and `PrepRecipe` objects
  ([\#112](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/112)).
- IMPROVED: QC summaries expose unambiguous `group_by` and `group`
  columns and reshape correctly oriented count data only once per call
  ([\#112](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/112)).

## dar 1.9.6

- NEW:
  [`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md)
  integrates bias-corrected LinDA linear and mixed-effects models with
  every compatible centralized contrast
  ([\#128](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/128)).
- IMPROVED: LinDA participates in the canonical tidy result contract and
  the Monte Carlo scientific-validation workflow with explicit log2
  fold-change semantics
  ([\#128](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/128)).

## dar 1.9.5

- NEW: The seven centralized DA engines accept structured `engine_args`
  for advanced arguments at explicit native-call stages
  ([\#126](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/126)).
- IMPROVED: Step export/import preserves nested advanced arguments and
  rejects values that cannot be reconstructed faithfully
  ([\#126](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/126)).

## dar 1.9.4

- NEW:
  [`tidy_results()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_results.md)
  exposes a stable, engine-labelled result contract while preserving the
  engine-specific tables stored in `PrepRecipe@results`
  ([\#127](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/127)).
- IMPROVED: Result consumers validate taxon-contrast-step keys and use
  explicit effect metrics and p-value mappings instead of
  method-specific aliases
  ([\#127](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/127)).

## dar 1.9.3

### BREAKING CHANGES

- CHANGED:
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  is now the single source of truth for differential abundance analysis
  targets and taxonomic resolution through its new `targets` and
  `tax_level` arguments
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).
- DEPRECATED:
  [`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
  [`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md),
  [`get_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_var.md),
  [`get_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_tax.md),
  the `var_info` and `tax_info` arguments to
  [`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md),
  and model-free DA execution now emit classified migration warnings.
  They remain functional during this first Bioconductor deprecation
  cycle
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).
- IMPROVED: Legacy centralized models, serialized recipes and imported
  step files are normalized to the new model representation without
  changing the original microbiome metadata
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).
- FIXED: Modeled exclusion and mutual plots preserve contrast and
  effect-direction keys, normalize adjusted p-value aliases without
  duplicate columns, and LEfSe executes every planned pairwise contrast
  for multilevel targets
  ([\#142](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/142)).

## dar 1.9.2

- NEW: Add a centralized statistical model with
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  and
  [`get_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_model.md),
  including fixed confounders, target-time interactions, supported
  random effects, explicit references, and a shared missing-value policy
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).
- NEW: Generate deterministic condition-within-time,
  time-within-condition and optional difference-in-differences
  contrasts, with a common result contract across DESeq2, ALDEx2,
  ANCOM-BC2, corncob, MaAsLin3, Wilcoxon and LEfSe
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).
- IMPROVED: Record executed and skipped methods in prepared recipes,
  validate engine capabilities before execution, and keep consensus
  results separated by contrast and effect direction
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).
- IMPROVED: Preserve centralized model specifications when exporting and
  importing recipe steps
  ([\#130](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/130)).

## dar 1.9.1

- IMPROVED: Enforce structural invariants for `Recipe` and inherited
  `PrepRecipe` objects through S4 validity checks
  ([\#124](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/124)).
- IMPROVED: Validate
  [`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md),
  [`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md),
  and preprocessing results early, with contextual errors that identify
  the failing step
  ([\#124](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/124)).

## dar 1.5.6

- CHANGED: Mark step_metagenomeseq as defunct and remove its
  implementation and related references.

## dar 1.5.5

- CHANGED: Migration from Maaslin2 to maaslin3

## dar 1.5.4

- FIXED: Resolved a crash occurring on servers with \>128 cores by
  limiting connection setup based on available R connections.
- CHANGED: Refactored parallel execution to use
  `with(..., local = TRUE)`. This ensures the internal parallel plan is
  temporary and does not overwrite the user’s global future
  configuration.

## dar 1.5.3

### Bug Fixes

- <https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/117>

## dar 1.1.3

- `setp_ancom` temporally deactivation

## dar 1.1.2

### Bug Fixes

- `data_import.Rmd`

## dar 0.99.13

### Bug Fixes

- Remove humann example form data import vignette

## dar 0.99.10

### New Features

- Reimplementing `step_corncob` after the return of corncob package to
  cran

## dar 0.99.9

### Bug Fixes

- Set workers parameter to 4 in order to avoid issues with BBS builds

## dar 0.99.8

### Improvements

- Reducing examples computation time

## dar 0.99.7

### Improvements

- Reducing examples computation time

## dar 0.99.6

### Improvements

- Suggest using
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  to install `dar` dependencies

## dar 0.99.5

### Improvements

- Reducing vignettes computation time

## dar 0.99.4

### Improvements

- Reducing tests and examples computation time

### Bug Fixes

- Fixing bug in Github Actions on Linux with `rlang` installation.

## dar 0.99.3

### New Features

- The `dar` package now accepts both `phyloseq` class objects and
  `TreeSummarizedExperiment` as inputs.
- The tutorial has been refocused to become a tutorial on how to import
  biom, qiime, mothur, metaphlan, and humann into
  `TreeSummarizedExperiment` and `phyloseq` class objects.
- The `Recipe` and `PrepRecipe` classes have been introduced, replacing
  the previous `recipe` and `prep_recipe` classes.
- The `subset` and `filter` operations have been updated to allow all
  steps of the `recipe` to be defined in a chainable manner.
- The functions `step_filter_by_abundance`, `step_filter_by_prevalence`,
  `step_filter_by_rarity`, and `step_filter_by_variance` have been added
  to enhance filtering functionality.

### Improvements

- The R version dependency has been updated to 4.4.0.
- The dependency on `data.table` has been removed.
- The re-export of `%>%` and `:=` has been removed. Now code examples
  and vignettes use `|>`.
- The `required_deps` function is no longer exported.
- The package now recommends more commonly used installation methods,
  such as
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  or
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).
- A warning message is now displayed whenever the `rarefy = TRUE` option
  is used, informing users that a fixed seed is being used and how it
  could impact their results.
- The package coverage has increased to 82.33%.

### Bug Fixes

- Unconventional package installation methods have been avoided, for
  example,
  [`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html).
- The setting of a seed within a function (`run_aldex`) has been
  addressed.
- The name of the data in the `R/data.R` documentation has been
  corrected from `NA`.

## dar 0.99.0

- Initial Bioconductor submission.
