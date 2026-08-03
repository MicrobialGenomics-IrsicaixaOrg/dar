# dar 1.9.5

* NEW: The seven centralized DA engines accept structured `engine_args` for
  advanced arguments at explicit native-call stages (#126).
* IMPROVED: Step export/import preserves nested advanced arguments and rejects
  values that cannot be reconstructed faithfully (#126).

# dar 1.9.4

* NEW: `tidy_results()` exposes a stable, engine-labelled result contract while
  preserving the engine-specific tables stored in `PrepRecipe@results` (#127).
* IMPROVED: Result consumers validate taxon-contrast-step keys and use explicit
  effect metrics and p-value mappings instead of method-specific aliases (#127).

# dar 1.9.3

## BREAKING CHANGES

* CHANGED: `add_model()` is now the single source of truth for differential
  abundance analysis targets and taxonomic resolution through its new
  `targets` and `tax_level` arguments (#142).
* DEPRECATED: `add_var()`, `add_tax()`, `get_var()`, `get_tax()`, the
  `var_info` and `tax_info` arguments to `recipe()`, and model-free DA
  execution now emit classified migration warnings. They remain functional
  during this first Bioconductor deprecation cycle (#142).
* IMPROVED: Legacy centralized models, serialized recipes and imported step
  files are normalized to the new model representation without changing the
  original microbiome metadata (#142).
* FIXED: Modeled exclusion and mutual plots preserve contrast and
  effect-direction keys, normalize adjusted p-value aliases without duplicate
  columns, and LEfSe executes every planned pairwise contrast for multilevel
  targets (#142).

# dar 1.9.2

* NEW: Add a centralized statistical model with `add_model()` and
  `get_model()`, including fixed confounders, target-time interactions,
  supported random effects, explicit references, and a shared missing-value
  policy (#130).
* NEW: Generate deterministic condition-within-time, time-within-condition and
  optional difference-in-differences contrasts, with a common result contract
  across DESeq2, ALDEx2, ANCOM-BC2, corncob, MaAsLin3, Wilcoxon and LEfSe
  (#130).
* IMPROVED: Record executed and skipped methods in prepared recipes, validate
  engine capabilities before execution, and keep consensus results separated
  by contrast and effect direction (#130).
* IMPROVED: Preserve centralized model specifications when exporting and
  importing recipe steps (#130).

# dar 1.9.1

* IMPROVED: Enforce structural invariants for `Recipe` and inherited
  `PrepRecipe` objects through S4 validity checks (#124).
* IMPROVED: Validate `add_var()`, `add_tax()`, and preprocessing results early,
  with contextual errors that identify the failing step (#124).

# dar 1.5.6

* CHANGED: Mark step_metagenomeseq as defunct and remove its implementation and related references.

# dar 1.5.5

* CHANGED: Migration from Maaslin2 to maaslin3

# dar 1.5.4

* FIXED: Resolved a crash occurring on servers with >128 cores by limiting connection setup based on available R connections.
* CHANGED: Refactored parallel execution to use `with(..., local = TRUE)`. This ensures the internal parallel plan is temporary and does not overwrite the user's global future configuration.

# dar 1.5.3

## Bug Fixes

- https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/117

# dar 1.1.3
- `setp_ancom` temporally deactivation

# dar 1.1.2

## Bug Fixes
- `data_import.Rmd`

# dar 0.99.13

## Bug Fixes
- Remove humann example form data import vignette

# dar 0.99.10

## New Features
- Reimplementing `step_corncob` after the return of corncob package to cran

# dar 0.99.9

## Bug Fixes
- Set workers parameter to 4 in order to avoid issues with BBS builds

# dar 0.99.8

## Improvements
- Reducing examples computation time

# dar 0.99.7

## Improvements
- Reducing examples computation time

# dar 0.99.6

## Improvements
- Suggest using `BiocManager::install()` to install `dar` dependencies

# dar 0.99.5

## Improvements
- Reducing vignettes computation time

# dar 0.99.4

## Improvements
- Reducing tests and examples computation time

## Bug Fixes
- Fixing bug in Github Actions on Linux with `rlang` installation. 

# dar 0.99.3

## New Features

- The `dar` package now accepts both `phyloseq` class objects and `TreeSummarizedExperiment` as inputs.
- The tutorial has been refocused to become a tutorial on how to import biom, qiime, mothur, metaphlan, and humann into `TreeSummarizedExperiment` and `phyloseq` class objects.
- The `Recipe` and `PrepRecipe` classes have been introduced, replacing the previous `recipe` and `prep_recipe` classes.
- The `subset` and `filter` operations have been updated to allow all steps of the `recipe` to be defined in a chainable manner.
- The functions `step_filter_by_abundance`, `step_filter_by_prevalence`, `step_filter_by_rarity`, and `step_filter_by_variance` have been added to enhance filtering functionality.

## Improvements

- The R version dependency has been updated to 4.4.0.
- The dependency on `data.table` has been removed.
- The re-export of `%>%` and `:=` has been removed. Now code examples and vignettes use `|>`.
- The `required_deps` function is no longer exported.
- The package now recommends more commonly used installation methods, such as `BiocManager::install()` or `install.packages()`.
- A warning message is now displayed whenever the `rarefy = TRUE` option is used, informing users that a fixed seed is being used and how it could impact their results.
- The package coverage has increased to 82.33%.

## Bug Fixes

- Unconventional package installation methods have been avoided, for example, `pak::pkg_install`.
- The setting of a seed within a function (`run_aldex`) has been addressed.
- The name of the data in the `R/data.R` documentation has been corrected from `NA`.

# dar 0.99.0

* Initial Bioconductor submission.
