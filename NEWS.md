# dar 1.0.0

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
