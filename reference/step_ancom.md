# ANCOM analysis

Determine taxa whose absolute abundances, per unit volume, of the
ecosystem (e.g., gut) are significantly different with changes in the
covariate of interest (e.g., group). The current version of ancombc2
function implements Analysis of Compositions of Microbiomes with Bias
Correction (ANCOM-BC2) in cross-sectional and repeated measurements
data. In addition to the two-group comparison, ANCOM-BC2 also supports
testing for continuous covariates and multi-group comparisons, including
the global test, pairwise directional test, Dunnett's type of test, and
trend test.

## Usage

``` r
step_ancom(
  rec,
  fix_formula = recipe_targets(rec),
  rand_formula = NULL,
  p_adj_method = "holm",
  prv_cut = 0.1,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  rarefy = FALSE,
  id = rand_id("ancom"),
  engine_args = list()
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- fix_formula:

  the character string expresses how the microbial absolute abundances
  for each taxon depend on the fixed effects in metadata. When
  specifying the fix_formula, make sure to include the group variable in
  the formula if it is not NULL.

- rand_formula:

  the character string expresses how the microbial absolute abundances
  for each taxon depend on the random effects in metadata. ANCOM-BC2
  follows the lmerTest package in formulating the random effects. See
  ?lmerTest::lmer for more details. Default is NULL.

- p_adj_method:

  character. method to adjust p-values. Default is "holm". Options
  include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
  "none". See ?stats::p.adjust for more details.

- prv_cut:

  a numerical fraction between 0 and 1. Taxa with prevalences less than
  prv_cut will be excluded in the analysis. For instance, suppose there
  are 100 samples, if a taxon has nonzero counts presented in less than
  10 samples, it will not be further analyzed. Default is 0.10.

- lib_cut:

  a numerical threshold for filtering samples based on library sizes.
  Samples with library sizes less than lib_cut will be excluded in the
  analysis. Default is 0, i.e. do not discard any sample.

- s0_perc:

  a numerical fraction between 0 and 1. Inspired by Significance
  Analysis of Microarrays (SAM) methodology, a small positive constant
  is added to the denominator of ANCOM-BC2 test statistic corresponding
  to each taxon to avoid the significance due to extremely small
  standard errors, especially for rare taxa. This small positive
  constant is chosen as s0_perc-th percentile of standard error values
  for each fixed effect. Default is 0.05 (5th percentile).

- group:

  character. The name of the group variable in metadata. group should be
  discrete. Specifying group is required for detecting structural zeros
  and performing multi-group comparisons (global test, pairwise
  directional test, Dunnett's type of test, and trend test). Default is
  NULL. If the group of interest contains only two categories, leave it
  as NULL.

- struc_zero:

  logical. Whether to detect structural zeros based on group. Default is
  FALSE. See Details for a more comprehensive discussion on structural
  zeros.

- neg_lb:

  logical. Whether to classify a taxon as a structural zero using its
  asymptotic lower bound. Default is FALSE.

- alpha:

  numeric. Level of significance. Default is 0.05.

- n_cl:

  numeric. The number of nodes to be forked. For details, see
  ?parallel::makeCluster. Default is 1 (no parallel computing).

- verbose:

  logical. Whether to generate verbose output during the ANCOM-BC2
  fitting process. Default is FALSE.

- global:

  logical. Whether to perform the global test. Default is FALSE.

- pairwise:

  logical. Whether to perform the pairwise directional test. Default is
  FALSE.

- dunnet:

  logical. Whether to perform the Dunnett's type of test. Default is
  FALSE.

- trend:

  logical. Whether to perform trend test. Default is FALSE.

- rarefy:

  Boolean indicating if OTU counts must be rarefyed. This rarefaction
  uses the standard R sample function to resample from the abundance
  values in the otu_table component of the first argument, physeq. Often
  one of the major goals of this procedure is to achieve parity in total
  number of counts between samples, as an alternative to other formal
  normalization procedures, which is why a single value for the
  sample.size is expected. If 'no_seed', rarefaction is performed
  without a set seed.

- id:

  A character string that is unique to this step to identify it.

- engine_args:

  Named lists of advanced arguments for the native `fit` stage.
  Arguments managed by dar or exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.4 * length(x))")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Genus 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Strudel 
#>      ◉ step_filter_taxa() id = filter_taxa__Cornish_pasty 
#> 
#> DA steps:
#> 

## Pass advanced arguments to ANCOM-BC2
step_ancom(
  rec,
  engine_args = list(
    fit = list(iter_control = list(tol = 1e-3, max_iter = 50))
  )
)
#> Warning: ! The centralized model overrides design arguments in 1 step.
#> ℹ ancom__Öçpoçmaq: fix_formula
#> ℹ Method-specific thresholds and preprocessing controls are unchanged.
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Genus 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Strudel 
#>      ◉ step_filter_taxa() id = filter_taxa__Cornish_pasty 
#> 
#> DA steps:
#> 
#>      ◉ step_ancom() id = ancom__Öçpoçmaq 

## Define step with default parameters and prep
rec <-
  step_ancom(rec) |>
  prep(parallel = FALSE)
#> Warning: ! The centralized model overrides design arguments in 1 step.
#> ℹ ancom__Pionono: fix_formula
#> ℹ Method-specific thresholds and preprocessing controls are unchanged.
#> Warning: The number of taxa used for estimating sample-specific biases is: 41
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Loading required package: foreach
#> Loading required package: rngtools
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
#> Warning: The number of taxa used for estimating sample-specific biases is: 41
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
#> Warning: The number of taxa used for estimating sample-specific biases is: 41
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.

rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 76 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: factor, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Genus 
#> 
#> Results:
#> 
#>      ✔ ancom__Pionono diff_taxa = 13 
#> 
#>      ℹ 21 taxon-contrast effects are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  step_ancom(rarefy = TRUE)
#> Warning: ! The centralized model overrides design arguments in 1 step.
#> ℹ ancom__Tu: fix_formula
#> ℹ Method-specific thresholds and preprocessing controls are unchanged.

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preporcessing steps:
#> 
#> 
#> DA steps:
#> 
#>      ◉ step_ancom() id = ancom__Tu 
```
