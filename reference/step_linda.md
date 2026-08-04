# LinDA differential-abundance analysis

`step_linda()` fits bias-corrected linear or linear mixed-effects models
with
[`MicrobiomeStat::linda()`](https://rdrr.io/pkg/MicrobiomeStat/man/linda.html).
It requires a centralized model defined by
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
and executes every contrast in that model's contrast plan.

## Usage

``` r
step_linda(
  rec,
  prev_filter = 0,
  mean_abund_filter = 0,
  max_abund_filter = 0,
  winsorize = TRUE,
  outlier_pct = 0.03,
  adaptive = TRUE,
  zero_handling = c("pseudo-count", "imputation"),
  pseudo_count = 0.5,
  corr_cut = 0.1,
  p_adj_method = "BH",
  alpha = 0.05,
  n_cpus = 1L,
  rarefy = FALSE,
  id = rand_id("linda"),
  engine_args = list()
)
```

## Arguments

- rec:

  A
  [Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
  object with a centralized model.

- prev_filter:

  Minimum feature prevalence retained for analysis.

- mean_abund_filter:

  Minimum mean relative abundance retained.

- max_abund_filter:

  Minimum maximum relative abundance retained.

- winsorize:

  Whether high-abundance outliers are winsorized.

- outlier_pct:

  Expected fraction of outliers used for winsorization.

- adaptive:

  Whether LinDA chooses zero handling from depth-covariate correlations.

- zero_handling:

  Zero treatment used when `adaptive = FALSE`: either `"pseudo-count"`
  or `"imputation"`.

- pseudo_count:

  Positive pseudo-count used by pseudo-count zero handling.

- corr_cut:

  Significance cutoff used by adaptive zero handling.

- p_adj_method:

  Multiple-testing correction accepted by
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- alpha:

  Adjusted p-value threshold used to classify significance.

- n_cpus:

  Number of cores used by LinDA mixed-effects models.

- rarefy:

  Whether counts are rarefied for this step.

- id:

  Unique identifier for this configured step.

- engine_args:

  Named lists of advanced arguments for the native `fit` stage.
  Arguments managed by dar or exposed above cannot be overridden.

## Value

A
[Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
object.

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy) |>
  add_model(
    ~ RiskGroup2,
    targets = "RiskGroup2",
    tax_level = "Species"
  ) |>
  step_linda(prev_filter = 0.1, adaptive = TRUE)

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
#> Preprocessing steps:
#> 
#> 
#> DA steps:
#> 
#>      ◉ step_linda() id = linda__Remonce 
```
