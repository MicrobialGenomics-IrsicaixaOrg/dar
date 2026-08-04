# DESeq2 analysis

Differential expression analysis based on the Negative Binomial (a.k.a.
Gamma-Poisson) distribution. This function performs a default analysis
through the steps: 1) estimation of size factors: `estimateSizeFactors`.
2) estimation of dispersion: `estimateDispersions`. 3) Negative Binomial
GLM fitting and Wald statistics: nbinomWaldTest. For complete details on
each step, see the manual pages of the respective functions. After the
DESeq function returns a DESeqDataSet object, results tables (log2 fold
changes and p-values) can be generated using the results function.
Shrunken LFC can then be generated using the lfcShrink function.

## Usage

``` r
step_deseq(
  rec,
  test = "Wald",
  fitType = "local",
  betaPrior = FALSE,
  type = "ashr",
  max_significance = 0.05,
  log2FC = 0,
  rarefy = FALSE,
  id = rand_id("deseq"),
  engine_args = list()
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- test:

  Either "Wald" or "LRT", which will then use either Wald significance
  tests (defined by nbinomWaldTest), or the likelihood ratio test on the
  difference in deviance between a full and reduced model formula
  (defined by nbinomLRT).

- fitType:

  either "parametric", "local", "mean", or "glmGamPoi" for the type of
  fitting of dispersions to the mean intensity. See estimateDispersions
  for description.

- betaPrior:

  whether or not to put a zero-mean normal prior on the non-intercept
  coefficients See nbinomWaldTest for description of the calculation of
  the beta prior. In versions \>=1.16, the default is set to FALSE, and
  shrunken LFCs are obtained afterwards using lfcShrink.

- type:

  "apeglm" is the adaptive Student's t prior shrinkage estimator from
  the 'apeglm' package; "ashr" is the adaptive shrinkage estimator from
  the 'ashr' package, using a fitted mixture of normals prior - see the
  Stephens (2016) reference below for citation; "normal" is the 2014
  DESeq2 shrinkage estimator using a Normal prior.

- max_significance:

  The q-value threshold for significance.

- log2FC:

  log2FC cutoff.

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

  Named lists of advanced arguments for the native `size_factors`,
  `fit`, `results`, or `shrink` stage. Arguments managed by dar or
  exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
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
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.4 * length(x))")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preprocessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Djevrek 
#>      ◉ step_filter_taxa() id = filter_taxa__Spina_santa 
#> 
#> DA steps:
#> 

## Pass advanced arguments to native DESeq2 stages
step_deseq(
  rec,
  engine_args = list(
    fit = list(minReplicatesForReplace = Inf),
    shrink = list(lfcThreshold = 1)
  )
)
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preprocessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Djevrek 
#>      ◉ step_filter_taxa() id = filter_taxa__Spina_santa 
#> 
#> DA steps:
#> 
#>      ◉ step_deseq() id = deseq__Streusel 

## Define step with default parameters and prep
rec <-
  step_deseq(rec) |>
  prep(parallel = FALSE)
#> Warning: Estimated rdf < 1.0; not estimating variance
#> Warning: Estimated rdf < 1.0; not estimating variance
#> Warning: Estimated rdf < 1.0; not estimating variance

rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 76 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Results:
#> 
#>      ✔ deseq__Rugelach diff_taxa = 4 
#> 
#>      ℹ 5 taxon-contrast effects are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  step_deseq(rarefy = TRUE)

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
#>      ◉ step_deseq() id = deseq__Lotus_seed_bun 
```
