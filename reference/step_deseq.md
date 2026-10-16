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
  for this Recipe. Execution requires
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  before
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md);
  the step can be configured before adding the model.

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

  Unique, non-empty identifier for this step. DA step IDs must not use
  the reserved internal names `taxa_id`, `contrast_id`, `comparison`,
  `contrast_type`, `var`, `effect`, `sum`, `.tax_order`,
  `.contrast_order` or `.effect_order`. Other identifiers, including
  names such as `.custom`, are retained unchanged.

- engine_args:

  Named lists of advanced arguments for the native `size_factors`,
  `fit`, `results`, or `shrink` stage. Arguments managed by dar or
  exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_rademu()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rademu.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")

rec_default <- step_deseq(rec)

# Pass advanced arguments to native DESeq2 stages.
rec_custom <- step_deseq(
  rec,
  engine_args = list(
    fit = list(minReplicatesForReplace = Inf),
    shrink = list(lfcThreshold = 1)
  )
)
```
