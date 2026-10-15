# ALDEx2 analysis

Differential abundance analysis using ALDEx2's Monte Carlo log-ratio
transformation and native statistical tests. With a centralized model,
fits a GLM for each planned contrast and reports its coefficient and
adjusted P value.

## Usage

``` r
step_aldex(
  rec,
  max_significance = 0.05,
  mc.samples = 128,
  denom = "all",
  rarefy = FALSE,
  id = rand_id("aldex"),
  engine_args = list()
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- max_significance:

  Adjusted P value cutoff for the GLM, or the Benjamini-Hochberg
  corrected Welch P value in model-free execution.

- mc.samples:

  The number of Monte Carlo instances to use to estimate the underlying
  distributions. The default is 128.

- denom:

  An any variable (all, iqlr, zero, lvha, median, user) indicating
  features to use as the denominator for the Geometric Mean calculation
  The default "all" uses the geometric mean abundance of all features.
  Using "median" returns the median abundance of all features. Using
  "iqlr" uses the features that are between the first and third quartile
  of the variance of the clr values across all samples. Using "zero"
  uses the non-zero features in each grop as the denominator. This
  approach is an extreme case where there are many nonzero features in
  one condition but many zeros in another. Using "lvha" uses features
  that have low variance (bottom quartile) and high relative abundance
  (top quartile in every sample). It is also possible to supply a vector
  of row indices to use as the denominator. Here, the experimentalist is
  determining a-priori which rows are thought to be invariant. In the
  case of RNA-seq, this could include ribosomal protein genes and and
  other house-keeping genes. This should be used with caution because
  the offsets may be different in the original data and in the data used
  by the function because features that are 0 in all samples are removed
  by aldex.clr. Centralized models require `"all"` or explicit numeric
  denominator indices.

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

  Named lists of advanced arguments for the native `clr` and `fit`
  stages. The former `effect` stage is not supported: centralized
  effects now come from the GLM. Arguments managed by dar or exposed
  above cannot be overridden.

## Value

An object of class `Recipe`

## Details

With
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md),
`effect` is the contrast-oriented coefficient from
[`ALDEx2::aldex.glm()`](https://rdrr.io/pkg/ALDEx2/man/aldex.glm.html),
averaged over Monte Carlo instances, and `effect_metric` is
`"aldex_glm_coefficient"`. It comes from the same GLM as the reported
adjusted P value. Under the default denominator, scale and Gaussian
identity-link settings, the response is CLR in base 2; this coefficient
is not an absolute-abundance log2 fold change or the standardized effect
from
[`ALDEx2::aldex.effect()`](https://rdrr.io/pkg/ALDEx2/man/aldex.effect.html).
Advanced native settings can change the response or link scale.
Recompute previously prepared centralized ALDEx2 results to obtain these
coefficients.

Deprecated model-free execution retains native two-group Welch tests and
standardized effects from
[`ALDEx2::aldex.effect()`](https://rdrr.io/pkg/ALDEx2/man/aldex.effect.html),
identified as `"standardized_clr_effect"` in tidy results.

## See also

Other Diff taxa steps:
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")

rec_default <- step_aldex(rec)

# Increase the number of Monte Carlo instances.
rec_more_instances <- step_aldex(rec, mc.samples = 256)
```
