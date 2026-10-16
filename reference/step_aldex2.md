# ALDEx2 analysis

Differential abundance analysis using ALDEx2's Monte Carlo log-ratio
transformation and native statistical tests. Requires
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md),
fits a GLM for each planned contrast and reports its coefficient and
adjusted P value.

## Usage

``` r
step_aldex2(
  rec,
  max_significance = 0.05,
  mc.samples = 128,
  denom = "all",
  rarefy = FALSE,
  id = rand_id("aldex2"),
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

- max_significance:

  Adjusted P value cutoff for the GLM.

- mc.samples:

  The number of Monte Carlo instances to use to estimate the underlying
  distributions. The default is 128.

- denom:

  Features used as the geometric-mean denominator: `"all"` (default) or
  explicit numeric row indices in the aggregated count matrix. Choose
  indices a priori and account for native removal of all-zero taxa.
  Data-dependent denominators are not supported with the shared design.

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
  stages. Effects come from the GLM. Arguments managed by dar or exposed
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

## See also

Other Diff taxa steps:
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
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

rec_default <- step_aldex2(rec)

# Increase the number of Monte Carlo instances.
rec_more_instances <- step_aldex2(rec, mc.samples = 256)
```
