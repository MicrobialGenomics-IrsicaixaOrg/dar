# MaAsLin3 analysis

MaAsLin 3 finds associations between microbiome meta-omics features and
complex metadata. It uses a unified framework to test for both abundance
(using linear models) and prevalence (using logistic regression).

## Usage

``` r
step_maaslin(
  rec,
  min_abundance = 0,
  min_prevalence = 0.1,
  min_variance = 0,
  normalization = "TSS",
  transform = "LOG",
  max_significance = 0.1,
  random_effects = NULL,
  correction = "BH",
  standardize = TRUE,
  reference = NULL,
  median_comparison_abundance = TRUE,
  rarefy = FALSE,
  id = rand_id("maaslin"),
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

- min_abundance:

  The minimum abundance for each feature.

- min_prevalence:

  The minimum percent of samples for which a feature is detected at
  minimum abundance.

- min_variance:

  Keep features with variance greater than.

- normalization:

  The normalization method to apply. Default: "TSS". Choices: "TSS",
  "CLR", "NONE".

- transform:

  The transform to apply. Default: "LOG" (Base 2). Choices: "LOG",
  "PLOG", "NONE".

- max_significance:

  The q-value threshold for significance.

- random_effects:

  The random effects for the model (vector of character strings).

- correction:

  The correction method for computing the q-value.

- standardize:

  Apply z-score so continuous metadata are on the same scale.

- reference:

  The factor to use as a reference for a variable with more than two
  levels provided as a string of 'variable,reference' semi-colon
  delimited.

- median_comparison_abundance:

  Test abundance coefficients against a null value corresponding to the
  median coefficient for a metadata variable across the features.
  Recommended for relative abundance (default: TRUE).

- rarefy:

  Boolean indicating if OTU counts must be rarefyed.

- id:

  Unique, non-empty identifier for this step. DA step IDs must not use
  the reserved internal names `taxa_id`, `contrast_id`, `comparison`,
  `contrast_type`, `var`, `effect`, `sum`, `.tax_order`,
  `.contrast_order` or `.effect_order`. Other identifiers, including
  names such as `.custom`, are retained unchanged.

- engine_args:

  Named lists of advanced arguments for the native `fit` or `contrast`
  stage. Arguments managed by dar or exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
[`step_rademu()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rademu.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")

rec_default <- step_maaslin(rec)

# Increase the minimum feature prevalence.
rec_filtered <- step_maaslin(rec, min_prevalence = 0.2)
```
