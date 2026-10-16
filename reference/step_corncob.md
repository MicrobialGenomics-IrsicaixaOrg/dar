# corncob analysis

Corncob an individual taxon regression model that uses abundance tables
and sample data. corncob is able to model differential abundance and
differential variability, and addresses each of the challenges presented
below:

## Usage

``` r
step_corncob(
  rec,
  phi.formula = stats::formula(~1),
  formula_null = NULL,
  phi.formula_null = NULL,
  link = "logit",
  phi.link = "logit",
  test = "Wald",
  boot = FALSE,
  B = 1000,
  filter_discriminant = TRUE,
  fdr_cutoff = 0.05,
  fdr = "fdr",
  log2FC = 0,
  rarefy = FALSE,
  id = rand_id("corncob"),
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

- phi.formula:

  An object of class formula without the response: a symbolic
  description of the model to be fitted to the dispersion.

- formula_null:

  Retained design argument, default `NULL`. Explicit values are
  overridden with a warning by the
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  contrasts.

- phi.formula_null:

  Formula for overdispersion under null, without response. Default
  `NULL`. Explicit values are overridden with a warning by the
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  contrasts.

- link:

  Link function for abundance covariates, defaults to "logit".

- phi.link:

  Link function for dispersion covariates, defaults to "logit".

- test:

  Hypothesis test, default `"Wald"`. Only Wald contrasts are supported;
  other tests make the step incompatible with the shared model.

- boot:

  Whether to request bootstrap inference, default `FALSE`. Bootstrap is
  unsupported and makes the step incompatible.

- B:

  Retained bootstrap iteration setting. Not used for Wald contrasts.

- filter_discriminant:

  Boolean. Defaults to TRUE. If FALSE, discriminant taxa will not be
  filtered out.

- fdr_cutoff:

  Integer. Defaults to 0.05. Desired type 1 error rate.

- fdr:

  Character. Defaults to "fdr". False discovery rate control method, see
  p.adjust for more options.

- log2FC:

  Minimum absolute effect on the reported coefficient scale. Despite the
  argument name, a logit-link coefficient is not a log2 fold change.
  Default zero.

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

  Named lists of advanced arguments for the native `fit` stage.
  Arguments managed by dar or exposed above cannot be overridden.

## Value

An object of class `Recipe`

## Details

- different sequencing depth

- excessive zeros from unobserved taxa

- high variability of empirical relative abundances (overdispersion)

- within-taxon correlation

- hypothesis testing with categorical and continuous covariates

Centralized models use non-bootstrap Wald contrasts. With
`engine_args = list(fit = list(robust = TRUE))`, inference uses
corncob's native sandwich covariance rather than the inverse Hessian.
This does not add support for random effects or repeated-measures
dependence.

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
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

rec_default <- step_corncob(rec)

# Use a stricter significance threshold.
rec_strict <- step_corncob(rec, fdr_cutoff = 0.01)
```
