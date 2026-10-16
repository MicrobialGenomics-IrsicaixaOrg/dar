# ALDEx3 fixed- and mixed-effects analysis

Fit planned contrasts with the official ALDEx3 engine. Fixed models use
native linear regression; mixed models use exact `lme4` fits, not the
approximate `blmm` engine. A centralized
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
is required.

## Usage

``` r
step_aldex3(
  rec,
  max_significance = 0.05,
  nsample = 2000,
  scale_model = "clr",
  scale_args = list(gamma = 0.5),
  p_adj_method = "BH",
  n_cores = 1,
  seed = 130L,
  rarefy = FALSE,
  id = rand_id("aldex3"),
  engine_args = list()
)
```

## Arguments

- rec:

  An unprepared Recipe object.

- max_significance:

  Native adjusted P value cutoff.

- nsample:

  Number of Monte Carlo draws passed to ALDEx3. Default 2000.

- scale_model:

  One of `"clr"`, `"tss"`, `"sample"`, `"coefficient"`.

- scale_args:

  Named list of scale-model parameters; see Details.

- p_adj_method:

  Native multiple-testing adjustment, default `"BH"`.

- n_cores:

  Number of cores for native mixed-model fitting, default one.

- seed:

  Integer Monte Carlo seed. Restores the caller's RNG state.

- rarefy:

  Whether to rarefy before fitting; default `FALSE`.

- id:

  Unique step identifier.

- engine_args:

  Named list with an optional `fit` list. Supports native `test`
  (`"t.HC3"`, `"t.HC0"`, `"t"`) for fixed models and `streamsize`.
  Validation is left to ALDEx3. The native test setting does not apply
  to mixed models. Data, model, scale, engine, seed and two-sided
  inference are managed by dar.

## Value

A Recipe with the ALDEx3 step added.

## Details

Supply sequence counts for inference. Taxonomic aggregation and common
preprocessing precede native Dirichlet sampling. No extra prevalence
filter or implicit "other" feature is added. Filtering changes the
composition and the features defining CLR scale assumptions. External
scale information must describe the retained community; removed counts
are not reconstructed. Native numerical parameters and input values are
passed to ALDEx3. Engine-specific validation and errors are left to that
package. dar checks its configuration and taxon, sample and contrast
identities.

Scale models are specified by name, not executable functions. `"clr"`
and `"tss"` accept `scale_args = list(gamma = 0.5)`, the native default
SD on the log2 coefficient scale. `gamma = 0` fixes the corresponding
scale assumption. `"sample"` requires named `s.mu` and `s.sd` vectors of
external log2-scale means and SDs, one per input sample ID.
`"coefficient"` requires named `c.mu` and `c.sd` vectors, one per column
of the resolved fixed-effects design. This interface supports
independent SDs, not covariance matrices or arbitrary scale functions.
The default is a modeling choice, not an estimate of scale uncertainty.
Coefficient coding and the units of numeric covariates affect these
priors; choose scale information appropriate to the study.

Scale uncertainty and priors always refer to the original resolved
design, even when a planned contrast requires releveling. Each contrast
uses the same seed and scale-model basis. The native model is refitted
for each contrast; mixed models can therefore be expensive.

Effects are mean Monte Carlo coefficients on the log2 abundance scale
conditional on the specified scale model
(`aldex3_log2_abundance_coefficient`). Sequence counts alone do not
measure absolute microbial load. P values and adjusted P values retain
ALDEx3's native directional Monte Carlo procedure, including adjustment
across taxa separately for each coefficient. They are not BH applied to
averaged P values, nor a joint adjustment over all planned contrasts.
Missing statistics remain missing. Execution settings and the effective
package version are recorded by
[`get_execution()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_execution.md).

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
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
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus") |>
  step_aldex3(scale_model = "clr", scale_args = list(gamma = 0.5))
# Use verified sequence counts for inference; prep() executes the recipe.
```
