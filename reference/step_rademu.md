# radEmu differential-abundance analysis

Fit planned contrasts with
[`radEmu::emuFit()`](https://statdivlab.github.io/radEmu/reference/emuFit.html).
A centralized
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
defines the fixed design and taxonomic analysis units. The default
robust score test is recommended by radEmu; Wald testing is available
explicitly.

## Usage

``` r
step_rademu(
  rec,
  test = c("score", "wald"),
  cluster = NULL,
  max_significance = 0.05,
  p_adj_method = "BH",
  constraint_param = 0.1,
  seed = 130L,
  rarefy = FALSE,
  id = rand_id("rademu"),
  engine_args = list()
)
```

## Arguments

- rec:

  An unprepared
  [Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
  object.

- test:

  Native test to report: `"score"` (default) or `"wald"`.

- cluster:

  Optional single metadata column identifying dependent observations.
  `NULL` uses native independent-observation inference.

- max_significance:

  Adjusted P value cutoff.

- p_adj_method:

  Multiple-testing adjustment passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- constraint_param:

  Native pseudo-Huber smoothing parameter, default 0.1.

- seed:

  Integer seed used for native fitting. Restores the caller's RNG.

- rarefy:

  Whether to rarefy before fitting; default `FALSE`.

- id:

  Unique step identifier.

- engine_args:

  Named list with an optional `fit` list of native arguments, for
  example `tolerance`, `maxit`, `maxit_null`, `alpha` and
  `remove_zero_comparison_pvals`. Solver controls can be passed directly
  in `fit`. Data, design, clustering, centring functions, fitted-model
  reuse, tested coordinates, test selection and result shape are managed
  by dar; `control` and `return_both_score_pvals` cannot be overridden.

## Value

A Recipe with the radEmu step added.

## Details

Effects are natural-log fold differences relative to the native
pseudo-Huber centre across analysis taxa. They are not unqualified
absolute-abundance changes or log2 fold changes. Each planned contrast
is reparameterized and fitted with its own native centring constraint.
This retains the native hypothesis rather than treating the nonlinear
centre as a linear transformation of previously centred coefficients.

Supply `cluster` to identify a metadata column containing cluster
labels. Labels are aligned with the common model cohort. Clustered
inference uses radEmu's GEE framework, not a mixed model. Random-effect
formulas such as `(1 | subject)` are incompatible and are never silently
translated. For repeated observations, use a fixed formula and
explicitly set `cluster = "subject"`. Missing cluster IDs are not
dropped privately.

dar adjusts the selected native P values across returned taxa separately
for each planned contrast. It does not adjust across all contrasts
jointly. Missing native statistics remain missing, including native
zero-comparison exclusions. Native convergence warnings are retained;
diagnostics and effective settings are available through
[`get_execution()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_execution.md).

The complete taxonomic analysis matrix is used in every fit. Counts are
not split by taxon or restricted to the two groups named in a contrast.
Score tests can be expensive because each taxon requires a restricted
fit. Native numerical validation and errors are left to radEmu.

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
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
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus") |>
  step_rademu()
# prep(rec) executes native score tests for every planned contrast and taxon.
```
