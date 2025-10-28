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
  formula_null = stats::formula(~1),
  phi.formula_null = stats::formula(~1),
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
  id = rand_id("corncob")
)

# S4 method for class 'Recipe'
step_corncob(
  rec,
  phi.formula = stats::formula(~1),
  formula_null = stats::formula(~1),
  phi.formula_null = stats::formula(~1),
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
  id = rand_id("corncob")
)

# S4 method for class 'PrepRecipe'
step_corncob(
  rec,
  phi.formula = stats::formula(~1),
  formula_null = stats::formula(~1),
  phi.formula_null = stats::formula(~1),
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
  id = rand_id("corncob")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- phi.formula:

  An object of class formula without the response: a symbolic
  description of the model to be fitted to the dispersion.

- formula_null:

  Formula for mean under null, without response.

- phi.formula_null:

  Formula for overdispersion under null, without response.

- link:

  Link function for abundance covariates, defaults to "logit".

- phi.link:

  Link function for dispersion covariates, defaults to "logit".

- test:

  Character. Hypothesis testing procedure to use. One of "Wald" or "LRT"
  (likelihood ratio test).

- boot:

  Boolean. Defaults to FALSE. Indicator of whether or not to use
  parametric bootstrap algorithm. (See pbWald and pbLRT).

- B:

  Optional integer. Number of bootstrap iterations. Ignored if boot is
  FALSE. Otherwise, defaults to 1000.

- filter_discriminant:

  Boolean. Defaults to TRUE. If FALSE, discriminant taxa will not be
  filtered out.

- fdr_cutoff:

  Integer. Defaults to 0.05. Desired type 1 error rate.

- fdr:

  Character. Defaults to "fdr". False discovery rate control method, see
  p.adjust for more options.

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

## Value

An object of class `Recipe`

## Details

- different sequencing depth

- excessive zeros from unobserved taxa

- high variability of empirical relative abundances (overdispersion)

- within-taxon correlation

- hypothesis testing with categorical and continuous covariates

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_metagenomeseq.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <- 
  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Ensaïmada 
#>      ◉ step_filter_taxa() id = filter_taxa__Mandelkubb 
#> 
#> DA steps:
#> 

## Define step with default parameters and prep
rec <- 
  step_corncob(rec) |>
  prep(parallel = FALSE)
  
rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 101 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Results:
#> 
#>      ✔ corncob__Pain_au_chocolat diff_taxa = 4 
#> 
#>      ℹ 4 taxa are present in all tested methods 
#> 
```
