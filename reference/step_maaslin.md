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
  for this Recipe.

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

  A character string that is unique to this step to identify it.

- engine_args:

  Named lists of advanced arguments for the native `fit` or `contrast`
  stage. Arguments managed by dar or exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
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
#>      ◉ step_subset_taxa() id = subset_taxa__Crocetta_of_Caltanissetta 
#>      ◉ step_filter_taxa() id = filter_taxa__Makroudh 
#> 
#> DA steps:
#> 

## Define step with default parameters and prep
rec <-
  step_maaslin(rec) |>
  prep(parallel = FALSE)

rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 76 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: factor, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Results:
#> 
#>      ✔ maaslin__Strudel diff_taxa = 1 
#> 
#>      ℹ 1 taxon-contrast effects are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  step_maaslin(rarefy = TRUE)

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
#>      ◉ step_maaslin() id = maaslin__Punschkrapfen 
```
