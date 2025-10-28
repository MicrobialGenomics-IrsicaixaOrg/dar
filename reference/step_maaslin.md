# MaAsLin2 analysis

MaAsLin2 finds associations between microbiome meta-omics features and
complex metadata in population-scale epidemiological studies. The
software includes multiple analysis methods (including support for
multiple covariates and repeated measures), filtering, normalization,
and transform options to customize analysis for your specific study.

## Usage

``` r
step_maaslin(
  rec,
  min_abundance = 0,
  min_prevalence = 0.1,
  min_variance = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.25,
  random_effects = NULL,
  correction = "BH",
  standardize = TRUE,
  reference = NULL,
  rarefy = FALSE,
  id = rand_id("maaslin")
)

# S4 method for class 'Recipe'
step_maaslin(
  rec,
  min_abundance = 0,
  min_prevalence = 0.1,
  min_variance = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.25,
  random_effects = NULL,
  correction = "BH",
  standardize = TRUE,
  reference = NULL,
  rarefy = FALSE,
  id = rand_id("maaslin")
)

# S4 method for class 'PrepRecipe'
step_maaslin(
  rec,
  min_abundance = 0,
  min_prevalence = 0.1,
  min_variance = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.25,
  random_effects = NULL,
  correction = "BH",
  standardize = TRUE,
  reference = NULL,
  rarefy = FALSE,
  id = rand_id("maaslin")
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
  "CLR", "CSS", "NONE", "TMM".

- transform:

  The transform to apply. Default: "LOG". Choices: "LOG", "LOGIT",
  "AST", "NONE".

- analysis_method:

  The analysis method to apply. Default: "LM". Choices: "LM", "CPLM",
  "ZICP", "NEGBIN", "ZINB".

- max_significance:

  The q-value threshold for significance.

- random_effects:

  The random effects for the model, comma-delimited for multiple
  effects.

- correction:

  The correction method for computing the q-value.

- standardize:

  Apply z-score so continuous metadata are on the same scale.

- reference:

  The factor to use as a reference for a variable with more than two
  levels provided as a string of 'variable,reference' semi-colon
  delimited for multiple variables.

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

## See also

Other Diff taxa steps:
[`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md),
[`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_metagenomeseq.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
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
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Puff_pastry 
#>      ◉ step_filter_taxa() id = filter_taxa__Kūčiukai 
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
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Results:
#> 
#>      ✔ maaslin__Paxlava diff_taxa = 65 
#> 
#>      ℹ 65 taxa are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_maaslin(rarefy = TRUE)

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Preporcessing steps:
#> 
#> 
#> DA steps:
#> 
#>      ◉ step_maaslin() id = maaslin__Öçpoçmaq 
```
