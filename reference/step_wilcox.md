# Wilcox analysis

Performs a wilcox test to determine features (be it Operational
Taxanomic Unit (OTU), species, etc.) that are differentially abundant
between two or more groups of multiple samples.

## Usage

``` r
step_wilcox(
  rec,
  norm_method = "compositional",
  max_significance = 0.05,
  p_adj_method = "BH",
  rarefy = FALSE,
  id = rand_id("wilcox"),
  engine_args = list()
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- norm_method:

  Transformation to apply. The options include: 'compositional' (ie
  relative abundance), 'Z', 'log10', 'log10p', 'hellinger', 'identity',
  'clr', 'alr', or any method from the vegan::decostand function.

- max_significance:

  The q-value threshold for significance.

- p_adj_method:

  Character. Specifying the method to adjust p-values for multiple
  comparisons. Default is “BH” (Benjamini-Hochberg procedure).

- rarefy:

  Boolean indicating if OTU counts must be rarefyed.

- id:

  A character string that is unique to this step to identify it.

- engine_args:

  Named lists of advanced arguments for the native `test` stage.
  Arguments managed by dar or exposed above cannot be overridden.

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
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea"))

## Define step with default parameters
rec <- step_wilcox(rec)
#> Registered S3 method overwritten by 'car':
#>   method           from
#>   na.action.merMod lme4
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
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Bizcocho 
#> 
#> DA steps:
#> 
#>      ◉ step_wilcox() id = wilcox__Apple_pie 
```
