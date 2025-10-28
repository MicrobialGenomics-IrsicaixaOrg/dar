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
  id = rand_id("wilcox")
)

# S4 method for class 'Recipe'
step_wilcox(
  rec,
  norm_method = "compositional",
  max_significance = 0.05,
  p_adj_method = "BH",
  rarefy = FALSE,
  id = rand_id("wilcox")
)

# S4 method for class 'PrepRecipe'
step_wilcox(
  rec,
  norm_method = "compositional",
  max_significance = 0.05,
  p_adj_method = "BH",
  rarefy = FALSE,
  id = rand_id("wilcox")
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
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_metagenomeseq.md)

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
#>      ◉ step_subset_taxa() id = subset_taxa__Kołacz 
#>      ◉ step_filter_taxa() id = filter_taxa__Franzbrötchen 
#> 
#> DA steps:
#> 

## Define step with default parameters and prep
rec <-
  step_wilcox(rec) |>
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
#>      ✔ wilcox__Briouat diff_taxa = 5 
#> 
#>      ℹ 5 taxa are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_wilcox(rarefy = TRUE)
#> ! Run wilcox without rarefaction is not recommended (id = wilcox__Papanași)

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
#>      ◉ step_wilcox() id = wilcox__Papanași 
```
