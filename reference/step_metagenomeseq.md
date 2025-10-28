# MetagenomeSeq analysis

metagenomeSeq is designed to determine features (be it Operational
Taxanomic Unit (OTU), species, etc.) that are differentially abundant
between two or more groups of multiple samples. metagenomeSeq is
designed to address the effects of both normalization and under-sampling
of microbial communities on disease association detection and the
testing of feature correlations.

## Usage

``` r
step_metagenomeseq(
  rec,
  zeroMod = NULL,
  useCSSoffset = TRUE,
  useMixedModel = FALSE,
  max_significance = 0.05,
  log2FC = 0,
  rarefy = FALSE,
  rm_zeros = 0,
  id = rand_id("metagenomeseq")
)

# S4 method for class 'Recipe'
step_metagenomeseq(
  rec,
  zeroMod = NULL,
  useCSSoffset = TRUE,
  useMixedModel = FALSE,
  max_significance = 0.05,
  log2FC = 0,
  rarefy = FALSE,
  rm_zeros = 0,
  id = rand_id("metagenomeseq")
)

# S4 method for class 'PrepRecipe'
step_metagenomeseq(
  rec,
  zeroMod = NULL,
  useCSSoffset = TRUE,
  useMixedModel = FALSE,
  max_significance = 0.05,
  log2FC = 0,
  rarefy = FALSE,
  rm_zeros = 0,
  id = rand_id("metagenomeseq")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- zeroMod:

  The zero model, the model to account for the change in the number of
  OTUs observed as a linear effect of the depth of coverage.

- useCSSoffset:

  Boolean, whether to include the default scaling parameters in the
  model or not.

- useMixedModel:

  Estimate the correlation between duplicate features or replicates
  using duplicateCorrelation.

- max_significance:

  The q-value threshold for significance.

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

- rm_zeros:

  Proportion of samples of the same categorical level with more than 0
  counts.

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
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

## Init Recipe
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.02 * length(x))")

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
#>      ◉ step_subset_taxa() id = subset_taxa__Roti_john 
#>      ◉ step_filter_taxa() id = filter_taxa__Roze_koek 
#> 
#> DA steps:
#> 

## Define step with default parameters and prep
rec <-
  step_metagenomeseq(rec, rm_zeros = 0.01) |>
  prep(parallel = FALSE)
#> Warning: NaNs produced

rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 291 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Phylum 
#> 
#> Results:
#> 
#>      ✔ metagenomeseq__Coussin_de_Lyon diff_taxa = 8 
#> 
#>      ℹ 8 taxa are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_metagenomeseq(rarefy = TRUE)

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
#>      ◉ step_metagenomeseq() id = metagenomeseq__Malsouka 
```
