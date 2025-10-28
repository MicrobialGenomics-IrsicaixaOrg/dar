# lefse analysis

Lefser is metagenomic biomarker discovery tool that is based on LEfSe
tool and is published by Huttenhower et al. 2011. Lefser is the R
implementation of the LEfSe method. Using statistical analyses, lefser
compares microbial populations of healthy and diseased subjects to
discover differencially expressed microorganisms. Lefser than computes
effect size, which estimates magnitude of differential expression
between the populations for each differentially expressed microorganism.
Subclasses of classes can also be assigned and used within the analysis.

## Usage

``` r
step_lefse(
  rec,
  kruskal.threshold = 0.05,
  wilcox.threshold = 0.05,
  lda.threshold = 2,
  subclassCol = NULL,
  assay = 1L,
  trim.names = FALSE,
  rarefy = TRUE,
  id = rand_id("lefse")
)

# S4 method for class 'Recipe'
step_lefse(
  rec,
  kruskal.threshold = 0.05,
  wilcox.threshold = 0.05,
  lda.threshold = 2,
  subclassCol = NULL,
  assay = 1L,
  trim.names = FALSE,
  rarefy = TRUE,
  id = rand_id("lefse")
)

# S4 method for class 'PrepRecipe'
step_lefse(
  rec,
  kruskal.threshold = 0.05,
  wilcox.threshold = 0.05,
  lda.threshold = 2,
  subclassCol = NULL,
  assay = 1L,
  trim.names = FALSE,
  rarefy = TRUE,
  id = rand_id("lefse")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- kruskal.threshold:

  numeric(1) The p-value for the Kruskal-Wallis Rank Sum Test (default
  0.05).

- wilcox.threshold:

  numeric(1) The p-value for the Wilcoxon Rank-Sum Test when 'blockCol'
  is present (default 0.05).

- lda.threshold:

  numeric(1) The effect size threshold (default 2.0).

- subclassCol:

  character(1) Optional column name in 'colData(expr)' indicating the
  blocks, usually a factor with two levels (e.g., 'c("adult",
  "senior")'; default NULL).

- assay:

  The i-th assay matrix in the ‘SummarizedExperiment' (’expr'; default
  1).

- trim.names:

  If 'TRUE' extracts the most specific taxonomic rank of organism.

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
#>      ◉ step_subset_taxa() id = subset_taxa__Canelé 
#>      ◉ step_filter_taxa() id = filter_taxa__Suncake 
#> 
#> DA steps:
#> 

## Define step with default parameters
rec <- step_lefse(rec) 
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
#>      ◉ step_subset_taxa() id = subset_taxa__Canelé 
#>      ◉ step_filter_taxa() id = filter_taxa__Suncake 
#> 
#> DA steps:
#> 
#>      ◉ step_lefse() id = lefse__Cuban_pastelito 

## Running lefse without rarefaction (not recommended)
rec <- 
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_lefse(rarefy = FALSE)
#> ! Run lefse without rarefaction is not recommended (id = lefse__Dutch_Baby_Pancake)
  
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
#>      ◉ step_lefse() id = lefse__Dutch_Baby_Pancake 
```
