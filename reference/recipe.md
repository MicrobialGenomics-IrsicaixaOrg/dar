# Create a Recipe for preprocessing data

A Recipe is a description of the steps to be applied to a data set in
order to prepare it for data analysis.

## Usage

``` r
recipe(
  microbiome_object = NULL,
  var_info = NULL,
  tax_info = NULL,
  steps = list()
)
```

## Arguments

- microbiome_object:

  Phyloseq-class object or TreeSummarizedExperiment-class object.

- var_info:

  A character string of column names corresponding to variables that
  will be used in any context.

- tax_info:

  A character string of taxonomic levels that will be used in any
  context.

- steps:

  list with steps.

## Value

An object of class `Recipe` with sub-objects:

- phyloseq:

  object of class `phyloseq` with taxa abundance information.

- var_info:

  A tibble that contains the current set of terms in the data set. This
  initially defaults to the same data contained in `var_info`.

- tax_info:

  A tibble that contains the current set of taxonomic levels that will
  be used in the analysis.

## Examples

``` r
data(metaHIV_phy)

## Define recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Phylum") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") |>
  step_metagenomeseq(rm_zeros = 0.01) |>
  step_maaslin()

## Prep recipe
da_results <- prep(rec)

## Consensus strategy
n_methods <- 2
da_results <- bake(da_results, count_cutoff = n_methods)

## Results
cool(da_results)
#> ℹ Bake for count_cutoff = 2
#> # A tibble: 1 × 2
#>   taxa_id taxa         
#>   <chr>   <chr>        
#> 1 Otu_96  Bacteroidetes

## You can also crate a recipe without var and tax info
rec <- recipe(metaHIV_phy)

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined variable of interest. Use add_var() to add it to Recipe! 
#>      ✖ undefined taxonomic level. Use add_tax() to add it to Recipe! 
#> 

## And define them later
rec <- rec |>
  add_var("RiskGroup2") |>
  add_tax("Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Genus 
#> 
#> 

## When trying to add an identical step to an existing one, the system
## returns an information message.
rec <- step_aldex(rec)
rec <- step_aldex(rec)
#> ! This step is already defined with the same parameters and will be skipped: 
#> rec %>% step_aldex(max_significance = 0.05, mc.samples = 128, denom = c('all'), rarefy = FALSE)

## The same with bake
da_results <- bake(da_results)
da_results <- bake(da_results)
#> ! This step is already defined with the same parameters and will be skipped: 
#> rec %>% step_bake(count_cutoff = NULL, weights = NULL, exclude = NULL)
```
