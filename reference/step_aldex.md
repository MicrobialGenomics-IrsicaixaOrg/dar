# ALDEx2 analysis

A differential abundance analysis for the comparison of two or more
conditions. For example, single-organism and meta-RNA-seq
high-throughput sequencing assays, or of selected and unselected values
from in-vitro sequence selections. Uses a Dirichlet-multinomial model to
infer abundance from counts, that has been optimized for three or more
experimental replicates. Infers sampling variation and calculates the
expected false discovery rate given the biological and sampling
variation using the Wilcox rank test or Welches t-test (aldex.ttest) or
the glm and Kruskal Wallis tests (aldex.glm). Reports both P and fdr
values calculated by the Benjamini Hochberg correction (Not supported in
`dar` package).

## Usage

``` r
step_aldex(
  rec,
  max_significance = 0.05,
  mc.samples = 128,
  denom = "all",
  rarefy = FALSE,
  id = rand_id("aldex")
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- max_significance:

  Benjamini-Hochberg corrected P value of Welch’s t test cutoff.

- mc.samples:

  The number of Monte Carlo instances to use to estimate the underlying
  distributions; since we are estimating central tendencies, 128 is
  usually sufficient, but larger numbers may be .

- denom:

  An any variable (all, iqlr, zero, lvha, median, user) indicating
  features to use as the denominator for the Geometric Mean calculation
  The default "all" uses the geometric mean abundance of all features.
  Using "median" returns the median abundance of all features. Using
  "iqlr" uses the features that are between the first and third quartile
  of the variance of the clr values across all samples. Using "zero"
  uses the non-zero features in each grop as the denominator. This
  approach is an extreme case where there are many nonzero features in
  one condition but many zeros in another. Using "lvha" uses features
  that have low variance (bottom quartile) and high relative abundance
  (top quartile in every sample). It is also possible to supply a vector
  of row indices to use as the denominator. Here, the experimentalist is
  determining a-priori which rows are thought to be invariant. In the
  case of RNA-seq, this could include ribosomal protein genes and and
  other house-keeping genes. This should be used with caution because
  the offsets may be different in the original data and in the data used
  by the function because features that are 0 in all samples are removed
  by aldex.clr.

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

The run_aldex function is a wrapper that performs log-ratio
transformation and statistical testing in a single line of code.
Specifically, this function: (a) generates Monte Carlo samples of the
Dirichlet distribution for each sample, (b) converts each instance using
a log-ratio transform, then (c) returns test results for two sample
(Welch's t, Wilcoxon) test. This function also estimates effect size for
two sample analyses.

## See also

Other Diff taxa steps:
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
#>      ◉ step_subset_taxa() id = subset_taxa__Cannoli_siciliani 
#>      ◉ step_filter_taxa() id = filter_taxa__Fig_roll 
#> 
#> DA steps:
#> 

## Define ALDEX step with default parameters and prep
rec <-
  step_aldex(rec) |>
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
#>      ✔ aldex__Samosa diff_taxa = 3 
#> 
#>      ℹ 3 taxa are present in all tested methods 
#> 

## Wearing rarefaction only for this step
rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_aldex(rarefy = TRUE)

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
#>      ◉ step_aldex() id = aldex__Linzer_torte 
```
