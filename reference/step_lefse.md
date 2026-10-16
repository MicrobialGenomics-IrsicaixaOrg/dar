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
  id = rand_id("lefse"),
  engine_args = list()
)
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe. Execution requires
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md)
  before
  [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md);
  the step can be configured before adding the model.

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

  Unique, non-empty identifier for this step. DA step IDs must not use
  the reserved internal names `taxa_id`, `contrast_id`, `comparison`,
  `contrast_type`, `var`, `effect`, `sum`, `.tax_order`,
  `.contrast_order` or `.effect_order`. Other identifiers, including
  names such as `.custom`, are retained unchanged.

- engine_args:

  Named lists of advanced arguments for the native `transform` or `fit`
  stage. Arguments managed by dar or exposed above cannot be overridden.

## Value

An object of class `Recipe`

## See also

Other Diff taxa steps:
[`step_aldex2()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex2.md),
[`step_aldex3()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex3.md),
[`step_ancom()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_ancom.md),
[`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md),
[`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md),
[`step_linda()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_linda.md),
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md),
[`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)

## Examples

``` r
data(metaHIV_phy)

rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum")

rec_default <- step_lefse(rec)

# Require a larger LDA score.
rec_higher_lda <- step_lefse(rec, lda.threshold = 3)
```
