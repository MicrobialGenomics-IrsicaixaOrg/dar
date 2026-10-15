# Data inputs and inspection

[`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
accepts a `phyloseq` or `TreeSummarizedExperiment` object. These
containers connect an abundance matrix to sample metadata and, when
available, taxonomy. Import files with the tools supplied by `phyloseq`
or `mia`, then inspect the resulting object before configuring an
analysis.

Importing data does not establish that their units, experimental design
or taxonomic annotations are appropriate for a particular
differential-abundance method.

## Start with an existing phyloseq object

The package includes `metaHIV_phy`, with its recorded sample metadata
and taxonomy. No statistical model is needed to inspect the data or
define preprocessing steps.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dar`](https://github.com/MicrobialGenomics-IrsicaixaOrg/dar)`)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"metaHIV_phy"``, package ``=`` ``"dar"``)`\
\
`rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``metaHIV_phy``)`\
`rec`\
`#> ── ``DAR Recipe`` ──────────────────────────────────────────────────────────────────`\
`#> Inputs:`\
`#> `\
`#>      ``ℹ`` phyloseq object with ``451`` taxa and ``156`` samples `\
`#>      ``✖`` undefined analysis target. Use ``add_model()`` to define the analysis! `\
`#>      ``✖`` undefined taxonomic level. Use ``add_model()`` to define the analysis!`

[`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
also accepts a `phyloseq` object whose stored abundance matrix has
samples in rows: the accessors below return one consistent orientation.
When constructing your own `phyloseq` object, set `taxa_are_rows` to
match the actual matrix, rather than changing this flag without
transposing the data.

## Inspect samples, taxonomy and abundances

Use the dar accessors explicitly when other attached packages export
functions with the same names.

\
`samples`` ``<-`` ``dar``::`[`sample_data`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)`(``rec``)`\
`taxonomy`` ``<-`` ``dar``::`[`tax_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)`(``rec``)`\
`abundances`` ``<-`` ``dar``::`[`otu_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)`(``rec``)`\
\
[`head`](https://rdrr.io/r/utils/head.html)`(``samples``[`[`c`](https://rdrr.io/r/base/c.html)`(``"sample_id"``, ``"RiskGroup2"``)``]``)`\
`#> ``# A tibble: 6 × 2`\
`#>   sample_id  RiskGroup2`\
`#>   ``<chr>``      ``<chr>``     `\
`#> ``1`` Sample_186 msm       `\
`#> ``2`` Sample_185 msm       `\
`#> ``3`` Sample_184 hts       `\
`#> ``4`` Sample_182 pwid      `\
`#> ``5`` Sample_181 msm       `\
`#> ``6`` Sample_170 msm`\
[`head`](https://rdrr.io/r/utils/head.html)`(``taxonomy``)`\
`#> ``# A tibble: 6 × 8`\
`#>   taxa_id Kingdom  Phylum         Class           Order     Family Genus Species`\
`#>   ``<chr>``   ``<chr>``    ``<chr>``          ``<chr>``           ``<chr>``     ``<chr>``  ``<chr>`` ``<chr>``  `\
`#> ``1`` Otu_1   Archaea  Euryarchaeota  Methanobacteria Methanob… Metha… Meth… Methan…`\
`#> ``2`` Otu_2   Archaea  Euryarchaeota  Methanobacteria Methanob… Metha… Meth… Methan…`\
`#> ``3`` Otu_3   Bacteria Actinobacteria Actinobacteria  Actinomy… Actin… Acti… Actino…`\
`#> ``4`` Otu_4   Bacteria Actinobacteria Actinobacteria  Actinomy… Actin… Acti… Actino…`\
`#> ``5`` Otu_5   Bacteria Actinobacteria Actinobacteria  Actinomy… Actin… Acti… Actino…`\
`#> ``6`` Otu_6   Bacteria Actinobacteria Actinobacteria  Actinomy… Actin… Acti… Actino…`\
`abundances``[`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``6``)``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``5``)``]`\
`#> ``# A tibble: 6 × 5`\
`#>   taxa_id Sample_186 Sample_185 Sample_184 Sample_182`\
`#>   ``<chr>``        ``<dbl>``      ``<dbl>``      ``<dbl>``      ``<dbl>`\
`#> ``1`` Otu_1         ``2``538          0       ``9``026      ``12``049`\
`#> ``2`` Otu_2            0          0          0        609`\
`#> ``3`` Otu_3            0          0          0          0`\
`#> ``4`` Otu_4            0          0          0          0`\
`#> ``5`` Otu_5            0          0          0          0`\
`#> ``6`` Otu_6            0          0          0          0`

The accessors expose the complete stored components, independently of
the target or taxonomic rank later selected with
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md):

- [`sample_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)
  returns `sample_id` followed by all sample metadata columns, in
  [`phyloseq::sample_names()`](https://rdrr.io/pkg/phyloseq/man/sample_names-methods.html)
  order.
- [`tax_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)
  returns `taxa_id` followed by all taxonomic ranks, in
  [`phyloseq::taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html)
  order.
- [`otu_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)
  returns one row per taxon, with `taxa_id` first and one column per
  sample. Taxon and sample order match the other two accessors.

\
[`stopifnot`](https://rdrr.io/r/base/stopifnot.html)`(`\
`  `[`identical`](https://rdrr.io/r/base/identical.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``abundances``)``[``-``1``]``, ``samples``$``sample_id``)``,`\
`  `[`identical`](https://rdrr.io/r/base/identical.html)`(``abundances``$``taxa_id``, ``taxonomy``$``taxa_id``)`\
`)`

Keep sample and taxon identifiers unique and consistent across input
files. Taxonomic labels are annotations, not replacements for feature
identifiers: several features can share a species name. Inspect missing
annotations and use standard rank names such as `Kingdom`, `Genus` and
`Species`.

If metadata or taxonomy is absent, the corresponding accessor returns an
ID-only tibble. That is sufficient for inspection, but it does not
supply the real covariates or taxonomic rank needed by a model. Add
genuine study metadata before calling
[`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md);
do not manufacture groups just to make an example run.

## Import a BIOM table

This example uses a file distributed with `phyloseq`, without
downloading data. The importer reads the sample metadata already present
in that file.

\
`biom_file`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(`\
`  ``"extdata"``, ``"rich_dense_otu_table.biom"``, package ``=`` ``"phyloseq"`\
`)`\
\
`phy`` ``<-`` ``phyloseq``::`[`import_biom`](https://rdrr.io/pkg/phyloseq/man/import_biom.html)`(`\
`  ``biom_file``,`\
`  parseFunction ``=`` ``phyloseq``::`[`parse_taxonomy_greengenes`](https://rdrr.io/pkg/phyloseq/man/parseTaxonomy-functions.html)\
`)`\
`biom_rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``phy``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``dar``::`[`sample_data`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)`(``biom_rec``)``)`\
`#> ``# A tibble: 6 × 5`\
`#>   sample_id BarcodeSequence LinkerPrimerSequence  BODY_SITE Description`\
`#>   ``<chr>``     ``<chr>``           ``<chr>``                 ``<chr>``     ``<chr>``      `\
`#> ``1`` Sample1   CGCTTATCGAGA    CATGCTGCCTCCCGTAGGAGT gut       human gut  `\
`#> ``2`` Sample2   CATACCAGTAGC    CATGCTGCCTCCCGTAGGAGT gut       human gut  `\
`#> ``3`` Sample3   CTCTCTACCTGT    CATGCTGCCTCCCGTAGGAGT gut       human gut  `\
`#> ``4`` Sample4   CTCTCGGCCTGT    CATGCTGCCTCCCGTAGGAGT skin      human skin `\
`#> ``5`` Sample5   CTCTCTACCAAT    CATGCTGCCTCCCGTAGGAGT skin      human skin `\
`#> ``6`` Sample6   CTAACTACCAAT    CATGCTGCCTCCCGTAGGAGT skin      human skin`\
[`head`](https://rdrr.io/r/utils/head.html)`(``dar``::`[`tax_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)`(``biom_rec``)``)`\
`#> ``# A tibble: 5 × 8`\
`#>   taxa_id  Kingdom  Phylum         Class              Order Family Genus Species`\
`#>   ``<chr>``    ``<chr>``    ``<chr>``          ``<chr>``              ``<chr>`` ``<chr>``  ``<chr>`` ``<chr>``  `\
`#> ``1`` GG_OTU_1 Bacteria Proteobacteria Gammaproteobacter… Ente… Enter… Esch… ``NA``     `\
`#> ``2`` GG_OTU_2 Bacteria Cyanobacteria  Nostocophycideae   Nost… Nosto… Doli… ``NA``     `\
`#> ``3`` GG_OTU_3 Archaea  Euryarchaeota  Methanomicrobia    Meth… Metha… Meth… ``NA``     `\
`#> ``4`` GG_OTU_4 Bacteria Firmicutes     Clostridia         Hala… Halan… Hala… Halana…`\
`#> ``5`` GG_OTU_5 Bacteria Proteobacteria Gammaproteobacter… Ente… Enter… Esch… ``NA`

The same file can be imported as a `TreeSummarizedExperiment`. Rows of a
TSE represent features, columns represent samples, `rowData()` holds
feature annotations and `colData()` holds sample metadata.

\
`tse`` ``<-`` ``mia``::`[`importBIOM`](https://microbiome.github.io/mia/reference/importBIOM.html)`(``biom_file``)`\
\
`# This bundled example contains these seven ranks, in this order.`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``SummarizedExperiment``::`[`rowData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``tse``)``)`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``"Kingdom"``, ``"Phylum"``, ``"Class"``, ``"Order"``, ``"Family"``, ``"Genus"``, ``"Species"`\
`)`\
\
`tse_rec`` ``<-`` `[`recipe`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)`(``tse``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``dar``::`[`tax_table`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)`(``tse_rec``)``)`\
`#> ``# A tibble: 5 × 8`\
`#>   taxa_id  Kingdom     Phylum            Class        Order Family Genus Species`\
`#>   ``<chr>``    ``<chr>``       ``<chr>``             ``<chr>``        ``<chr>`` ``<chr>``  ``<chr>`` ``<chr>``  `\
`#> ``1`` GG_OTU_1 k__Bacteria p__Proteobacteria c__Gammapro… o__E… f__En… g__E… s__    `\
`#> ``2`` GG_OTU_2 k__Bacteria p__Cyanobacteria  c__Nostocop… o__N… f__No… g__D… s__    `\
`#> ``3`` GG_OTU_3 k__Archaea  p__Euryarchaeota  c__Methanom… o__M… f__Me… g__M… s__    `\
`#> ``4`` GG_OTU_4 k__Bacteria p__Firmicutes     c__Clostrid… o__H… f__Ha… g__H… s__Hal…`\
`#> ``5`` GG_OTU_5 k__Bacteria p__Proteobacteria c__Gammapro… o__E… f__En… g__E… s__`

The rank relabeling above is specific to this file. Inspect the
annotations in your own data instead of assigning rank names by position
without checking their meaning. A TSE input is converted internally to
`phyloseq`; if the intended abundance table is not the default assay,
select it explicitly with
`mia::convertToPhyloseq(..., assay.type = "your_assay")` before creating
the recipe.

For other formats, consult the importers’ documentation rather than
treating all formats as equivalent:

- [`phyloseq::import_qiime()`](https://rdrr.io/pkg/phyloseq/man/import_qiime.html)
  reads legacy QIIME tables;
  [`mia::importQIIME2()`](https://microbiome.github.io/mia/reference/importQIIME2.html)
  reads QIIME 2 artifacts.
- [`phyloseq::import_mothur()`](https://rdrr.io/pkg/phyloseq/man/import_mothur.html)
  and
  [`mia::importMothur()`](https://microbiome.github.io/mia/reference/importMothur.html)
  read the supported mothur files.
- [`mia::importMetaPhlAn()`](https://microbiome.github.io/mia/reference/importMetaPhlAn.html)
  reads MetaPhlAn abundance profiles.

For example, open the installed help with
[`help("importQIIME2", package = "mia")`](https://microbiome.github.io/mia/reference/importQIIME2.html)
or
[`help("import_biom", package = "phyloseq")`](https://rdrr.io/pkg/phyloseq/man/import_biom.html).

## Check the measurement scale

An abundance matrix can contain raw read counts, proportions,
percentages or other transformed values. These are not interchangeable
inputs. In particular, MetaPhlAn relative-abundance profiles are not raw
read counts. Neither importing them nor converting their container makes
them suitable for a count-based model; rounding percentages does not
recover read counts.

Check the provenance and assay units, not just whether values happen to
be integers. Retain raw counts when they are available and choose a
method whose input requirements match the measurements. Missing values
are also distinct from observed zeros.

For QC and filtering, continue with [Quality control and
preprocessing](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/filtering_subsetting.md).
For formulas, target variables and engine compatibility, see
[Statistical models and
contrasts](https://microbialgenomics-irsicaixaorg.github.io/dar/articles/statistical_models.md).

## Session info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.1 (2026-06-24)`\
`#> Platform: x86_64-pc-linux-gnu`\
`#> Running under: Ubuntu 24.04.4 LTS`\
`#> `\
`#> Matrix products: default`\
`#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 `\
`#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0`\
`#> `\
`#> locale:`\
`#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              `\
`#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    `\
`#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   `\
`#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 `\
`#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            `\
`#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       `\
`#> `\
`#> time zone: UTC`\
`#> tzcode source: system (glibc)`\
`#> `\
`#> attached base packages:`\
`#> [1] stats     graphics  grDevices utils     datasets  methods   base     `\
`#> `\
`#> other attached packages:`\
`#> [1] dar_1.9.12`\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>   [1] DBI_1.3.0                       gridExtra_2.3.1                `\
`#>   [3] permute_0.9-10                  rlang_1.3.0                    `\
`#>   [5] magrittr_2.0.5                  ade4_1.7-24                    `\
`#>   [7] scater_1.41.2                   otel_0.2.0                     `\
`#>   [9] matrixStats_1.5.0               compiler_4.6.1                 `\
`#>  [11] DelayedMatrixStats_1.35.0       mgcv_1.9-4                     `\
`#>  [13] systemfonts_1.3.2               vctrs_0.7.3                    `\
`#>  [15] reshape2_1.4.5                  stringr_1.6.0                  `\
`#>  [17] pkgconfig_2.0.3                 crayon_1.5.3                   `\
`#>  [19] fastmap_1.2.0                   backports_1.5.1                `\
`#>  [21] XVector_0.53.0                  scuttle_1.23.2                 `\
`#>  [23] utf8_1.2.6                      rmarkdown_2.32                 `\
`#>  [25] ggbeeswarm_0.7.3                DirichletMultinomial_1.55.1    `\
`#>  [27] ragg_1.5.2                      purrr_1.2.2                    `\
`#>  [29] bluster_1.23.1                  MultiAssayExperiment_1.39.1    `\
`#>  [31] xfun_0.61                       cachem_1.1.0                   `\
`#>  [33] beachmat_2.29.3                 jsonlite_2.0.0                 `\
`#>  [35] biomformat_1.41.0               DelayedArray_0.39.8            `\
`#>  [37] BiocParallel_1.47.0             irlba_2.4.1                    `\
`#>  [39] parallel_4.6.1                  cluster_2.1.8.3                `\
`#>  [41] R6_2.6.1                        bslib_0.12.0                   `\
`#>  [43] stringi_1.8.9                   RColorBrewer_1.1-3             `\
`#>  [45] GenomicRanges_1.65.4            jquerylib_0.1.4                `\
`#>  [47] Rcpp_1.1.2                      Seqinfo_1.3.2                  `\
`#>  [49] SummarizedExperiment_1.43.0     iterators_1.0.14               `\
`#>  [51] knitr_1.52                      DECIPHER_3.9.3                 `\
`#>  [53] IRanges_2.47.5                  mia_1.21.11                    `\
`#>  [55] Matrix_1.7-6                    splines_4.6.1                  `\
`#>  [57] igraph_2.3.4                    tidyselect_1.2.1               `\
`#>  [59] viridis_0.6.5                   abind_1.4-8                    `\
`#>  [61] yaml_2.3.12                     vegan_2.7-6                    `\
`#>  [63] TreeSummarizedExperiment_2.21.0 codetools_0.2-20               `\
`#>  [65] lattice_0.23-1                  tibble_3.3.1                   `\
`#>  [67] plyr_1.8.9                      Biobase_2.73.2                 `\
`#>  [69] treeio_1.37.1                   S7_0.2.2                       `\
`#>  [71] evaluate_1.0.5                  desc_1.4.3                     `\
`#>  [73] survival_3.8-12                 Biostrings_2.81.9              `\
`#>  [75] pillar_1.11.1                   phyloseq_1.57.0                `\
`#>  [77] MatrixGenerics_1.25.0           checkmate_2.3.4                `\
`#>  [79] foreach_1.5.2                   stats4_4.6.1                   `\
`#>  [81] generics_0.1.4                  S4Vectors_0.51.10              `\
`#>  [83] ggplot2_4.0.3                   sparseMatrixStats_1.25.0       `\
`#>  [85] scales_1.4.0                    tidytree_0.4.8                 `\
`#>  [87] glue_1.8.1                      lazyeval_0.2.3                 `\
`#>  [89] tools_4.6.1                     BiocNeighbors_2.7.3            `\
`#>  [91] data.table_1.18.6.1             ScaledMatrix_1.21.0            `\
`#>  [93] fs_2.1.0                        grid_4.6.1                     `\
`#>  [95] tidyr_1.3.2                     ape_5.8-1                      `\
`#>  [97] SingleCellExperiment_1.35.2     nlme_3.1-171                   `\
`#>  [99] beeswarm_0.4.0                  BiocSingular_1.29.2            `\
`#> [101] vipor_0.4.7                     cli_3.6.6                      `\
`#> [103] rsvd_1.0.5                      rappdirs_0.3.4                 `\
`#> [105] textshaping_1.0.5               viridisLite_0.4.3              `\
`#> [107] S4Arrays_1.13.2                 dplyr_1.2.1                    `\
`#> [109] gtable_0.3.6                    yulab.utils_0.2.5              `\
`#> [111] sass_0.4.10                     digest_0.6.39                  `\
`#> [113] BiocGenerics_0.59.12            ggrepel_0.9.8                  `\
`#> [115] SparseArray_1.13.4              decontam_1.33.0                `\
`#> [117] htmlwidgets_1.6.4               farver_2.1.2                   `\
`#> [119] htmltools_0.5.9                 pkgdown_2.2.1                  `\
`#> [121] multtest_2.69.0                 lifecycle_1.0.5                `\
`#> [123] MASS_7.3-66`
