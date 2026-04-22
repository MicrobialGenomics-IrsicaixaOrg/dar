# Converting Common Data Formats to Phyloseq and TreeSummarizedExperiment

The “dar” package is a versatile and user-friendly tool designed to
accept inputs in a variety of formats. It primarily utilizes the
`phyloseq` format but also supports the `TreeSummarizedExperiment`
format. This flexibility allows users to conduct differential abundance
analysis smoothly, irrespective of their initial data format. To
facilitate this, a detailed guide is available to aid users in
converting other prevalent data formats, such as `biome`, `mothur`,
`metaphlan`, and more, into the necessary `phyloseq` or
`TreeSummarizedExperiment` formats.

``` r

suppressPackageStartupMessages(library(mia))
suppressPackageStartupMessages(library(phyloseq))
```

## Importing Data from `biome` Format

The `biome` format is a commonly used format in bioinformatics to
represent microbiome sequencing data. Here’s how you can import data in
`biome` format to both `phyloseq` and `TreeSummarizedExperiment.`

### To Phyloseq

To convert data from the `biome` format to the `phyloseq` format, you
can use the
[`phyloseq::import_biom()`](https://rdrr.io/pkg/phyloseq/man/import_biom.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Example of a rich dense biom file
rich_dense_biom  <-
  system.file("extdata", "rich_dense_otu_table.biom",  package = "phyloseq")

# Import biom as a phyloseq-class object
phy <- phyloseq::import_biom(
  rich_dense_biom,  
  parseFunction = parse_taxonomy_greengenes
)

phy
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 5 taxa and 6 samples ]
#> sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
#> tax_table()   Taxonomy Table:    [ 5 taxa by 7 taxonomic ranks ]

# Print sample_data
phyloseq::sample_data(phy)
#>         BarcodeSequence  LinkerPrimerSequence BODY_SITE Description
#> Sample1    CGCTTATCGAGA CATGCTGCCTCCCGTAGGAGT       gut   human gut
#> Sample2    CATACCAGTAGC CATGCTGCCTCCCGTAGGAGT       gut   human gut
#> Sample3    CTCTCTACCTGT CATGCTGCCTCCCGTAGGAGT       gut   human gut
#> Sample4    CTCTCGGCCTGT CATGCTGCCTCCCGTAGGAGT      skin  human skin
#> Sample5    CTCTCTACCAAT CATGCTGCCTCCCGTAGGAGT      skin  human skin
#> Sample6    CTAACTACCAAT CATGCTGCCTCCCGTAGGAGT      skin  human skin

# Print tax_table
phyloseq::tax_table(phy)
#> Taxonomy Table:     [5 taxa by 7 taxonomic ranks]:
#>          Kingdom    Phylum           Class                 Order              
#> GG_OTU_1 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacteriales"
#> GG_OTU_2 "Bacteria" "Cyanobacteria"  "Nostocophycideae"    "Nostocales"       
#> GG_OTU_3 "Archaea"  "Euryarchaeota"  "Methanomicrobia"     "Methanosarcinales"
#> GG_OTU_4 "Bacteria" "Firmicutes"     "Clostridia"          "Halanaerobiales"  
#> GG_OTU_5 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacteriales"
#>          Family               Genus            Species                       
#> GG_OTU_1 "Enterobacteriaceae" "Escherichia"    NA                            
#> GG_OTU_2 "Nostocaceae"        "Dolichospermum" NA                            
#> GG_OTU_3 "Methanosarcinaceae" "Methanosarcina" NA                            
#> GG_OTU_4 "Halanaerobiaceae"   "Halanaerobium"  "Halanaerobiumsaccharolyticum"
#> GG_OTU_5 "Enterobacteriaceae" "Escherichia"    NA

# Recipe init
rec <- dar::recipe(phy, var_info = "BODY_SITE", tax_info = "Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 5 taxa and 6 samples 
#>      ℹ variable of interes BODY_SITE (class: character, levels: gut, skin) 
#>      ℹ taxonomic level Genus
```

### To TreeSummarizedExperiment

To convert data from the `biome` format to the
`TreeSummarizedExperiment` format, you can use the
[`mia::importBIOM()`](https://microbiome.github.io/mia/reference/importBIOM.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Example of a rich dense biom file
rich_dense_biom  <-
  system.file("extdata", "rich_dense_otu_table.biom",  package = "phyloseq")

# Import biom as a phyloseq-class object
tse <- mia::importBIOM(rich_dense_biom)

tse
#> class: TreeSummarizedExperiment 
#> dim: 5 6 
#> metadata(0):
#> assays(1): counts
#> rownames(5): GG_OTU_1 GG_OTU_2 GG_OTU_3 GG_OTU_4 GG_OTU_5
#> rowData names(7): taxonomy1 taxonomy2 ... taxonomy6 taxonomy7
#> colnames(6): Sample1 Sample2 ... Sample5 Sample6
#> colData names(4): BarcodeSequence LinkerPrimerSequence BODY_SITE
#>   Description
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# Print sample_data
colData(tse)
#> DataFrame with 6 rows and 4 columns
#>         BarcodeSequence  LinkerPrimerSequence   BODY_SITE Description
#>             <character>           <character> <character> <character>
#> Sample1    CGCTTATCGAGA CATGCTGCCTCCCGTAGGAGT         gut   human gut
#> Sample2    CATACCAGTAGC CATGCTGCCTCCCGTAGGAGT         gut   human gut
#> Sample3    CTCTCTACCTGT CATGCTGCCTCCCGTAGGAGT         gut   human gut
#> Sample4    CTCTCGGCCTGT CATGCTGCCTCCCGTAGGAGT        skin  human skin
#> Sample5    CTCTCTACCAAT CATGCTGCCTCCCGTAGGAGT        skin  human skin
#> Sample6    CTAACTACCAAT CATGCTGCCTCCCGTAGGAGT        skin  human skin

# Print tax_table
rowData(tse)
#> DataFrame with 5 rows and 7 columns
#>            taxonomy1         taxonomy2              taxonomy3
#>          <character>       <character>            <character>
#> GG_OTU_1 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#> GG_OTU_2 k__Bacteria  p__Cyanobacteria    c__Nostocophycideae
#> GG_OTU_3  k__Archaea  p__Euryarchaeota     c__Methanomicrobia
#> GG_OTU_4 k__Bacteria     p__Firmicutes          c__Clostridia
#> GG_OTU_5 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#>                     taxonomy4             taxonomy5         taxonomy6
#>                   <character>           <character>       <character>
#> GG_OTU_1 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
#> GG_OTU_2        o__Nostocales        f__Nostocaceae g__Dolichospermum
#> GG_OTU_3 o__Methanosarcinales f__Methanosarcinaceae g__Methanosarcina
#> GG_OTU_4   o__Halanaerobiales   f__Halanaerobiaceae  g__Halanaerobium
#> GG_OTU_5 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
#>                       taxonomy7
#>                     <character>
#> GG_OTU_1                    s__
#> GG_OTU_2                    s__
#> GG_OTU_3                    s__
#> GG_OTU_4 s__Halanaerobiumsacc..
#> GG_OTU_5                    s__

# Change the column names of the tax_table
colnames(rowData(tse)) <- 
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rowData(tse)
#> DataFrame with 5 rows and 7 columns
#>              Kingdom            Phylum                  Class
#>          <character>       <character>            <character>
#> GG_OTU_1 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#> GG_OTU_2 k__Bacteria  p__Cyanobacteria    c__Nostocophycideae
#> GG_OTU_3  k__Archaea  p__Euryarchaeota     c__Methanomicrobia
#> GG_OTU_4 k__Bacteria     p__Firmicutes          c__Clostridia
#> GG_OTU_5 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#>                         Order                Family             Genus
#>                   <character>           <character>       <character>
#> GG_OTU_1 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
#> GG_OTU_2        o__Nostocales        f__Nostocaceae g__Dolichospermum
#> GG_OTU_3 o__Methanosarcinales f__Methanosarcinaceae g__Methanosarcina
#> GG_OTU_4   o__Halanaerobiales   f__Halanaerobiaceae  g__Halanaerobium
#> GG_OTU_5 o__Enterobacteriales f__Enterobacteriaceae    g__Escherichia
#>                         Species
#>                     <character>
#> GG_OTU_1                    s__
#> GG_OTU_2                    s__
#> GG_OTU_3                    s__
#> GG_OTU_4 s__Halanaerobiumsacc..
#> GG_OTU_5                    s__

# Recipe init
rec <- dar::recipe(tse, var_info = "BODY_SITE", tax_info = "Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 5 taxa and 6 samples 
#>      ℹ variable of interes BODY_SITE (class: character, levels: gut, skin) 
#>      ℹ taxonomic level Genus
```

## Importing Data from `qiime` Format

The `qiime` format is another commonly used format in bioinformatics for
microbiome sequencing data. Here’s how you can import data in `qiime`
format to both `Phyloseq` and `TreeSummarizedExperiment.`

### To Phyloseq

To convert data from the `qiime` format to the `Phyloseq` format, you
can use the
[`phyloseq::import_qiime()`](https://rdrr.io/pkg/phyloseq/man/import_qiime.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Import QIIME data
phy_qiime <- phyloseq::import_qiime(
  otufilename = system.file("extdata", "GP_otu_table_rand_short.txt.gz", package = "phyloseq"),
  mapfilename = system.file("extdata", "master_map.txt", package = "phyloseq"),
  treefilename = system.file("extdata", "GP_tree_rand_short.newick.gz", package = "phyloseq")
)
#> Processing map file...
#> Processing otu/tax file...
#> Reading file into memory prior to parsing...
#> Detecting first header line...
#> Header is on line 2  
#> Converting input file to a table...
#> Defining OTU table... 
#> Parsing taxonomy table...
#> Processing phylogenetic tree...
#>  /__w/_temp/Library/phyloseq/extdata/GP_tree_rand_short.newick.gz ...

phy_qiime
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 500 taxa and 26 samples ]
#> sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 500 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 500 tips and 499 internal nodes ]

# Recipe init
rec <- dar::recipe(phy_qiime, var_info = "SampleType", tax_info = "Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 500 taxa and 26 samples 
#>      ℹ variable of interes SampleType (class: character, levels: Feces, Freshwater, Freshwater (creek), Mock, Ocean, Sediment (estuary), Skin, Soil, Tongue) 
#>      ℹ taxonomic level Genus
```

### To TreeSummarizedExperiment

To convert data from the `qiime` format to the
`TreeSummarizedExperiment` format, you can use the
[`mia::importQIIME2()`](https://microbiome.github.io/mia/reference/importQIIME2.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Import QIIME data to tse
tse_qiime <- mia::importQIIME2(
  featureTableFile = system.file("extdata", "table.qza", package = "mia"),
  taxonomyTableFile = system.file("extdata", "taxonomy.qza", package = "mia"),
  sampleMetaFile = system.file("extdata", "sample-metadata.tsv", package = "mia"),
  refSeqFile = system.file("extdata", "refseq.qza", package = "mia"),
  phyTreeFile = system.file("extdata", "tree.qza", package = "mia")
)

tse_qiime

# Recipe init
rec <- dar::recipe(tse_qiime, var_info = "body.site", tax_info = "genus")

rec
```

## Importing Data from `mothur` Format

The `mothur` format is another commonly used format in bioinformatics
for microbiome sequencing data. Here’s how you can import data in
`mothur` format to both `Phyloseq` and `TreeSummarizedExperiment.`

### To Phyloseq

To convert data from the `mothur` format to the `Phyloseq` format, you
can use the
[`phyloseq::import_mothur()`](https://rdrr.io/pkg/phyloseq/man/import_mothur.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Import Mothur data
phy_mothur <- phyloseq::import_mothur(
  mothur_list_file = system.file("extdata", "esophagus.fn.list.gz", package = "phyloseq"),
  mothur_group_file = system.file("extdata", "esophagus.good.groups.gz", package = "phyloseq"),
  mothur_tree_file = system.file("extdata", "esophagus.tree.gz", package = "phyloseq")
)

phy_mothur
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 591 taxa and 3 samples ]
#> phy_tree()    Phylogenetic Tree: [ 591 tips and 590 internal nodes ]

# Recipe init
rec <- dar::recipe(phy_mothur)

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 591 taxa and 3 samples 
#>      ✖ undefined variable of interest. Use add_var() to add it to Recipe! 
#>      ✖ undefined taxonomic level. Use add_tax() to add it to Recipe!
```

### To TreeSummarizedExperiment

To convert data from the `mothur` format to the
`TreeSummarizedExperiment` format, you can use the
[`mia::importMothur()`](https://microbiome.github.io/mia/reference/importMothur.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Import Mothur data to TreeSummarizedExperiment
tse_mothur <- mia::importMothur(
  sharedFile = system.file("extdata", "mothur_example.shared", package = "mia"),
  taxonomyFile = system.file("extdata", "mothur_example.cons.taxonomy", package = "mia"),
  designFile = system.file("extdata", "mothur_example.design", package = "mia")
) |> methods::as("TreeSummarizedExperiment")

tse_mothur
#> class: TreeSummarizedExperiment 
#> dim: 100 100 
#> metadata(0):
#> assays(1): counts
#> rownames(100): Otu001 Otu002 ... Otu099 Otu100
#> rowData names(8): OTU Size ... Family Genus
#> colnames(100): Sample1 Sample2 ... Sample99 Sample100
#> colData names(7): group sex ... numOtus Group
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# Recipe init
rec <- dar::recipe(tse_mothur, var_info = "drug", tax_info = "Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 100 taxa and 100 samples 
#>      ℹ variable of interes drug (class: character, levels: A, B) 
#>      ℹ taxonomic level Genus
```

## Importing Data from `metaphlan` Format

The `metaphlan` format is another commonly used format in bioinformatics
for microbiome sequencing data. Here’s how you can import data in
`metaphlan` format to `TreeSummarizedExperiment.`

### To TreeSummarizedExperiment

To convert data from the `metaphlan` format to the
`TreeSummarizedExperiment` format, you can use the
[`mia::importMetaPhlAn()`](https://microbiome.github.io/mia/reference/importMetaPhlAn.html)
function. Here’s a step-by-step example of how to perform this
conversion:

``` r

# Importing data from Metaphlan
tse_metaphlan <- mia::importMetaPhlAn(
  file = system.file("extdata", "merged_abundance_table.txt", package = "mia")
)

# Recipe init
tse_metaphlan <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = SummarizedExperiment::assay(tse_metaphlan)),
  rowData = SummarizedExperiment::rowData(tse_metaphlan),
  colData = SummarizedExperiment::colData(tse_metaphlan) |> 
    as.data.frame() |> 
    dplyr::mutate(condition = rep(c("A", "B"), times = 3))
)

rec <- dar::recipe(tse_metaphlan, var_info = "condition", tax_info = "Genus")

rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 16 taxa and 6 samples 
#>      ℹ variable of interes condition (class: character, levels: A, B) 
#>      ℹ taxonomic level Genus
```

## Conclusion

In this guide, we have explored various methods for importing microbiome
sequencing data from different formats into `Phyloseq` and
`TreeSummarizedExperiment.` We’ve covered the `biome`, `qiime`,
`mothur`, `metaphlan`, and `humann` formats, providing step-by-step
examples for each.

The flexibility of these tools allows for a smooth transition between
different data formats, making it easier to conduct your analysis
irrespective of the initial data format. By following the steps outlined
in this guide, you should be able to successfully convert your data and
carry out further differential abundance analysis.

Remember, the specific details of your data may require you to adjust
the parameters in the import functions. Always inspect your data after
conversion to ensure it has been imported correctly.

## Session info

``` r

devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.2 (2025-10-31)
#>  os       Ubuntu 24.04.3 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       UTC
#>  date     2026-04-22
#>  pandoc   3.8.2.1 @ /usr/bin/ (via rmarkdown)
#>  quarto   1.7.32 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package                  * version  date (UTC) lib source
#>  abind                      1.4-8    2024-09-12 [1] RSPM (R 4.5.0)
#>  ade4                       1.7-24   2026-03-21 [1] RSPM (R 4.5.0)
#>  ape                        5.8-1    2024-12-16 [1] RSPM (R 4.5.0)
#>  assertthat                 0.2.1    2019-03-21 [1] RSPM (R 4.5.0)
#>  backports                  1.5.1    2026-04-03 [1] RSPM (R 4.5.0)
#>  beachmat                   2.26.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  beeswarm                   0.4.0    2021-06-01 [1] RSPM (R 4.5.0)
#>  Biobase                  * 2.70.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocBaseUtils              1.12.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocGenerics             * 0.56.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocNeighbors              2.4.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocParallel               1.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  BiocSingular               1.26.1   2025-11-17 [1] Bioconductor 3.22 (R 4.5.2)
#>  biomformat                 1.38.3   2026-03-16 [1] Bioconductor 3.22 (R 4.5.2)
#>  Biostrings               * 2.78.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  bluster                    1.20.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  bslib                      0.10.0   2026-01-26 [2] RSPM (R 4.5.0)
#>  ca                         0.71.1   2020-01-24 [1] RSPM (R 4.5.0)
#>  cachem                     1.1.0    2024-05-16 [2] RSPM (R 4.5.0)
#>  cellranger                 1.1.0    2016-07-27 [1] RSPM (R 4.5.0)
#>  checkmate                  2.3.4    2026-02-03 [1] RSPM (R 4.5.0)
#>  cli                        3.6.6    2026-04-09 [2] RSPM (R 4.5.0)
#>  cluster                    2.1.8.2  2026-02-05 [3] RSPM (R 4.5.0)
#>  coda                       0.19-4.1 2024-01-31 [1] RSPM (R 4.5.0)
#>  codetools                  0.2-20   2024-03-31 [3] CRAN (R 4.5.2)
#>  crayon                     1.5.3    2024-06-20 [2] RSPM (R 4.5.0)
#>  dar                        1.7.1    2026-04-22 [1] Bioconductor
#>  data.table                 1.18.2.1 2026-01-27 [1] RSPM (R 4.5.0)
#>  DBI                        1.3.0    2026-02-25 [1] RSPM (R 4.5.0)
#>  DECIPHER                   3.6.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  decontam                   1.30.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  DelayedArray               0.36.1   2026-03-31 [1] Bioconductor 3.22 (R 4.5.2)
#>  DelayedMatrixStats         1.32.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  dendextend                 1.19.1   2025-07-15 [1] RSPM (R 4.5.0)
#>  desc                       1.4.3    2023-12-10 [2] RSPM (R 4.5.0)
#>  devtools                   2.5.1    2026-04-16 [2] RSPM (R 4.5.0)
#>  digest                     0.6.39   2025-11-19 [2] RSPM (R 4.5.0)
#>  DirichletMultinomial       1.52.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  dplyr                      1.2.1    2026-04-03 [1] RSPM (R 4.5.0)
#>  ellipsis                   0.3.3    2026-04-04 [2] RSPM (R 4.5.0)
#>  emmeans                    2.0.3    2026-04-09 [1] RSPM (R 4.5.0)
#>  estimability               1.5.1    2024-05-12 [1] RSPM (R 4.5.0)
#>  evaluate                   1.0.5    2025-08-27 [2] RSPM (R 4.5.0)
#>  farver                     2.1.2    2024-05-13 [1] RSPM (R 4.5.0)
#>  fastmap                    1.2.0    2024-05-15 [2] RSPM (R 4.5.0)
#>  fillpattern                1.0.3    2026-02-13 [1] RSPM (R 4.5.0)
#>  foreach                    1.5.2    2022-02-02 [1] RSPM (R 4.5.0)
#>  fs                         2.1.0    2026-04-18 [2] RSPM (R 4.5.0)
#>  generics                 * 0.1.4    2025-05-09 [1] RSPM (R 4.5.0)
#>  GenomicRanges            * 1.62.1   2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
#>  ggbeeswarm                 0.7.3    2025-11-29 [1] RSPM (R 4.5.0)
#>  ggnewscale                 0.5.2    2025-06-20 [1] RSPM (R 4.5.0)
#>  ggplot2                    4.0.2    2026-02-03 [1] RSPM (R 4.5.0)
#>  ggrepel                    0.9.8    2026-03-17 [1] RSPM (R 4.5.0)
#>  ggtext                     0.1.2    2022-09-16 [1] RSPM (R 4.5.0)
#>  glue                       1.8.1    2026-04-17 [2] RSPM (R 4.5.0)
#>  gridExtra                  2.3      2017-09-09 [1] RSPM (R 4.5.0)
#>  gridtext                   0.1.6    2026-02-19 [1] RSPM (R 4.5.0)
#>  gtable                     0.3.6    2024-10-25 [1] RSPM (R 4.5.0)
#>  heatmaply                  1.6.0    2025-07-12 [1] RSPM (R 4.5.0)
#>  hms                        1.1.4    2025-10-17 [1] RSPM (R 4.5.0)
#>  htmltools                  0.5.9    2025-12-04 [2] RSPM (R 4.5.0)
#>  htmlwidgets                1.6.4    2023-12-06 [2] RSPM (R 4.5.0)
#>  httr                       1.4.8    2026-02-13 [1] RSPM (R 4.5.0)
#>  igraph                     2.3.0    2026-04-21 [1] RSPM (R 4.5.0)
#>  IRanges                  * 2.44.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  irlba                      2.3.7    2026-01-30 [1] RSPM (R 4.5.0)
#>  iterators                  1.0.14   2022-02-05 [1] RSPM (R 4.5.0)
#>  jquerylib                  0.1.4    2021-04-26 [2] RSPM (R 4.5.0)
#>  jsonlite                   2.0.0    2025-03-27 [2] RSPM (R 4.5.0)
#>  knitr                      1.51     2025-12-20 [2] RSPM (R 4.5.0)
#>  lattice                    0.22-9   2026-02-09 [3] RSPM (R 4.5.0)
#>  lazyeval                   0.2.3    2026-04-04 [1] RSPM (R 4.5.0)
#>  lifecycle                  1.0.5    2026-01-08 [2] RSPM (R 4.5.0)
#>  magrittr                   2.0.5    2026-04-04 [2] RSPM (R 4.5.0)
#>  MASS                       7.3-65   2025-02-28 [3] CRAN (R 4.5.2)
#>  Matrix                     1.7-5    2026-03-21 [3] RSPM (R 4.5.0)
#>  MatrixGenerics           * 1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  matrixStats              * 1.5.0    2025-01-07 [1] RSPM (R 4.5.0)
#>  memoise                    2.0.1    2021-11-26 [2] RSPM (R 4.5.0)
#>  mgcv                       1.9-4    2025-11-07 [3] RSPM (R 4.5.0)
#>  mia                      * 1.18.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  multcomp                   1.4-30   2026-03-09 [1] RSPM (R 4.5.0)
#>  MultiAssayExperiment     * 1.36.2   2026-03-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  multtest                   2.66.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  mvtnorm                    1.3-7    2026-04-15 [1] RSPM (R 4.5.0)
#>  nlme                       3.1-169  2026-03-27 [3] RSPM (R 4.5.0)
#>  otel                       0.2.0    2025-08-29 [2] RSPM (R 4.5.0)
#>  parallelly                 1.47.0   2026-04-17 [1] RSPM (R 4.5.0)
#>  patchwork                  1.3.2    2025-08-25 [1] RSPM (R 4.5.0)
#>  permute                    0.9-10   2026-02-06 [1] RSPM (R 4.5.0)
#>  phyloseq                 * 1.54.2   2026-03-02 [1] Bioconductor 3.22 (R 4.5.2)
#>  pillar                     1.11.1   2025-09-17 [2] RSPM (R 4.5.0)
#>  pkgbuild                   1.4.8    2025-05-26 [2] RSPM (R 4.5.0)
#>  pkgconfig                  2.0.3    2019-09-22 [2] RSPM (R 4.5.0)
#>  pkgdown                    2.2.0    2025-11-06 [2] RSPM (R 4.5.0)
#>  pkgload                    1.5.1    2026-04-01 [2] RSPM (R 4.5.0)
#>  plotly                     4.12.0   2026-01-24 [1] RSPM (R 4.5.0)
#>  plyr                       1.8.9    2023-10-02 [1] RSPM (R 4.5.0)
#>  purrr                      1.2.2    2026-04-10 [2] RSPM (R 4.5.0)
#>  R6                         2.6.1    2025-02-15 [2] RSPM (R 4.5.0)
#>  ragg                       1.5.2    2026-03-23 [2] RSPM (R 4.5.0)
#>  rappdirs                   0.3.4    2026-01-17 [2] RSPM (R 4.5.0)
#>  rbiom                      2.2.1    2025-06-27 [1] RSPM (R 4.5.0)
#>  RColorBrewer               1.1-3    2022-04-03 [1] RSPM (R 4.5.0)
#>  Rcpp                       1.1.1-1  2026-04-16 [2] RSPM (R 4.5.0)
#>  readr                      2.2.0    2026-02-19 [1] RSPM (R 4.5.0)
#>  readxl                     1.4.5    2025-03-07 [1] RSPM (R 4.5.0)
#>  registry                   0.5-1    2019-03-05 [1] RSPM (R 4.5.0)
#>  reshape2                   1.4.5    2025-11-12 [1] RSPM (R 4.5.0)
#>  rlang                      1.2.0    2026-04-06 [2] RSPM (R 4.5.0)
#>  rmarkdown                  2.31     2026-03-26 [2] RSPM (R 4.5.0)
#>  rsvd                       1.0.5    2021-04-16 [1] RSPM (R 4.5.0)
#>  S4Arrays                   1.10.1   2025-12-01 [1] Bioconductor 3.22 (R 4.5.2)
#>  S4Vectors                * 0.48.1   2026-04-05 [1] Bioconductor 3.22 (R 4.5.2)
#>  S7                         0.2.2    2026-04-22 [1] CRAN (R 4.5.2)
#>  sandwich                   3.1-1    2024-09-15 [1] RSPM (R 4.5.0)
#>  sass                       0.4.10   2025-04-11 [2] RSPM (R 4.5.0)
#>  ScaledMatrix               1.18.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  scales                     1.4.0    2025-04-24 [1] RSPM (R 4.5.0)
#>  scater                     1.38.1   2026-03-20 [1] Bioconductor 3.22 (R 4.5.2)
#>  scuttle                    1.20.0   2025-10-30 [1] Bioconductor 3.22 (R 4.5.2)
#>  Seqinfo                  * 1.0.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  seriation                  1.5.8    2025-08-20 [1] RSPM (R 4.5.0)
#>  sessioninfo                1.2.3    2025-02-05 [2] RSPM (R 4.5.0)
#>  SingleCellExperiment     * 1.32.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  slam                       0.1-55   2024-11-13 [1] RSPM (R 4.5.0)
#>  SparseArray                1.10.10  2026-03-30 [1] Bioconductor 3.22 (R 4.5.2)
#>  sparseMatrixStats          1.22.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  stringi                    1.8.7    2025-03-27 [2] RSPM (R 4.5.0)
#>  stringr                    1.6.0    2025-11-04 [2] RSPM (R 4.5.0)
#>  SummarizedExperiment     * 1.40.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  survival                   3.8-6    2026-01-16 [3] RSPM (R 4.5.0)
#>  systemfonts                1.3.2    2026-03-05 [2] RSPM (R 4.5.0)
#>  textshaping                1.0.5    2026-03-06 [2] RSPM (R 4.5.0)
#>  TH.data                    1.1-5    2025-11-17 [1] RSPM (R 4.5.0)
#>  tibble                     3.3.1    2026-01-11 [2] RSPM (R 4.5.0)
#>  tidyr                      1.3.2    2025-12-19 [1] RSPM (R 4.5.0)
#>  tidyselect                 1.2.1    2024-03-11 [1] RSPM (R 4.5.0)
#>  tidytree                   0.4.7    2026-01-08 [1] RSPM (R 4.5.0)
#>  treeio                     1.34.0   2025-10-30 [1] Bioconductor 3.22 (R 4.5.2)
#>  TreeSummarizedExperiment * 2.18.0   2025-10-30 [1] Bioconductor 3.22 (R 4.5.2)
#>  TSP                        1.2.7    2026-03-23 [1] RSPM (R 4.5.0)
#>  tzdb                       0.5.0    2025-03-15 [1] RSPM (R 4.5.0)
#>  UpSetR                     1.4.0    2019-05-22 [1] RSPM (R 4.5.0)
#>  usethis                    3.2.1    2025-09-06 [2] RSPM (R 4.5.0)
#>  vctrs                      0.7.3    2026-04-11 [2] RSPM (R 4.5.0)
#>  vegan                      2.7-3    2026-03-04 [1] RSPM (R 4.5.0)
#>  vipor                      0.4.7    2023-12-18 [1] RSPM (R 4.5.0)
#>  viridis                    0.6.5    2024-01-29 [1] RSPM (R 4.5.0)
#>  viridisLite                0.4.3    2026-02-04 [1] RSPM (R 4.5.0)
#>  webshot                    0.5.5    2023-06-26 [1] RSPM (R 4.5.0)
#>  withr                      3.0.2    2024-10-28 [2] RSPM (R 4.5.0)
#>  xfun                       0.57     2026-03-20 [2] RSPM (R 4.5.0)
#>  xml2                       1.5.2    2026-01-17 [2] RSPM (R 4.5.0)
#>  xtable                     1.8-8    2026-02-22 [2] RSPM (R 4.5.0)
#>  XVector                  * 0.50.0   2025-10-29 [1] Bioconductor 3.22 (R 4.5.2)
#>  yaml                       2.3.12   2025-12-10 [2] RSPM (R 4.5.0)
#>  yulab.utils                0.2.4    2026-02-02 [1] RSPM (R 4.5.0)
#>  zoo                        1.8-15   2025-12-15 [1] RSPM (R 4.5.0)
#> 
#>  [1] /__w/_temp/Library
#>  [2] /usr/local/lib/R/site-library
#>  [3] /usr/local/lib/R/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
