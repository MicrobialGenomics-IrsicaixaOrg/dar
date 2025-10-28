# Loads Phyloseq data

This function returns a validated Phyloseq object by loading it directly
from a file with the .rds extension. Alternatively, this function can
also take three text files as input that will be used to construct and
validate the Phyloseq object: - Counts matrix with the otu_id in the
first column. - Taxonomic annotation matrix with the otu_id in the first
column. - Metadata annotation with sample_id in the first column.

## Usage

``` r
read_data(data_path)

validate_otu(otu)

validate_sample_data(sample_data)

validate_tax_table(tax_table)

validate_phyloseq(phy, slots = c("sample_data", "tax_table"))

read_phyloseq(file_path)

read_file(file_path, ext = c(".txt|.csv|.tsv"))
```

## Arguments

- data_path:

  List of length 1 or 3, with the paths of the input files.

## Value

a phylseq-class object

## Examples

``` r
# From a phyloseq object saved with .rds extension.
system.file("extdata", "metaHIV_phy.rds", package = "dar") |>
  read_data()
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 451 taxa and 156 samples ]
#> sample_data() Sample Data:       [ 156 samples by 3 sample variables ]
#> tax_table()   Taxonomy Table:    [ 451 taxa by 7 taxonomic ranks ]

# From the three components of a phyloseq object saved as a plain text.
data_path <- c(
  system.file("extdata", "metaHIV_counts.txt", package = "dar"),
  system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
  system.file("extdata", "metaHIV_taxas.txt", package = "dar")
)

read_data(data_path)
#> Rows: 451 Columns: 157
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr   (1): otu_id
#> dbl (156): Sample_186, Sample_185, Sample_184, Sample_182, Sample_181, Sampl...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 156 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (4): sample_id, RiskGroup2, Cluster, HIV_Status
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 451 Columns: 8
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (8): otu_id, Kingdom, Phylum, Class, Order, Family, Genus, Species
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 451 taxa and 156 samples ]
#> sample_data() Sample Data:       [ 156 samples by 3 sample variables ]
#> tax_table()   Taxonomy Table:    [ 451 taxa by 7 taxonomic ranks ]
```
