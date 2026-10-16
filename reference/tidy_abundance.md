# Extract abundances in the shared taxonomic analysis view

`tidy_abundance()` returns one row per taxon and sample. When a
centralized model is defined, taxa are aggregated through its selected
taxonomic rank, using the same identities as differential-abundance
results. Without a centralized model, the stored features are returned
without aggregation.

## Usage

``` r
tidy_abundance(rec, taxa_ids = NULL, transform = "identity", metadata = NULL)
```

## Arguments

- rec:

  A
  [Recipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
  or
  [PrepRecipe](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
  object.

- taxa_ids:

  Optional character vector of taxon identifiers to retain. IDs refer to
  the shared analysis view, not absorbed features. `NULL` retains every
  analysis taxon; `character(0)` returns no rows. Output keeps analysis
  taxon order rather than selector order.

- transform:

  Either `"identity"` for stored abundances or `"compositional"` for
  within-sample fractions. Fractions use all shared analysis taxa as
  their denominator before `taxa_ids` selection. Zero-depth samples have
  missing fractions and produce a warning.

- metadata:

  Character vector of sample metadata columns to include. `NULL`
  includes every column; `character(0)` includes none. Supply the
  original, unprefixed column names.

## Value

A tibble with `taxa_id`, `taxa`, `sample_id`, `abundance`, complete
namespaced taxonomy and selected namespaced sample metadata. `taxa` uses
the configured taxonomic label, or feature IDs when no single rank is
configured. Sample metadata types are preserved. Rows follow analysis
taxon order, then sample order.

## Details

This is a descriptive view of the data currently stored in the recipe. A
prepared recipe includes common preprocessing and the shared model
cohort; an unprepared recipe does not execute its configured steps.
Neither view reproduces method-specific filtering, rarefaction or
transformations, or adjusted model predictions.

Taxonomic columns are prefixed with `taxonomy_` and sample metadata
columns with `metadata_`, so annotations cannot overwrite identifiers or
abundances. Full lineage through the selected rank is retained; lower
ranks removed by aggregation remain missing. Taxa sharing a label but
not a lineage retain distinct `taxa_id` values.

## Examples

``` r
data(metaHIV_phy)
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus")
abundances <- tidy_abundance(rec, metadata = "RiskGroup2")
head(abundances)
#> # A tibble: 6 × 12
#>   taxa_id taxa              sample_id abundance taxonomy_Kingdom taxonomy_Phylum
#>   <chr>   <chr>             <chr>         <dbl> <chr>            <chr>          
#> 1 Otu_1   Methanobrevibact… Sample_1…      2538 Archaea          Euryarchaeota  
#> 2 Otu_1   Methanobrevibact… Sample_1…         0 Archaea          Euryarchaeota  
#> 3 Otu_1   Methanobrevibact… Sample_1…      9026 Archaea          Euryarchaeota  
#> 4 Otu_1   Methanobrevibact… Sample_1…     12049 Archaea          Euryarchaeota  
#> 5 Otu_1   Methanobrevibact… Sample_1…     22126 Archaea          Euryarchaeota  
#> 6 Otu_1   Methanobrevibact… Sample_1…     15374 Archaea          Euryarchaeota  
#> # ℹ 6 more variables: taxonomy_Class <chr>, taxonomy_Order <chr>,
#> #   taxonomy_Family <chr>, taxonomy_Genus <chr>, taxonomy_Species <chr>,
#> #   metadata_RiskGroup2 <chr>

fractions <- tidy_abundance(
  rec, transform = "compositional", metadata = character(0)
)
```
