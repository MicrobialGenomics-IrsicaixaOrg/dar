# Import steps from json file

Import steps from json file

## Usage

``` r
import_steps(
  rec,
  file,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16")
)
```

## Arguments

- rec:

  A Recipe object.

- file:

  Path to the input file.

- parallel:

  if FALSE, no palatalization. if TRUE, parallel execution using future
  and furrr packages.

- workers:

  Number of workers for palatalization.

## Value

recipe-class object

## Examples

``` r
data(metaHIV_phy)

## Initialize the Recipe with a phyloseq object
rec <- recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> 

## Import steps
json_file <- system.file("extdata", "test.json", package = "dar")
rec <- import_steps(rec, json_file)
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Statistical model:
#> 
#>      ℹ ~RiskGroup2 
#> 
#> Preprocessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Viennoiserie 
#>      ◉ step_filter_taxa() id = filter_taxa__Karakudamono 
#> 
#> DA steps:
#> 
#>      ◉ step_maaslin() id = maaslin__Welsh_cake 
#>      ◉ step_deseq() id = deseq__Coussin_de_Lyon 

## If the json file contains 'bake', the Recipe is automatically prepared.
json_file <- system.file("extdata", "test_bake.json", package = "dar")
rec <-
  recipe(metaHIV_phy) |>
  add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
  import_steps(json_file)
#> ! bakes found in imported recipe
#> ℹ running `prep()`

rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 101 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: factor, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Results:
#> 
#>      ✔ maaslin__Welsh_cake diff_taxa = 36 
#>      ✔ deseq__Coussin_de_Lyon diff_taxa = 49 
#> 
#>      ℹ 29 taxon-contrast effects are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Carac 
cool(rec)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 29 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_102 Prevo… RiskGroup2… RiskGroup… main          Risk… up                2
#>  2 Otu_115 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  3 Otu_115 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  4 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  5 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  6 Otu_129 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  7 Otu_255 Rumin… RiskGroup2… RiskGroup… main          Risk… up                2
#>  8 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                2
#>  9 Otu_261 Dorea… RiskGroup2… RiskGroup… main          Risk… up                2
#> 10 Otu_262 Dorea… RiskGroup2… RiskGroup… main          Risk… up                2
#> # ℹ 19 more rows
#> # ℹ 1 more variable: methods <chr>
```
