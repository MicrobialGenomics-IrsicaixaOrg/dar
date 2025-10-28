# Import steps from json file

Import steps from json file

## Usage

``` r
import_steps(rec, file, parallel = TRUE, workers = future::availableCores())
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
rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
rec
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
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
#> Preporcessing steps:
#> 
#>      ◉ step_subset_taxa() id = subset_taxa__Boyoz 
#>      ◉ step_filter_taxa() id = filter_taxa__Charlotte 
#> 
#> DA steps:
#> 
#>      ◉ step_maaslin() id = maaslin__ChaSan 
#>      ◉ step_metagenomeseq() id = metagenomeseq__Zlebia 
#>      ◉ step_deseq() id = deseq__Linzer_torte 

## If the json file contains 'bake', the Recipe is automatically prepared. 
json_file <- system.file("extdata", "test_bake.json", package = "dar")
rec <- 
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  import_steps(json_file)
#> ! bakes found in imported recipe
#> ℹ running prep()
  
rec
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 101 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Species 
#> 
#> Results:
#> 
#>      ✔ maaslin__ChaSan diff_taxa = 86 
#>      ✔ metagenomeseq__Zlebia diff_taxa = 69 
#>      ✔ deseq__Linzer_torte diff_taxa = 53 
#> 
#>      ℹ 46 taxa are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Trdelník 
cool(rec)
#> ℹ Bake for count_cutoff = 3
#> # A tibble: 46 × 2
#>    taxa_id taxa                        
#>    <chr>   <chr>                       
#>  1 Otu_96  Prevotella_copri            
#>  2 Otu_102 Prevotella_sp_AM42_24       
#>  3 Otu_340 Holdemanella_biformis       
#>  4 Otu_332 Catenibacterium_mitsuokai   
#>  5 Otu_108 Prevotella_sp_CAG_520       
#>  6 Otu_115 Alistipes_finegoldii        
#>  7 Otu_119 Alistipes_putredinis        
#>  8 Otu_45  Slackia_isoflavoniconvertens
#>  9 Otu_82  Barnesiella_intestinihominis
#> 10 Otu_113 Prevotella_stercorea        
#> # ℹ 36 more rows
```
