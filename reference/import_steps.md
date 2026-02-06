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
#>      ◉ step_subset_taxa() id = subset_taxa__Bear_claw 
#>      ◉ step_filter_taxa() id = filter_taxa__Spanakopita 
#> 
#> DA steps:
#> 
#>      ◉ step_maaslin() id = maaslin__Eccles_cake 
#>      ◉ step_deseq() id = deseq__Belekoy 

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
#>      ✔ maaslin__Eccles_cake diff_taxa = 41 
#>      ✔ deseq__Belekoy diff_taxa = 53 
#> 
#>      ℹ 27 taxa are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Makmur 
cool(rec)
#> ℹ Bake for count_cutoff = 2
#> # A tibble: 27 × 2
#>    taxa_id taxa                        
#>    <chr>   <chr>                       
#>  1 Otu_78  Bacteroides_uniformis       
#>  2 Otu_88  Odoribacter_splanchnicus    
#>  3 Otu_119 Alistipes_putredinis        
#>  4 Otu_129 Parabacteroides_merdae      
#>  5 Otu_125 Parabacteroides_distasonis  
#>  6 Otu_96  Prevotella_copri            
#>  7 Otu_82  Barnesiella_intestinihominis
#>  8 Otu_51  Bacteroides_dorei           
#>  9 Otu_52  Bacteroides_eggerthii       
#> 10 Otu_332 Catenibacterium_mitsuokai   
#> # ℹ 17 more rows
```
