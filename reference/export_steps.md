# Export step parameters as json.

Export step parameters as json.

## Usage

``` r
export_steps(rec, file_name)
```

## Arguments

- rec:

  A Recipe object.

- file_name:

  The path and file name of the optout file.

## Value

invisible

## Examples

``` r
data(metaHIV_phy)

## Create a Recipe with steps
rec <- 
  recipe(metaHIV_phy, "RiskGroup2", "Species") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") |>
  step_filter_by_prevalence(0.4) |>
  step_maaslin()
 
## Prep Recipe   
rec <- prep(rec, parallel = TRUE)

## Export to json file
export_steps(rec, tempfile(fileext = ".json"))
```
