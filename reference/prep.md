# Performs all the steps defined in a Recipe

For a Recipe with at least one preprocessing or DA operation run the
steps in a convenient order.

## Usage

``` r
prep(
  rec,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16"),
  force = FALSE
)

# S4 method for class 'Recipe'
prep(
  rec,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16"),
  force = FALSE
)
```

## Arguments

- rec:

  A `Recipe` object. and furrr packages.

- parallel:

  if FALSE, no palatalization. if TRUE, parallel execution using future
  and furrr packages.

- workers:

  Number of workers for palatalization.

- force:

  Force the reexecution of all steps. This remove previous results.

## Value

A `PrepRecipe` object.

## Examples

``` r
data(metaHIV_phy)

## Define Recipe
rec <-
  recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Class") |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.03 * length(x))") |>
  step_maaslin()

## Prep Recipe
da_results <- prep(rec)

## If you try

## Consensus strategy
n_methods <- 2
da_results <- bake(da_results, count_cutoff = n_methods)
da_results
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Class 
#> 
#> Results:
#> 
#>      ✔ maaslin__Plăcintă diff_taxa = 53 
#> 
#>      ℹ 53 taxa are present in all tested methods 
#> 
#> Bakes:
#> 
#>      ◉ 1 -> count_cutoff: 2, weights: NULL, exclude: NULL, id: bake__Welsh_cake 

## If you try to run prep on an object of class PrepRecipe it returns an
## error.
err <- testthat::expect_error(prep(da_results))
err
#> <error/rlang_error>
#> Error in `prep()`:
#> ! The input Recipe has already been prep!
#> ℹ To force the rerun of all steps plese run prep(rec, force = T)
#> ---
#> Backtrace:
#>      ▆
#>   1. └─pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
#>   2.   └─pkgdown::build_site(...)
#>   3.     └─pkgdown:::build_site_local(...)
#>   4.       └─pkgdown::build_reference(...)
#>   5.         ├─pkgdown:::unwrap_purrr_error(...)
#>   6.         │ └─base::withCallingHandlers(...)
#>   7.         └─purrr::map(...)
#>   8.           └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#>   9.             ├─purrr:::with_indexed_errors(...)
#>  10.             │ └─base::withCallingHandlers(...)
#>  11.             ├─purrr:::call_with_cleanup(...)
#>  12.             └─pkgdown (local) .f(.x[[i]], ...)
#>  13.               ├─base::withCallingHandlers(...)
#>  14.               └─pkgdown:::data_reference_topic(...)
#>  15.                 └─pkgdown:::run_examples(...)
#>  16.                   └─pkgdown:::highlight_examples(code, topic, env = env)
#>  17.                     └─downlit::evaluate_and_highlight(...)
#>  18.                       └─evaluate::evaluate(code, child_env(env), new_device = TRUE, output_handler = output_handler)
#>  19.                         ├─base::withRestarts(...)
#>  20.                         │ └─base (local) withRestartList(expr, restarts)
#>  21.                         │   ├─base (local) withOneRestart(withRestartList(expr, restarts[-nr]), restarts[[nr]])
#>  22.                         │   │ └─base (local) doWithOneRestart(return(expr), restart)
#>  23.                         │   └─base (local) withRestartList(expr, restarts[-nr])
#>  24.                         │     └─base (local) withOneRestart(expr, restarts[[1L]])
#>  25.                         │       └─base (local) doWithOneRestart(return(expr), restart)
#>  26.                         ├─evaluate:::with_handlers(...)
#>  27.                         │ ├─base::eval(call)
#>  28.                         │ │ └─base::eval(call)
#>  29.                         │ └─base::withCallingHandlers(...)
#>  30.                         ├─base::withVisible(eval(expr, envir))
#>  31.                         └─base::eval(expr, envir)
#>  32.                           └─base::eval(expr, envir)
#>  33.                             ├─testthat::expect_error(prep(da_results))
#>  34.                             │ └─testthat:::expect_condition_matching_(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             ├─dar::prep(da_results)
#>  40.                             └─dar::prep(da_results)

## You can force the overwrite with:
prep(rec, force = TRUE)
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Class 
#> 
#> Results:
#> 
#>      ✔ maaslin__Plăcintă diff_taxa = 53 
#> 
#>      ℹ 53 taxa are present in all tested methods 
#> 

## This function can operate in parallel thanks to future and furrr packages.
prep(rec, parallel = TRUE, workers = 2)
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 278 taxa and 156 samples 
#>      ℹ variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
#>      ℹ taxonomic level Class 
#> 
#> Results:
#> 
#>      ✔ maaslin__Plăcintă diff_taxa = 53 
#> 
#>      ℹ 53 taxa are present in all tested methods 
#> 
```
