# Extract results from defined bake

Extract results from defined bake

## Usage

``` r
cool(rec, bake = 1)

# S4 method for class 'Recipe'
cool(rec, bake = 1)

# S4 method for class 'PrepRecipe'
cool(rec, bake = 1)
```

## Arguments

- rec:

  A `Recipe` object.

- bake:

  Name or index of the bake to extract.

## Value

tbl_df

## Examples

``` r
data(test_prep_rec)

## First we need to add bakes (extraction strategies) to the PrepRecipe.
rec <- bake(test_prep_rec)

## Finally we can extract the results with the cool method
cool(rec)
#> ℹ Bake for count_cutoff = 3
#> # A tibble: 23 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_96  Prevotella     
#>  2 Otu_340 Holdemanella   
#>  3 Otu_332 Catenibacterium
#>  4 Otu_119 Alistipes      
#>  5 Otu_45  Slackia        
#>  6 Otu_82  Barnesiella    
#>  7 Otu_78  Bacteroides    
#>  8 Otu_88  Odoribacter    
#>  9 Otu_259 Coprococcus    
#> 10 Otu_35  Collinsella    
#> # ℹ 13 more rows

## By default cool extracts the results of the first bake. If we have more
## bakes we can extract the one that you want with the bake parameter.
rec <- bake(rec, count_cutoff = 1)
cool(rec, 2)
#> ℹ Bake for count_cutoff = 1
#> # A tibble: 88 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_96  Prevotella     
#>  2 Otu_102 Prevotella     
#>  3 Otu_340 Holdemanella   
#>  4 Otu_332 Catenibacterium
#>  5 Otu_108 Prevotella     
#>  6 Otu_115 Alistipes      
#>  7 Otu_119 Alistipes      
#>  8 Otu_45  Slackia        
#>  9 Otu_82  Barnesiella    
#> 10 Otu_113 Prevotella     
#> # ℹ 78 more rows

## bake and cool methods needs a prep-Recipe. If you pass a non-PrepRecipe
## the output is an error.
data(test_rec)
err <- testthat::expect_error(cool(test_rec))
err
#> <error/rlang_error>
#> Error in `cool()`:
#> ! This function needs a PrepRecipe!
#> • Run prep(rec) and then try with cool()
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
#>  33.                             ├─testthat::expect_error(cool(test_rec))
#>  34.                             │ └─testthat:::expect_condition_matching(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             ├─dar::cool(test_rec)
#>  40.                             └─dar::cool(test_rec)
```
