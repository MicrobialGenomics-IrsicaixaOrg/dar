# Define consensus strategies from a Recipe

For a prep Recipe adds a consensus strategy to use for result
extraction.

## Usage

``` r
bake(
  rec,
  count_cutoff = NULL,
  weights = NULL,
  exclude = NULL,
  id = rand_id("bake")
)
```

## Arguments

- rec:

  A PrepRecipe object. The step will be added to the sequence of
  operations for this Recipe.

- count_cutoff:

  Indicates the minimum number of methods in which an OTU must be
  present (Default: NULL). If `count_cutoff` is `NULL`, it is set to
  `length(steps_ids(rec, "da")) - length(exclude)`.

- weights:

  Named numeric vector with the weight for each method.

- exclude:

  Method ids to exclude.

- id:

  A character string that is unique to this step to identify it.

## Value

An object of class `PrepRecipe`

## Examples

``` r
data(test_prep_rec)
rec <- test_prep_rec

## Default bake extracts common OTUs in all DA tested methods
## (In this case the Recipe contains 3 methods)
res <- bake(rec)
cool(res)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 16 × 2
#>    taxa_id taxa            
#>    <chr>   <chr>           
#>  1 Otu_78  Bacteroides     
#>  2 Otu_88  Odoribacter     
#>  3 Otu_119 Alistipes       
#>  4 Otu_125 Parabacteroides 
#>  5 Otu_96  Prevotella      
#>  6 Otu_82  Barnesiella     
#>  7 Otu_94  Paraprevotella  
#>  8 Otu_332 Catenibacterium 
#>  9 Otu_255 Blautia         
#> 10 Otu_35  Collinsella     
#> 11 Otu_262 Dorea           
#> 12 Otu_365 Mitsuokella     
#> 13 Otu_259 Coprococcus     
#> 14 Otu_330 Ruthenibacterium
#> 15 Otu_307 Faecalibacterium
#> 16 Otu_433 Haemophilus     

## bake() requires a PrepRecipe
data(test_rec)
err <- testthat::expect_error(bake(test_rec))
err
#> <error/dar_error_invalid_input>
#> Error in `check_prep_recipe()`:
#> ✖ The argument `rec` must be a prepared <PrepRecipe>.
#> ! You supplied an unprepped <Recipe>.
#> ℹ Please run `prep()` on your recipe before using this function.
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
#>  33.                             ├─testthat::expect_error(bake(test_rec))
#>  34.                             │ └─testthat:::expect_condition_matching_(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             └─dar::bake(test_rec)
#>  40.                               └─dar:::check_prep_recipe(rec)

## We can use the parameter `count_cutoff` to select those OTUs
## shared by at least two methods
res <- bake(rec, count_cutoff = 2)
cool(res)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 16 × 2
#>    taxa_id taxa            
#>    <chr>   <chr>           
#>  1 Otu_78  Bacteroides     
#>  2 Otu_88  Odoribacter     
#>  3 Otu_119 Alistipes       
#>  4 Otu_125 Parabacteroides 
#>  5 Otu_96  Prevotella      
#>  6 Otu_82  Barnesiella     
#>  7 Otu_94  Paraprevotella  
#>  8 Otu_332 Catenibacterium 
#>  9 Otu_255 Blautia         
#> 10 Otu_35  Collinsella     
#> 11 Otu_262 Dorea           
#> 12 Otu_365 Mitsuokella     
#> 13 Otu_259 Coprococcus     
#> 14 Otu_330 Ruthenibacterium
#> 15 Otu_307 Faecalibacterium
#> 16 Otu_433 Haemophilus     

## Furthermore, we can exclude methods from the consensus strategy
## with the `exclude` parameter.
res <- bake(rec, exclude = steps_ids(rec, "da")[1])
cool(res)
#> ℹ Baking with count_cutoff = 1
#> Excluding "maaslin__Welsh_cake".
#> # A tibble: 27 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_34  Olsenella      
#>  2 Otu_35  Collinsella    
#>  3 Otu_45  Slackia        
#>  4 Otu_78  Bacteroides    
#>  5 Otu_82  Barnesiella    
#>  6 Otu_88  Odoribacter    
#>  7 Otu_94  Paraprevotella 
#>  8 Otu_96  Prevotella     
#>  9 Otu_119 Alistipes      
#> 10 Otu_125 Parabacteroides
#> # ℹ 17 more rows

## Finally, we can use the `weights` parameter to weigh each method.
weights <- c(2, 1)
names(weights) <- steps_ids(rec, "da")
res <- bake(rec, weights = weights)
cool(res)
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 41 × 2
#>    taxa_id taxa           
#>    <chr>   <chr>          
#>  1 Otu_78  Bacteroides    
#>  2 Otu_88  Odoribacter    
#>  3 Otu_119 Alistipes      
#>  4 Otu_129 Parabacteroides
#>  5 Otu_125 Parabacteroides
#>  6 Otu_96  Prevotella     
#>  7 Otu_82  Barnesiella    
#>  8 Otu_94  Paraprevotella 
#>  9 Otu_79  Bacteroides    
#> 10 Otu_51  Bacteroides    
#> # ℹ 31 more rows
```
