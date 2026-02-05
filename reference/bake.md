# Define consensus strategies from a Recipe

For a prep Recipe adds a consensus strategies to use for result
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

# S4 method for class 'PrepRecipe'
bake(
  rec,
  count_cutoff = NULL,
  weights = NULL,
  exclude = NULL,
  id = rand_id("bake")
)

# S4 method for class 'Recipe'
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

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- count_cutoff:

  Indicates the minimum number of methods in which an OTU must be
  present (Default: NULL). If count_cutoff is NULL count_cutoff is equal
  to `length(steps_ids(rec, "da")) - length(exclude)`

- weights:

  Named vector with the ponderation value for each method.

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

## bake and cool methods needs a PrepRecipe. If you pass a non-PrepRecipe
## the output is an error.
data(test_rec)
err <- testthat::expect_error(bake(test_rec))
err
#> <error/rlang_error>
#> Error in `bake()`:
#> ! This function needs a PrepRecipe!
#> • Run prep(rec) and then try with bake()
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
#>  39.                             ├─dar::bake(test_rec)
#>  40.                             └─dar::bake(test_rec)

## We can use the parameter `cout_cutoff` to for example select those OTUs
## shared with at least two methods
res <- bake(rec, count_cutoff = 2)
cool(res)
#> ℹ Bake for count_cutoff = 2
#> # A tibble: 34 × 2
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
#> # ℹ 24 more rows

## Furthermore, we can exclude methods from the consensus strategy with the 
## `exclude` parameter.
res <- bake(rec, exclude = steps_ids(rec, "da")[1])
cool(res)
#> ℹ Bake for count_cutoff = 2
#> ℹ Results from maaslin__ChaSan are excluded
#> # A tibble: 34 × 2
#>    taxa_id taxa                 
#>    <chr>   <chr>                
#>  1 Otu_371 Megasphaera          
#>  2 Otu_361 Phascolarctobacterium
#>  3 Otu_392 Parasutterella       
#>  4 Otu_214 Clostridium          
#>  5 Otu_365 Mitsuokella          
#>  6 Otu_332 Catenibacterium      
#>  7 Otu_34  Olsenella            
#>  8 Otu_119 Alistipes            
#>  9 Otu_397 Bilophila            
#> 10 Otu_407 Anaerobiospirillum   
#> # ℹ 24 more rows

## Finally, we can use the `weights` parameter to weigh each method.
weights <- c(2, 1, 1)
names(weights) <- steps_ids(rec, "da")
res <- bake(rec, weights = weights)
cool(res)
#> ℹ Bake for count_cutoff = 3
#> # A tibble: 34 × 2
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
#> # ℹ 24 more rows
```
