# Overlap of significant OTUs between tested methods.

Overlap of significant OTUs between tested methods.

## Usage

``` r
overlap_df(rec, steps = steps_ids(rec, "da"), type = "all")

# S4 method for class 'Recipe'
overlap_df(rec, steps = steps_ids(rec, "da"), type = "all")

# S4 method for class 'PrepRecipe'
overlap_df(rec, steps = steps_ids(rec, "da"), type = "all")
```

## Arguments

- rec:

  A `Recipe` object.

- steps:

  Character vector with step_ids to take in account.

- type:

  Indicates whether to use all taxa ("all") or only those that are
  differentially abundant in at least one method ("da"). Default as
  "all".

## Value

df

## Examples

``` r
data(test_prep_rec)

## Running the function returns a UpSet plot ordered by frequency.
df <- overlap_df(test_prep_rec, steps_ids(test_prep_rec, "da"))
head(df)
#>                      maaslin__Eccles_cake deseq__Belekoy
#> maaslin__Eccles_cake            1.0000000      0.6435644
#> deseq__Belekoy                  0.6435644      1.0000000

## If you want to exclude a method for the plot, you can remove it with the
## step parameter. In the following example we eliminate from the graph the
## results of maaslin
overlap_df(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#>                deseq__Belekoy
#> deseq__Belekoy              1

## overlap_df function needs a prep-Recipe. If you pass a a non-prep
## Recipe the output is an error.
data(test_rec)
err <- testthat::expect_error(overlap_df(test_rec))
err
#> <error/rlang_error>
#> Error in `overlap_df()`:
#> ! This function needs a PrepRecipe!
#> • Run prep(rec) and then try with overlap_df()
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
#>  33.                             ├─testthat::expect_error(overlap_df(test_rec))
#>  34.                             │ └─testthat:::expect_condition_matching_(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             ├─dar::overlap_df(test_rec)
#>  40.                             └─dar::overlap_df(test_rec)
```
