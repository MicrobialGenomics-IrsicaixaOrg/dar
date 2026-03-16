# Mutual finding plot

Plots number of differentially abundant features mutually found by
defined number of methods, colored by the differential abundance
direction and separated by comparison.

## Usage

``` r
mutual_plt(
  rec,
  count_cutoff = NULL,
  comparisons = NULL,
  steps = steps_ids(rec, type = "da"),
  top_n = 20
)
```

## Arguments

- rec:

  A prepped `Recipe` object.

- count_cutoff:

  Indicates the minimum number of methods in which an OTU must be
  present (Default: NULL). If count_cutoff is NULL count_cutoff is equal
  to `length(steps_ids(rec, "da")) * 2 / 3`.

- comparisons:

  By default, this function plots all comparisons. However, if the user
  indicates the comparison or comparisons of interest, only the selected
  ones will be plotted.

- steps:

  Character vector with step_ids to take in account. Default all "da"
  methods.

- top_n:

  Maximum number of taxa to represent. Default: 20.

## Value

ggplot2

## Examples

``` r
data(test_prep_rec)

## Running the function returns a tile plot,
mutual_plt(test_prep_rec)
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## The count_cutoff indicates the minimum number of methods in which an OTU
## must be present. By default the value is equal to
## length(steps_ids(rec, "da")) * 2 / 3 but it is customizable.
mutual_plt(
  test_prep_rec, 
  count_cutoff = length(steps_ids(test_prep_rec, "da"))
)
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## A single comparisons can be plotted through the comparison parameter.
mutual_plt(test_prep_rec, comparisons = c("hts_msm"))
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## If you want to exclude a method for the plot, you can remove it with the
## step parameter. In the following example we eliminate from the graph the
## results of maaslin.
mutual_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## mutual_plt function needs a PrepRecipe. If you pass a a non-PrepRecipe
## the output is an error.
data(test_rec)
err <- testthat::expect_error(mutual_plt(test_rec))
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
#>  33.                             ├─testthat::expect_error(mutual_plt(test_rec))
#>  34.                             │ └─testthat:::expect_condition_matching_(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             └─dar::mutual_plt(test_rec)
#>  40.                               └─dar:::check_prep_recipe(rec)
```
