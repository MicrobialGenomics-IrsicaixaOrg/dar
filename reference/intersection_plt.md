# Plot results using UpSet plot

Plot results using UpSet plot

## Usage

``` r
intersection_plt(
  rec,
  steps = steps_ids(rec, "da"),
  ordered_by = c("freq", "degree"),
  font_size = 2
)

# S4 method for class 'Recipe'
intersection_plt(rec, steps, font_size)

# S4 method for class 'PrepRecipe'
intersection_plt(
  rec,
  steps = steps_ids(rec, "da"),
  ordered_by = c("freq", "degree"),
  font_size = 2
)
```

## Arguments

- rec:

  A `Recipe` object.

- steps:

  Character vector with step_ids to take in account.

- ordered_by:

  How the intersections in the matrix should be ordered by. Options
  include frequency (entered as "freq"), degree, or both in any order.

- font_size:

  Size of the font.

## Value

UpSet plot

## Examples

``` r
data(test_prep_rec)

## Running the function returns a UpSet plot ordered by frequency.
intersection_plt(test_prep_rec)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.


## Alternatively, you can order the plot by degree
intersection_plt(test_prep_rec, ordered_by = "degree")


## If you want to exclude a method for the plot, you can remove it with the
## step parameter. In the following example we eliminate from the graph the
## results of maaslin
intersection_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])


## intersection_plt function needs a PrepRecipe. If you pass a a non-prep
## Recipe the output is an error.
data(test_rec)
err <- testthat::expect_error(intersection_plt(test_rec))
err
#> <error/rlang_error>
#> Error in `.local()`:
#> ! This function needs a PrepRecipe!
#> • Run prep(rec) and then try with intersection_plt()
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
#>  33.                             ├─testthat::expect_error(intersection_plt(test_rec))
#>  34.                             │ └─testthat:::expect_condition_matching(...)
#>  35.                             │   └─testthat:::quasi_capture(...)
#>  36.                             │     ├─testthat (local) .capture(...)
#>  37.                             │     │ └─base::withCallingHandlers(...)
#>  38.                             │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#>  39.                             ├─dar::intersection_plt(test_rec)
#>  40.                             └─dar::intersection_plt(test_rec)
#>  41.                               └─dar (local) .local(rec, steps, font_size)
```
