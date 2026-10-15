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
  top_n = 20,
  target = NULL,
  contrast_id = NULL
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

- target:

  Optional modeled target to include.

- contrast_id:

  Optional modeled contrast identifier to include.

## Value

ggplot2, or invisible `NULL` with an informational message when no
significant taxa are available for the selected hypotheses.

## Examples

``` r
data(test_prep_rec)

mutual_plt(test_prep_rec)
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## Require agreement among all methods
mutual_plt(
  test_prep_rec,
  count_cutoff = length(steps_ids(test_prep_rec, "da"))
)


## Show a single comparison
mutual_plt(test_prep_rec, comparisons = "hts_msm")
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.


## Exclude the first method from the plot
mutual_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#> ℹ `count_cutoff` set to 1.
#> ! Taxa present in selected methods are greater than the cutoff `top_n` = 20.
#> ℹ The top 20 significant taxa with the greatest overlap between methods will be
#>   used.
```
