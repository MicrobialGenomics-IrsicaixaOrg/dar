# Plot the number of shared DA OTUs between methods.

Plot the number of shared DA OTUs between methods.

## Usage

``` r
exclusion_plt(
  rec,
  steps = steps_ids(rec, "da"),
  target = NULL,
  contrast_id = NULL
)
```

## Arguments

- rec:

  A prepped `Recipe` object.

- steps:

  Character vector with step_ids to take in account.

- target:

  Optional modeled target to include.

- contrast_id:

  Optional modeled contrast identifier to include.

## Value

ggplot2-class object

## Examples

``` r
data(test_prep_rec)

exclusion_plt(test_prep_rec)


## Exclude the first method from the plot
exclusion_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
```
