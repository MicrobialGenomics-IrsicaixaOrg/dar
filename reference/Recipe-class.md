# Recipe-class object

A Recipe is a description of the steps to be applied to a data set in
order to prepare it for data analysis.

## Usage

``` r
# S4 method for class 'PrepRecipe'
show(object)
```

## Arguments

- object:

  A Recipe object.

## Value

Recipe-class object

## Slots

- `phyloseq`:

  Phyloseq-class object.

- `var_info`:

  A tibble that contains the current set of terms in the data set. This
  initially defaults to the same data contained in `var_info`.

- `tax_info`:

  A tibble that contains the current set of taxonomic levels that will
  be used in the analysis.

- `steps`:

  List of step-class objects that will be used by DA.
