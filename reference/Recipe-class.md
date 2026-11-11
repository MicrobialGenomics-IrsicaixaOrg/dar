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

- `steps`:

  List of step-class objects that will be used by DA.

- `model`:

  Optional centralized statistical model specification created by
  [`add_model()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_model.md).
