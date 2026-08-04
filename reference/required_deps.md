# Returns required pakcages for Recipe object

Returns required pakcages for Recipe object

## Usage

``` r
required_deps(rec, steps = NULL)
```

## Arguments

- rec:

  A `Recipe` object

- steps:

  Optional list of steps whose dependencies should be checked.

## Value

character

## Examples

``` r
data(test_rec)

## The function returns instructions to install any uninstalled dependencies 
## needed to run the Recipe steps
dar:::required_deps(test_rec)

## The function also works with PrepRecipe-class objects
data(test_prep_rec)
dar:::required_deps(test_prep_rec)
```
