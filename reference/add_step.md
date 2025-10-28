# Add a New Operation to the Current Recipe

`add_step` adds a step to the last location in the Recipe. `add_check`
does the same for checks.

## Usage

``` r
add_step(rec, object)

# S4 method for class 'Recipe'
add_step(rec, object)

# S4 method for class 'PrepRecipe'
add_step(rec, object)
```

## Arguments

- rec:

  A
  [`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md).

- object:

  A step or check object.

## Value

A updated
[`Recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
with the new operation in the last slot.
