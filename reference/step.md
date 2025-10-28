# Overall Wrappers to Make New `step_X` or `check_Y` Objects

`step` sets the class of the `step` and `check` is for checks.

## Usage

``` r
step(subclass, ..., .prefix = "step_")

check(subclass, ..., .prefix = "check_")
```

## Arguments

- subclass:

  A character string for the resulting class. For example, if
  `subclass = "blah"` the step object that is returned has class
  `step_blah` or `check_blah` depending on the context.

- ...:

  All arguments to the operator that should be returned.

- .prefix:

  Prefix to the subclass created.

## Value

An updated step or check with the new class.
