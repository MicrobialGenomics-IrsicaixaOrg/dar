# Defunct functions in dar

These functions are retained only to provide migration guidance.
`phy_qc()` is defunct; use
[`recipe_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe_qc.md)
with an explicit `group_by` when grouped metrics are required.
`step_metagenomeseq()` was removed after being defunct for a complete
Bioconductor release cycle; use
[`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md)
or another supported differential-abundance engine.

## Usage

``` r
phy_qc(rec, ...)
```

## Arguments

- rec:

  A recipe object. This argument is no longer evaluated.

- ...:

  Additional arguments. These arguments are no longer evaluated.

## Value

This function does not return a value; it raises a defunct error.
