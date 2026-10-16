# PrepRecipe for metaHIV_phy data

A prepared version of
[test_rec](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/test_rec.md),
fitted through the centralized model API at genus level. Results retain
planned pairwise `RiskGroup2` contrasts and method-specific effects.
Original stored features remain available through
[`otu_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md);
[`tidy_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tidy_abundance.md)
returns their shared genus-level view. Unavailable native estimates
remain missing;
[`effect_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/effect_plt.md)
omits them rather than treating them as zero effects. Generated
reproducibly by `data-raw/metaHIV_recipes.R`.

## Usage

``` r
data("test_prep_rec")
```

## Format

A PrepRecipe object.
