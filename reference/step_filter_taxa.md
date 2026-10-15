# Filter taxa based on across-sample OTU abundance criteria

Add a preprocessing step that applies one safely serializable predicate
to each taxon's abundance vector across samples. Filtering is performed
during
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md);
taxa for which the predicate returns `FALSE` are removed.

## Usage

``` r
step_filter_taxa(rec, .f, id = rand_id("filter_taxa"))
```

## Arguments

- rec:

  A Recipe object. The step will be added to the sequence of operations
  for this Recipe.

- .f:

  A function with exactly one argument named `x`, without a default, or
  a character string containing that function expression. It must return
  a single logical value for the abundance vector `x`. Only the
  restricted expression syntax described below is accepted; lists of
  functions are not supported.

- id:

  A character string that is unique to this step to identify it.

## Value

A `Recipe` with the filtering step added. The constructor does not
filter the data or return a logical vector. During
[`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md),
the step uses
[`phyloseq::filter_taxa()`](https://rdrr.io/pkg/phyloseq/man/filter_taxa.html)
with `prune = TRUE`; there is no `prune` argument in this wrapper.

## Details

The predicate body must be a single expression, without braces,
assignments, external variables or user-defined helper calls. Write
thresholds as literal numbers inside the expression, even when using the
function directly in R. These restrictions keep recipes safely
exportable and importable; they are checked when the step is added, not
just on export.

Accepted expressions use `x`, finite numeric or non-missing logical
literals, parentheses, arithmetic operators (`+`, `-`, `*`, `/`, `^`,
`%%`, `%/%`), comparisons (`>`, `>=`, `<`, `<=`, `==`, `!=`), logical
operators (`&`, `|`, `!`), and calls to `sum`, `length`, `mean`,
`median`, `min`, `max`, `var`, `sd`, `any`, `all` or `is.na`. The only
supported named argument is `na.rm` in `sum`, `mean`, `median`, `min`,
`max`, `var`, `sd`, `any` and `all`. Namespace-qualified calls are not
accepted.

## See also

Other filter phy steps:
[`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md),
[`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md),
[`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md),
[`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md)

## Examples

``` r
data(metaHIV_phy)

## Keep taxa present in at least 3% of samples
recipe(metaHIV_phy) |>
  step_filter_taxa(.f = function(x) sum(x > 0) >= 0.03 * length(x)) |>
  prep()
#> ── DAR Results ─────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 279 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ✖ undefined taxonomic level. Use add_model() to define the analysis! 
#> Results:
#> 

## The same predicate can be supplied as a string
recipe(metaHIV_phy) |>
  step_filter_taxa(.f = "function(x) sum(x > 0) >= 0.03 * length(x)")
#> ── DAR Recipe ──────────────────────────────────────────────────────────────────
#> Inputs:
#> 
#>      ℹ phyloseq object with 451 taxa and 156 samples 
#>      ✖ undefined analysis target. Use add_model() to define the analysis! 
#>      ✖ undefined taxonomic level. Use add_model() to define the analysis! 
#> Preprocessing steps:
#> 
#>      ◉ step_filter_taxa() id = filter_taxa__Pithivier 
#> 
#> DA steps:
#> 
```
