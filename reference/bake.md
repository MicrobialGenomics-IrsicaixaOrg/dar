# Define consensus strategies from a Recipe

For a prep Recipe adds a consensus strategy to use for result
extraction.

## Usage

``` r
bake(
  rec,
  count_cutoff = NULL,
  weights = NULL,
  exclude = NULL,
  id = rand_id("bake")
)
```

## Arguments

- rec:

  A PrepRecipe object. The step will be added to the sequence of
  operations for this Recipe.

- count_cutoff:

  Indicates the minimum number of methods in which an OTU must be
  present (Default: NULL). If `count_cutoff` is `NULL`, it is set to
  `length(steps_ids(rec, "da")) - length(exclude)`.

- weights:

  Named numeric vector with the weight for each method.

- exclude:

  Method ids to exclude.

- id:

  A character string that is unique to this step to identify it.

## Value

An object of class `PrepRecipe`

## Examples

``` r
data(test_prep_rec)

## Require agreement among all differential-abundance methods
cool(bake(test_prep_rec))
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 16 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  2 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  3 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  4 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                2
#>  5 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… up                2
#>  6 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… down              2
#>  7 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… up                2
#>  8 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… down              2
#>  9 Otu_35  Colli… RiskGroup2… RiskGroup… main          Risk… up                2
#> 10 Otu_354 Firmi… RiskGroup2… RiskGroup… main          Risk… up                2
#> 11 Otu_78  Bacte… RiskGroup2… RiskGroup… main          Risk… down              2
#> 12 Otu_78  Bacte… RiskGroup2… RiskGroup… main          Risk… up                2
#> 13 Otu_82  Barne… RiskGroup2… RiskGroup… main          Risk… down              2
#> 14 Otu_88  Odori… RiskGroup2… RiskGroup… main          Risk… down              2
#> 15 Otu_88  Odori… RiskGroup2… RiskGroup… main          Risk… up                2
#> 16 Otu_96  Prevo… RiskGroup2… RiskGroup… main          Risk… up                2
#> # ℹ 1 more variable: methods <chr>

## Retain taxa detected by at least one method
cool(bake(test_prep_rec, count_cutoff = 1))
#> ℹ Baking with count_cutoff = 1
#> # A tibble: 46 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_1   Metha… RiskGroup2… RiskGroup… main          Risk… up                1
#>  2 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              2
#>  3 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                2
#>  4 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              2
#>  5 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… up                1
#>  6 Otu_231 Eubac… RiskGroup2… RiskGroup… main          Risk… up                1
#>  7 Otu_255 Blaut… RiskGroup2… RiskGroup… main          Risk… up                1
#>  8 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                2
#>  9 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… down              1
#> 10 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… up                2
#> # ℹ 36 more rows
#> # ℹ 1 more variable: methods <chr>

## Exclude a method from the consensus
method_ids <- steps_ids(test_prep_rec, "da")
cool(bake(test_prep_rec, exclude = method_ids[1]))
#> ℹ Baking with count_cutoff = 1
#> Excluding "maaslin__Welsh_cake".
#> # A tibble: 36 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_1   Metha… RiskGroup2… RiskGroup… main          Risk… up                1
#>  2 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              1
#>  3 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                1
#>  4 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              1
#>  5 Otu_231 Eubac… RiskGroup2… RiskGroup… main          Risk… up                1
#>  6 Otu_255 Blaut… RiskGroup2… RiskGroup… main          Risk… up                1
#>  7 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                1
#>  8 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… up                1
#>  9 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… down              1
#> 10 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… up                1
#> # ℹ 26 more rows
#> # ℹ 1 more variable: methods <chr>

## Give the first method twice the weight of the others
weights <- setNames(rep(1, length(method_ids)), method_ids)
weights[1] <- 2
cool(bake(test_prep_rec, weights = weights))
#> ℹ Baking with count_cutoff = 2
#> # A tibble: 26 × 9
#>    taxa_id taxa   contrast_id comparison contrast_type var   effect method_count
#>    <chr>   <chr>  <chr>       <chr>      <chr>         <chr> <chr>         <dbl>
#>  1 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… down              3
#>  2 Otu_119 Alist… RiskGroup2… RiskGroup… main          Risk… up                3
#>  3 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… down              3
#>  4 Otu_125 Parab… RiskGroup2… RiskGroup… main          Risk… up                2
#>  5 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… up                3
#>  6 Otu_259 Copro… RiskGroup2… RiskGroup… main          Risk… down              2
#>  7 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… up                3
#>  8 Otu_262 Dorea  RiskGroup2… RiskGroup… main          Risk… down              3
#>  9 Otu_298 Intes… RiskGroup2… RiskGroup… main          Risk… up                2
#> 10 Otu_332 Caten… RiskGroup2… RiskGroup… main          Risk… up                3
#> # ℹ 16 more rows
#> # ℹ 1 more variable: methods <chr>
```
