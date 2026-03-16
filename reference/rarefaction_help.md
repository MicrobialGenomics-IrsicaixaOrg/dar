# Information about the Rarefaction process

Prints detailed information regarding the use of rarefaction in
microbiome analysis, its implications, and relevant literature.

## Usage

``` r
rarefaction_help()
```

## Value

invisible

## Examples

``` r
rarefaction_help()
#> ── Rarefaction in Microbiome Analysis ──────────────────────────────────────────
#> 
#> ℹ What is Rarefaction? It is a process that randomly subsamples the data to a specified depth to account for differences in sequencing depth between samples.
#> 
#> ! The Controversy:
#> ℹ It can lead to an irreversible loss of information.
#> ℹ It can lead to false positives in differential abundance testing.
#> 
#> ✔ Reproducibility: If performed with a set seed, it ensures reproducible results, but removes randomness in the subsampling process.
#> ℹ More info: <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2>
#> 
```
