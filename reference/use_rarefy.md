# Perform Rarefaction on Phyloseq Object

This function performs rarefaction on a phyloseq object if the `rarefy`
parameter is set to TRUE. Rarefaction is a process that randomly
subsamples the data to a specified depth. This is done to account for
differences in sequencing depth between samples. However, this process
is not without controversy. Rarefaction can lead to loss of information
and can also lead to false positives in differential abundance testing.
For more information, see
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2

## Usage

``` r
use_rarefy(phy, rarefy)
```

## Arguments

- phy:

  A phyloseq object.

- rarefy:

  A logical value indicating whether to perform rarefaction. If
  'no_seed', rarefaction is performed without a set seed. If FALSE, no
  rarefaction is performed.

## Value

A phyloseq object after rarefaction if `rarefy` is TRUE or "no_seed",
otherwise the original phyloseq object is returned.

## Examples

``` r
data(metaHIV_phy)

## With seed
# phy_rarefied <- dar:::use_rarefy(metaHIV_phy, TRUE)

## Witout seed
# phy_rarefied <- dar:::use_rarefy(metaHIV_phy, "no_seed")
```
