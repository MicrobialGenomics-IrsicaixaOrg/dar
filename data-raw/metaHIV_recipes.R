library(dar)
data("metaHIV_phy")

## metaHIV_prep_recipe
test_rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Genus") %>%
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) %>%
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") %>%
  step_maaslin() %>%
  step_metagenomeseq(rm_zeros = 0.01) %>%
  step_deseq()

test_rec

test_prep_rec <- prep(test_rec, parallel = TRUE)
test_prep_rec

usethis::use_data(test_rec, test_prep_rec, overwrite = TRUE)

export_steps(test_prep_rec, "inst/extdata/test.json")

rec <- bake(test_prep_rec)
export_steps(rec, "inst/extdata/test_bake.json")
