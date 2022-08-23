library(dar)
data("metaHIV_phy")

## metaHIV_prep_recipe
test_rec <-
  recipe(metaHIV_phy, "RiskGroup2", "Genus") %>%
  step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
  step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") %>%
  step_maaslin() %>%
  step_metagenomeseq() %>%
  step_deseq()

test_prep_rec <- prep(test_rec, parallel = TRUE)

usethis::use_data(test_rec, test_prep_rec, overwrite = TRUE)
