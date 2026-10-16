library(dar)
data("metaHIV_phy")
set.seed(130)

## Keep the configured IDs stable when rebuilding the bundled example.
test_rec <- recipe(metaHIV_phy) |>
  add_model(
    ~ RiskGroup2, targets = "RiskGroup2", tax_level = "Genus",
    reference = c(RiskGroup2 = "hts")
  ) |>
  step_subset_taxa(
    tax_level = "Kingdom", taxa = c("Bacteria", "Archaea"),
    id = "subset_taxa__Viennoiserie"
  ) |>
  step_filter_taxa(
    .f = "function(x) sum(x > 0) >= (0.3 * length(x))",
    id = "filter_taxa__Karakudamono"
  ) |>
  step_maaslin(id = "maaslin__Welsh_cake") |>
  step_deseq(id = "deseq__Coussin_de_Lyon")

test_rec

test_prep_rec <- prep(test_rec, parallel = FALSE)
test_prep_rec

## Validate the fitted units before replacing the example data.
analysis_taxa <- dar:::analysis_tax_table(test_prep_rec)
results <- tidy_results(test_prep_rec)
expected_labels <- analysis_taxa$taxa[
  match(results$taxa_id, analysis_taxa$taxa_id)
]
stopifnot(
  identical(get_model(test_prep_rec)$tax_level, "Genus"),
  all(results$taxa_id %in% analysis_taxa$taxa_id),
  !anyDuplicated(results[c("taxa_id", "contrast_id", "step_id")]),
  identical(results$taxa, expected_labels)
)

usethis::use_data(test_rec, test_prep_rec, overwrite = TRUE)

export_steps(test_prep_rec, "inst/extdata/test.json")

rec <- bake(test_prep_rec, id = "bake__bundled")
export_steps(rec, "inst/extdata/test_bake.json")
