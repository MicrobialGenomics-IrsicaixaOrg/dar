benchmark_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
repository_root <- normalizePath(file.path(benchmark_root, "..", ".."), mustWork = TRUE)
pkgload::load_all(repository_root, quiet = TRUE)
source(file.path(benchmark_root, "R", "runtime.R"), local = .GlobalEnv)
benchmark_load(benchmark_root)

# This is a result-contract fixture, not a fitted statistical method.
benchmark_prepared_fixture <- function(simulation, id, selected = TRUE, direction = 1) {
  rec <- benchmark_engine_recipe(simulation$recipe, id)
  tax <- dar:::analysis_tax_table(rec)[1L, , drop = FALSE]
  plan <- simulation$contrast_plan[1L, , drop = FALSE]
  raw <- tibble::tibble(
    taxa_id = tax$taxa_id, taxa = tax$taxa, contrast_id = plan$contrast_id,
    comparison = plan$comparison, contrast_type = plan$contrast_type, var = plan$var,
    effect = as.numeric(direction), signif = selected,
    padj = 0.01, pval = 0.005, p = 0.005, q = 0.01, pvalue = 0.005,
    we.ep = 0.005, we.eBH = 0.01, qval = 0.01
  )
  prepared <- dar:::prep_recipe(rec, stats::setNames(list(list(model = raw)), id), list(),
    execution = list(contrasts = plan, dropped_samples = character(), executed_steps = id,
      skipped_steps = tibble::tibble(step_id = character(), engine = character(),
                                    compatible = logical(), reason = character())))
  list(status = "success", prepared = prepared,
       results = as.data.frame(dar::tidy_results(prepared)))
}
