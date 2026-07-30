# Function add_step() @ L81

    Code
      rec <- step_maaslin(rec, id = "test_2")
    Message
      ! This step is already defined with the same parameters and will be skipped.
      i `rec %>% step_maaslin(min_abundance = 0, min_prevalence = 0.1, min_variance = 0, normalization = c('TSS'), transform = c('LOG'), max_significance = 0.1, random_effects = NULL, correction = c('BH'), standardize = TRUE, reference = NULL, median_comparison_abundance = TRUE, rarefy = FALSE)`

---

    Code
      prepped_dup <- bake(prepped)
    Message
      ! This step is already defined with the same parameters and will be skipped.
      i `rec %>% step_bake(count_cutoff = NULL, weights = NULL, exclude = NULL)`

