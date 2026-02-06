# Function step() @ L36

    Code
      step_maaslin(step_maaslin(recipe(metaHIV_phy, "RiskGroup2", "Species"), id = "panallet"))
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_maaslin(min_abundance = 0, min_prevalence = 0.1, min_variance = 0, normalization = c('TSS'), transform = c('LOG'), max_significance = 0.1, random_effects = NULL, correction = c('BH'), standardize = TRUE, reference = NULL, median_comparison_abundance = TRUE, rarefy = FALSE)
    Output
      -- DAR Recipe ------------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 451 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Species 
      
      Preporcessing steps:
      
      
      DA steps:
      
           (*) step_maaslin() id = panallet 

---

    Code
      expect_s4_class(step_maaslin(step_maaslin(recipe(metaHIV_phy, "RiskGroup2",
        "Species"))), "Recipe")
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_maaslin(min_abundance = 0, min_prevalence = 0.1, min_variance = 0, normalization = c('TSS'), transform = c('LOG'), max_significance = 0.1, random_effects = NULL, correction = c('BH'), standardize = TRUE, reference = NULL, median_comparison_abundance = TRUE, rarefy = FALSE)

---

    Code
      expect_s4_class(bake(bake(test_prep_rec)), "PrepRecipe")
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_bake(count_cutoff = NULL, weights = NULL, exclude = NULL)

