# Function step() @ L36

    Code
      step_metagenomeseq(step_metagenomeseq(recipe(metaHIV_phy, "RiskGroup2",
        "Species"), id = "panallet"))
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_metagenomeseq(zeroMod = NULL, useCSSoffset = TRUE, useMixedModel = FALSE, max_significance = 0.05, log2FC = 0, rarefy = FALSE, rm_zeros = 0)
    Output
      -- DAR Recipe ------------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 451 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Species 
      
      Preporcessing steps:
      
      
      DA steps:
      
           (*) step_metagenomeseq() id = panallet 

---

    Code
      expect_s4_class(step_metagenomeseq(step_metagenomeseq(recipe(metaHIV_phy,
        "RiskGroup2", "Species"))), "Recipe")
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_metagenomeseq(zeroMod = NULL, useCSSoffset = TRUE, useMixedModel = FALSE, max_significance = 0.05, log2FC = 0, rarefy = FALSE, rm_zeros = 0)
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_metagenomeseq(zeroMod = NULL, useCSSoffset = TRUE, useMixedModel = FALSE, max_significance = 0.05, log2FC = 0, rarefy = FALSE, rm_zeros = 0)

---

    Code
      expect_s4_class(bake(bake(test_prep_rec)), "PrepRecipe")
    Message
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_bake(count_cutoff = NULL, weights = NULL, exclude = NULL)
      ! This step is already defined with the same parameters and will be skipped: 
      rec %>% step_bake(count_cutoff = NULL, weights = NULL, exclude = NULL)

