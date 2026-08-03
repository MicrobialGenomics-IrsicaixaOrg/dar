# Function steps_ids() @ L426

    Code
      print(test_prep_rec)
    Output
      -- DAR Results -----------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 101 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Genus 
      
      Results:
      
           v maaslin__Welsh_cake diff_taxa = 41 
           v deseq__Coussin_de_Lyon diff_taxa = 27 
      
           i 16 taxa are present in all tested methods 
      

# Function import_steps() @ L604

    Code
      suppressWarnings(import_steps(recipe(metaHIV_phy, "RiskGroup2", "Class"),
      system.file("extdata", "test_bake.json", package = "dar"), parallel = FALSE))
    Message
      ! bakes found in imported recipe
      i running `prep()`
    Output
      -- DAR Results -----------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 101 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Class 
      
      Results:
      
           v maaslin__Welsh_cake diff_taxa = 41 
           v deseq__Coussin_de_Lyon diff_taxa = 5 
      
           i 4 taxa are present in all tested methods 
      
      Bakes:
      
           (*) 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Carac 

