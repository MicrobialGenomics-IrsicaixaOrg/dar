# Function steps_ids() @ L252

    Code
      print(test_prep_rec)
    Output
      -- DAR Results -----------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 101 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Genus 
      
      Results:
      
           v maaslin__Eccles_cake diff_taxa = 41 
           v deseq__Belekoy diff_taxa = 27 
      
           i 16 taxa are present in all tested methods 
      

# Function export_steps() @ L330

    Code
      readr::read_lines(file)
    Output
       [1] "{"                                                               
       [2] "   \"tax_level\": [c(\"Kingdom\")],"                             
       [3] "   \"taxa\": [c(\"Bacteria\", \"Archaea\")],"                    
       [4] "   \"id\": [c(\"subset_taxa__Bear_claw\")],"                     
       [5] "}"                                                               
       [6] "{"                                                               
       [7] "   \".f\": [c(\"function(x) sum(x > 0) >= (0.3 * length(x))\")],"
       [8] "   \"id\": [c(\"filter_taxa__Spanakopita\")],"                   
       [9] "}"                                                               
      [10] "{"                                                               
      [11] "   \"min_abundance\": 0,"                                        
      [12] "   \"min_prevalence\": 0.1,"                                     
      [13] "   \"min_variance\": 0,"                                         
      [14] "   \"normalization\": [c(\"TSS\")],"                             
      [15] "   \"transform\": [c(\"LOG\")],"                                 
      [16] "   \"max_significance\": 0.1,"                                   
      [17] "   \"random_effects\": ,"                                        
      [18] "   \"correction\": [c(\"BH\")],"                                 
      [19] "   \"standardize\": TRUE,"                                       
      [20] "   \"reference\": ,"                                             
      [21] "   \"median_comparison_abundance\": TRUE,"                       
      [22] "   \"rarefy\": FALSE,"                                           
      [23] "   \"id\": [c(\"maaslin__Eccles_cake\")],"                       
      [24] "}"                                                               
      [25] "{"                                                               
      [26] "   \"test\": [c(\"Wald\")],"                                     
      [27] "   \"fitType\": [c(\"local\")],"                                 
      [28] "   \"betaPrior\": FALSE,"                                        
      [29] "   \"type\": [c(\"ashr\")],"                                     
      [30] "   \"max_significance\": 0.05,"                                  
      [31] "   \"log2FC\": 0,"                                               
      [32] "   \"rarefy\": FALSE,"                                           
      [33] "   \"id\": [c(\"deseq__Belekoy\")],"                             
      [34] "}"                                                               

# Function import_steps() @ L397

    Code
      import_steps(recipe(metaHIV_phy, "RiskGroup2", "Class"), system.file("extdata",
        "test_bake.json", package = "dar"))
    Message
      ! bakes found in imported recipe
      i running prep()
    Output
      -- DAR Results -----------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 101 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Class 
      
      Results:
      
           v maaslin__Eccles_cake diff_taxa = 41 
           v deseq__Belekoy diff_taxa = 5 
      
           i 4 taxa are present in all tested methods 
      
      Bakes:
      
           (*) 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Makmur 

