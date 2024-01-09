# Function steps_ids() @ L253

    Code
      print(test_prep_rec)
    Output
      -- DAR Results -----------------------------------------------------------------
      Inputs:
      
           i phyloseq object with 101 taxa and 156 samples 
           i variable of interes RiskGroup2 (class: character, levels: hts, msm, pwid) 
           i taxonomic level Genus 
      
      Results:
      
           v maaslin__ChaSan diff_taxa = 86 
           v metagenomeseq__Zlebia diff_taxa = 32 
           v deseq__Linzer_torte diff_taxa = 27 
      
           i 23 taxa are present in all tested methods 
      

# Function export_steps() @ L331

    Code
      readr::read_lines(file)
    Output
       [1] "{"                                                               
       [2] "   \"tax_level\": [c(\"Kingdom\")],"                             
       [3] "   \"taxa\": [c(\"Bacteria\", \"Archaea\")],"                    
       [4] "   \"id\": [c(\"subset_taxa__Boyoz\")],"                         
       [5] "}"                                                               
       [6] "{"                                                               
       [7] "   \".f\": [c(\"function(x) sum(x > 0) >= (0.3 * length(x))\")],"
       [8] "   \"id\": [c(\"filter_taxa__Charlotte\")],"                     
       [9] "}"                                                               
      [10] "{"                                                               
      [11] "   \"min_abundance\": 0,"                                        
      [12] "   \"min_prevalence\": 0.1,"                                     
      [13] "   \"min_variance\": 0,"                                         
      [14] "   \"normalization\": [c(\"TSS\")],"                             
      [15] "   \"transform\": [c(\"LOG\")],"                                 
      [16] "   \"analysis_method\": [c(\"LM\")],"                            
      [17] "   \"max_significance\": 0.25,"                                  
      [18] "   \"random_effects\": ,"                                        
      [19] "   \"correction\": [c(\"BH\")],"                                 
      [20] "   \"standardize\": TRUE,"                                       
      [21] "   \"reference\": ,"                                             
      [22] "   \"rarefy\": FALSE,"                                           
      [23] "   \"id\": [c(\"maaslin__ChaSan\")],"                            
      [24] "}"                                                               
      [25] "{"                                                               
      [26] "   \"zeroMod\": ,"                                               
      [27] "   \"useCSSoffset\": TRUE,"                                      
      [28] "   \"useMixedModel\": FALSE,"                                    
      [29] "   \"max_significance\": 0.05,"                                  
      [30] "   \"log2FC\": 0,"                                               
      [31] "   \"rarefy\": FALSE,"                                           
      [32] "   \"rm_zeros\": 0.01,"                                          
      [33] "   \"id\": [c(\"metagenomeseq__Zlebia\")],"                      
      [34] "}"                                                               
      [35] "{"                                                               
      [36] "   \"test\": [c(\"Wald\")],"                                     
      [37] "   \"fitType\": [c(\"local\")],"                                 
      [38] "   \"betaPrior\": FALSE,"                                        
      [39] "   \"type\": [c(\"ashr\")],"                                     
      [40] "   \"max_significance\": 0.05,"                                  
      [41] "   \"log2FC\": 0,"                                               
      [42] "   \"rarefy\": FALSE,"                                           
      [43] "   \"id\": [c(\"deseq__Linzer_torte\")],"                        
      [44] "}"                                                               

# Function import_steps() @ L398

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
      
           v maaslin__ChaSan diff_taxa = 86 
           v metagenomeseq__Zlebia diff_taxa = 6 
           v deseq__Linzer_torte diff_taxa = 5 
      
           i 4 taxa are present in all tested methods 
      
      Bakes:
      
           (*) 1 -> count_cutoff: NULL, weights: NULL, exclude: NULL, id: bake__Trdelník 

