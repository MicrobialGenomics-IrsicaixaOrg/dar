# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:216 (file:line)
  
data(test_prep_rec)
rec <- test_prep_rec
expect_equal(
  steps_ids(rec), 
  c("subset_taxa__Boyoz",
    "filter_taxa__Charlotte", 
    "maaslin__ChaSan", 
    "metagenomeseq__Zlebia", 
    "deseq__Linzer_torte"    
   )
)
expect_equal(
  steps_ids(rec, "da"), 
  c("maaslin__ChaSan", "metagenomeseq__Zlebia", "deseq__Linzer_torte")
)
expect_equal(
  steps_ids(rec, "prepro"), 
  c("subset_taxa__Boyoz", "filter_taxa__Charlotte")
)
expect_error(steps_ids(rec, "das"))
expect_type(steps_ids(rec), "character")
