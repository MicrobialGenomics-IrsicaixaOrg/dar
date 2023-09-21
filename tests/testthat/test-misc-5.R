# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:209 (file:line)
  
data(test_prep_rec)
rec <- test_prep_rec
expect_equal(
  steps_ids(rec), 
  c("subset_taxa__Jalebi",
    "filter_taxa__Palmier", 
    "maaslin__Coussin_de_Lyon", 
    "metagenomeseq__Nazook", 
    "deseq__Pan_dulce"    
   )
)
expect_equal(
  steps_ids(rec, "da"), 
  c("maaslin__Coussin_de_Lyon", "metagenomeseq__Nazook", "deseq__Pan_dulce")
)
expect_equal(
  steps_ids(rec, "prepro"), 
  c("subset_taxa__Jalebi", "filter_taxa__Palmier")
)
expect_error(steps_ids(rec, "das"))
expect_type(steps_ids(rec), "character")
