# File created by roxut; edit the function definition file, not this file

# Test found in misc.R:217 (file:line)
  
data(test_prep_rec)
rec <- test_prep_rec
expect_equal(
  steps_ids(rec), 
  c("subset_taxa__Semla",
    "filter_taxa__Danish_pastry", 
    "maaslin__Tortell", 
    "metagenomeseq__Pogača", 
    "deseq__Tortita_negra"    
   )
)
expect_equal(
  steps_ids(rec, "da"), 
  c("maaslin__Tortell", "metagenomeseq__Pogača", "deseq__Tortita_negra")
)
expect_equal(
  steps_ids(rec, "prepro"), 
  c("subset_taxa__Semla", "filter_taxa__Danish_pastry")
)
expect_error(steps_ids(rec, "das"))
expect_type(steps_ids(rec), "character")
