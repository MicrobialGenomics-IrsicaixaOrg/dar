# File created by roxut; edit the function definition file, not this file

# Test found in bake.R:20 (file:line)
  
data(test_prep_rec)

## Test internal check
rec <- test_prep_rec
weights <- c(2, 1)
names(weights) <- steps_ids(rec, "da")[1:2]
res <- bake(rec, weights = weights)
expect_error(cool(res))
