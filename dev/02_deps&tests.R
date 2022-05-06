#########################################
####      TRACK OF DEV CHANGES      #####
#########################################

## DEPENDENCIES
usethis::use_pipe()
usethis::use_tidy_eval()

## TESTS
usethis::use_test("test-true")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "metaHIV_phy_data", open = FALSE)
