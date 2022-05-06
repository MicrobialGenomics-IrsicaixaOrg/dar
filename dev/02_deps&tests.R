#########################################
####      TRACK OF DEV CHANGES      #####
#########################################

## DEPENDENCIES
usethis::use_pipe()
usethis::use_tidy_eval()
usethis::use_package("phyloseq")
usethis::use_package("data.table")
usethis::use_package("dplyr")
usethis::use_package("glue")
usethis::use_package("methods")
usethis::use_package("purrr")
usethis::use_package("stringr")
usethis::use_package("tibble")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "metaHIV_phy_data", open = FALSE)
