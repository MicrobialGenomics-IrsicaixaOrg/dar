########################################
####       INIT REPO PACKAGE       #####
########################################

## Use git ----
usethis::use_git()

## Fill the DESCRIPTION ----
## Add meta data about your application
usethis::use_description()

## Create Common Files ----
## See ?usethis for more information
usethis::use_mit_license(copyright_holder = "Francesc Català-Moll")
usethis::use_readme_rmd(open = FALSE)
usethis::use_code_of_conduct(contact = "fcatala@iriscaixa.es")
usethis::use_lifecycle_badge("Experimental")
usethis::use_news_md(open = FALSE)
usethis::use_github_links()
usethis::use_testthat()
