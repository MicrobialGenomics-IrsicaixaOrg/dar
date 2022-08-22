## code to prepare `pastry` dataset goes here

pastry_df <-
    "https://raw.githubusercontent.com/prasertcbs/basic-dataset/master/pastry.csv" %>%
    readr::read_csv(show_col_types = FALSE) %>%
    dplyr::mutate(id = stringr::str_replace_all(Name, " ", "_"), .before = 1)

usethis::use_data(pastry_df, overwrite = TRUE, internal = TRUE)
