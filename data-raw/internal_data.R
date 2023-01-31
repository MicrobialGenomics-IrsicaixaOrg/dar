## code to prepare `pastry` dataset goes here

pastry_df <-
    "https://raw.githubusercontent.com/prasertcbs/basic-dataset/master/pastry.csv" %>%
    readr::read_csv(show_col_types = FALSE) %>%
    dplyr::mutate(id = stringr::str_remove_all(Name, "[(].*|[']")) %>%
    dplyr::mutate(id = stringr::str_replace_all(id, " ", "_"), .before = 1)

usethis::use_data(pastry_df, overwrite = TRUE, internal = TRUE)

