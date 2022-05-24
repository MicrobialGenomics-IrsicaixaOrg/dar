#' Make a random identification field for steps
#'
#' @export
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#' @keywords internal
rand_id <- function(prefix = "step", len = 5) {
  candidates <- c(letters, LETTERS, paste(0:9))
  paste(prefix,
        paste0(sample(candidates, len, replace = TRUE), collapse = ""),
        sep = "_"
  )
}

#' Generate all unique contrasts between levels of a categorical variable.
#'
#' @param var categorical variable
#' @param phy phyloseq object
#' @param as_list boolean indicating if output must be returned as a list.
#' @param n_cut minimum of observations by level.
#'
#' @return tibble or list
get_comparisons <- function(var, phy, as_list = TRUE, n_cut = 1) {
  dat <-
    phyloseq::sample_data(phy) %>%
    to_tibble("sample_id") %>%
    dplyr::count(!!dplyr::sym(var)) %>%
    dplyr::filter(.data$n >= .env$n_cut) %>%
    dplyr::pull(!!dplyr::sym(var)) %>%
    as.character() %>%
    sort() %>%
    utils::combn(m = 2) %>%
    t() %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    stats::setNames(c("x", "y"))

  if (as_list) {
    dat <- purrr::map(seq_len(nrow(dat)), ~ c(dat[[.x, 1]], dat[[.x, 2]]))
  }
  dat
}

to_tibble <- function(df, id_name = "otu_id") {
  df %>%
    data.frame() %>%
    tibble::as_tibble(rownames = id_name)
}


step_to_expr <- function(step) {
  params <-
    step %>%
    purrr::discard(names(.) == "id") %>%
    purrr::map2_chr(names(.), ~ {
      if (is.null(.x)) { return(  glue::glue("{.y} = NULL"))}
      if (is.character(.x)) { return(glue::glue("{.y} = '{.x}'")) }
      if (inherits(.x, "formula")) { return(paste0(.y, " = ", paste0(.x, collapse = ""))) }
      glue::glue("{.y} = {.x}")
    }) %>%
    stringr::str_c(collapse = ", ")

  method <-
    step["id"] %>%
    stringr::str_remove_all(".{6}$")

  glue::glue("rec %>% dar::run_{method}({params})")
}


#' Title
#'
#' @param rec
#' @param steps
#'
#' @return
#' @export
#'
#' @examples
find_intersections <- function(rec, steps = steps_ids(rec, "da")) {
  intersection_df(rec, steps) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = -1) %>%
    dplyr::filter(value == 1) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(
      step_ids = purrr::map_chr(name, ~ .x) %>% stringr::str_c(collapse = ", "),
      sum_methods = sum(value)
    ) %>%
    dplyr::right_join(tax_table(rec), ., by = "taxa_id")
}

#' Title
#'
#' @param rec
#' @param steps
#' @param ordered_by
#'
#' @return
#' @export
#'
#' @examples
intersection_plt <- function(rec,
                             steps = steps_ids(rec, "da"),
                             ordered_by = c("freq", "degree")) {
  UpSetR::upset(
    data = intersection_df(rec, steps),
    nsets = length(rec@results),
    sets.bar.color = "#56B4E9",
    order.by = ordered_by
  )
}

#' Title
#'
#' @param rec
#' @param type
#'
#' @return
#' @export
#'
#' @examples
steps_ids <- function(rec, type = "all") {
  if (!type %in% c("all", "da", "prepro")) {
    rlang::abort(c(x = "asdf"))
  }

  out <- purrr::map_chr(rec@steps, ~ .x[["id"]])
  switch(
    type,
    "all" = out,
    "da" = purrr::discard(out, stringr::str_detect(out, "subset|filter")),
    "prepro" = purrr::keep(out, stringr::str_detect(out, "subset|filter"))
  )
}

#' Title
#'
#' @param rec
#' @param steps
#'
#' @return
#' @export
#'
#' @examples
intersection_df <- function(rec, steps = steps_ids(rec, "da")) {
  names(rec@results) %>%
    purrr::keep(. %in% steps) %>%
    purrr::set_names() %>%
    purrr::map_dfc( ~ {
      taxa <- rec@results[[.x]][[1]][["taxa_id"]]
      rownames(rec@phyloseq@otu_table) %>%
        tibble::tibble(taxa_id = .) %>%
        dplyr::mutate(!!.x := dplyr::if_else(taxa_id %in% taxa, 1, 0)) %>%
        dplyr::select(!!.x)
    }) %>%
    dplyr::mutate(taxa_id = rownames(rec@phyloseq@otu_table), .before = 1) %>%
    as.data.frame()
}

