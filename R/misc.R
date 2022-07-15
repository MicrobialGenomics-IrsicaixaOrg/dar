#' Make a random identification field for steps
#'
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#'
#' @export
#' @return character vector
#' @examples
#' rand_id("step")
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
#' @keywords internal
#' @examples \dontrun{
#' get_comparisons("RiskGroup2", metaHIV_phy)
#' }
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

#' Wrapper to convert phyloseq slots to tibble
#'
#' @param df output of `otu_table()`, `sample_data()` or `tax_table()` phyloseq functions.
#' @param id_name Name of the new column generated from rownames
#'
#' @return tibble
#' @keywords internal
#' @examples \dontrun{
#' to_tibble(phyloseq::otu_table(metaHIV_phy))
#' }
to_tibble <- function(df, id_name = "otu_id") {
  df %>%
    data.frame() %>%
    tibble::as_tibble(rownames = id_name)
}

#' Extracts parameters from steps and makes a character vector with the expression to
#' evaluate
#'
#' @param step object of class step
#'
#' @return character vector
#' @keywords internal
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

  glue::glue("rec %>% dar:::run_{method}({params})")
}


#' Finds common OTU between method results
#'
#' @param rec A recipe object.
#' @param steps character vector with step ids to take in account
#'
#' @return tibble
#' @export
find_intersections <- function(rec, steps = steps_ids(rec, "da")) {
  intersection_df(rec, steps) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = -1) %>%
    dplyr::filter(.data$value == 1) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(
      step_ids = purrr::map_chr(.data$name, ~ .x) %>% stringr::str_c(collapse = ", "),
      sum_methods = sum(.data$value)
    ) %>%
    dplyr::right_join(tax_table(rec), ., by = "taxa_id") %>%
    dplyr::arrange(-.data$sum_methods)
}

#' Get step_ids from recipe
#'
#' @param rec A recipe object.
#' @param type character vector indicating the type class. Options `c("all", "da",
#'   "prepro")`.
#'
#' @return character vector
#' @export
steps_ids <- function(rec, type = "all") {
  if (!type %in% c("all", "da", "prepro")) {
    rlang::abort(c(
      "Incorrect step type!",
      i = glue::glue(
        "Please use one of: {crayon::bgMagenta('c(all, da, prepro)')}"
      )
    ))
  }

  out <- purrr::map_chr(rec@steps, ~ .x[["id"]])
  switch(
    type,
    "all" = out,
    "da" = purrr::discard(out, stringr::str_detect(out, "subset|filter")),
    "prepro" = purrr::keep(out, stringr::str_detect(out, "subset|filter"))
  )
}

#' @keywords internal
tick <- function() {
  crayon::green(cli::symbol$tick)
}

#' @keywords internal
cross <- function() {
  crayon::red(cli::symbol$cross)
}

#' @keywords internal
info <- function() {
  crayon::blue(cli::symbol$info)
}

#' @keywords internal
dot <- function() {
  crayon::blue(cli::symbol$circle_filled)
}

#' Export step parameters as json.
#'
#' @param rec A recipe object.
#' @param file_name The path and file name of the optout file.
#'
#' @importFrom glue double_quote
#' @export
export_steps <- function(rec, file_name) {
  to_cat <-
    rec@steps %>%
    purrr::map_chr(~ {
      params <-
        names(.x) %>%
        purrr::map_chr(function(.y) {
          msg <- .x[[.y]]
          if (is.character(.x[[.y]]) | is.factor(.x[[.y]])) {
            msg <- double_quote(.x[[.y]])
          }
          stringr::str_c("   ", double_quote(.y), ": ", paste0(msg, collapse = ""), ",")
        }) %>%
        stringr::str_c(collapse = "\n")

      stringr::str_c("{\n", params, "\n}")
    })

  cat(to_cat, file = file_name)
}

#' Import steps from json file
#'
#' @param rec A recipe object.
#' @param file Path to the input file
#'
#' @export
import_steps <- function(rec, file) {
  lines <-
    readr::read_lines(file) %>%
    purrr::discard(stringr::str_detect(., "[{]|[}]"))

  id_idx <-  which(stringr::str_detect(lines, "id"))
  for (i in seq_along(id_idx)) {
    low_idx <- 1
    if (i != 1) {
      low_idx <- id_idx[i - 1] + 1
    }

    fun_name <-
      lines[id_idx[i]] %>%
      stringr::str_remove_all(".*: |,|\\\"") %>%
      stringr::str_remove_all(".{6}$")

    lines[low_idx:id_idx[i]] %>%
      purrr::map_chr(function(.x) {
        params <-
          .x %>%
          stringr::str_squish() %>%
          stringr::str_split(pattern = ": ") %>%
          unlist()

        param <-
          params[[1]] %>%
          stringr::str_remove_all("\\\"")

        value <-
          params[[2]] %>%
          stringr::str_squish() %>%
          stringr::str_remove_all(",$")

        if (stringr::str_detect(value, "%in%")) {
          value <-
            stringr::str_remove_all(value, "\\\"") %>%
            stringr::str_replace_all("\\\\", "\\\"")

          value <- encodeString(value, quote = "'")
        }

        if (stringr::str_count(value) == 0) {
          value = expression(NULL)
        }

        stringr::str_c(param, " = ", paste0(value, sep = ""))

      }) %>% stringr::str_c(collapse = ", ") %>%
      stringr::str_c("rec <<- step_", fun_name, "(rec, ", ., ")") %>%
      parse(text = .) %>%
      eval()
  }
  rec
}

