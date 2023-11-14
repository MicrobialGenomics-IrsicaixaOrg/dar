#' Make a random identification field for steps
#'
#' @param prefix A single character string
#' @return A character string with the prefix and random letters separated by
#'   and underscore.
#'
#' @export
#' @autoglobal
#' @return character vector
#' @tests testthat
#' set.seed(123)
#' expect_equal(rand_id(), "step__Filo")
#' 
#' @examples
#' rand_id("step")
rand_id <- function(prefix = "step") {
  candidate <- sample(x = pastry_df$id, size = 1)
  paste(prefix, candidate, sep = "__")
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
#' @autoglobal
#' 
#' @tests testthat
#' data("metaHIV_phy")
#' phy <- metaHIV_phy
#' var <- "RiskGroup2"
#' 
#' is_list <- get_comparisons(var, phy)
#' is_tibb <- get_comparisons(var, phy, as_list = FALSE)
#' filtered <- get_comparisons(var, phy, n_cut = 30)
#' 
#' expect_equal(class(is_list), "list")
#' expect_s3_class(is_tibb, "tbl_df")
#' expect_length(is_list, 3)
#' expect_equal(nrow(is_tibb), 3)
#' expect_length(filtered, 1)
#' 
#' @examples 
#' data(test_rec)
#' get_phy(test_rec) %>% 
#'   dar:::get_comparisons("RiskGroup2", .)
get_comparisons <- function(var, phy, as_list = TRUE, n_cut = 1) {
  dat <-
    phyloseq::sample_data(phy) %>%
    to_tibble("sample_id") %>%
    dplyr::count(!!dplyr::sym(var)) %>%
    dplyr::filter(n >= .env$n_cut) %>%
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
#' @param df output of `otu_table()`, `sample_data()` or `tax_table()` phyloseq
#'   functions.
#' @param id_name Name of the new column generated from rownames
#'
#' @return tibble
#' @keywords internal
#' @autoglobal
#' 
#' @tests testthat
#' data("metaHIV_phy")
#' phy <- metaHIV_phy
#' df <- phyloseq::otu_table(phy)
#' 
#' test_1 <- to_tibble(df)
#' test_2 <- to_tibble(df, id_name = "random_id")
#' 
#' expect_length(test_1, 157)
#' expect_s3_class(test_1, "tbl_df")
#' expect_s3_class(test_2, "tbl_df")
#' expect_equal(names(test_1)[1], "otu_id")
#' expect_equal(names(test_2)[1], "random_id")
#' 
#' @examples
#' data(test_rec)
#' otu_table <- 
#'   get_phy(test_rec) %>% 
#'   phyloseq::otu_table()
#' 
#' dar:::to_tibble(otu_table)
to_tibble <- function(df, id_name = "otu_id") {
  df %>%
    data.frame(check.names = FALSE) %>%
    tibble::as_tibble(rownames = id_name)
}

#' Extracts parameters from steps and makes a character vector with the
#' expression to evaluate
#'
#' @param step object of class step
#'
#' @return character vector
#' @keywords internal
#' @autoglobal
#' @tests testthat
#' data(test_prep_rec)
#' exprs <- test_prep_rec@steps %>% purrr::map_chr(step_to_expr)
#' expect_length(exprs, 5)
#' expect_true(all(stringr::str_detect(exprs, "dar:::run_")))
step_to_expr <- function(step) {
  params <-
    step %>%
    purrr::discard(names(.) == "id") %>%
    purrr::map2_chr(names(.), ~ {
      if (is.null(.x)) { return(  glue::glue("{.y} = NULL"))}
      if (is.character(.x)) { return(glue::glue("{.y} = '{.x}'")) }
      if (inherits(.x, "formula")) { 
        return(paste0(.y, " = ", paste0(.x, collapse = ""))) 
      }
      if (.y == "weights" & is(step, "step_bake") & !is.null(.x)) {
        text <- .x %>% 
          purrr::map2_chr(names(.), ~ { paste0(.y, " = ", paste0(.x)) }) %>% 
          stringr::str_c(collapse = ", ")
        return(glue::glue("{.y} = c({text})"))
      }
      
      glue::glue("{.y} = {.x}")
    }) %>%
    stringr::str_c(collapse = ", ")

  method <-
    step["id"] %>%
    stringr::str_remove_all("__.*")

  glue::glue("rec %>% dar:::run_{method}({params})")
}


#' Finds common OTU between method results
#'
#' @param rec A recipe object.
#' @param steps character vector with step ids to take in account
#'
#' @return tibble
#' @export
#' @autoglobal
#' 
#' @tests testthat
#' data(test_prep_rec)
#' res_1 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")
#' )
#' res_2 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-1]
#' )
#' res_3 <- find_intersections(
#'   test_prep_rec, steps = steps_ids(test_prep_rec, type = "da")[-2]
#' )
#' 
#' expect_equal(nrow(res_1), 88)
#' expect_equal(nrow(res_2), 36)
#' expect_equal(nrow(res_3), 86)
#' expect_s3_class(res_1, "tbl_df")
#' 
#' @examples 
#' data(test_prep_rec)
#' 
#' ## From a prep-recipe we can extract a tibble with all intersections
#' intersections <- find_intersections(test_prep_rec)
#' intersections
#' 
#' ## Additionally, we can exclude some methods form the table
#' intersections <- find_intersections(
#'   test_prep_rec, 
#'   steps = steps_ids(test_prep_rec, "da")[-1]
#' )
#' 
#' intersections
find_intersections <- function(rec, steps = steps_ids(rec, "da")) {
  intersection_df(rec, steps) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = -1) %>%
    dplyr::filter(value == 1) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(
      step_ids = 
        purrr::map_chr(name, ~ .x) %>% 
        stringr::str_c(collapse = ", "),
      sum_methods = sum(value)
    ) %>%
    dplyr::right_join(tax_table(rec), ., by = "taxa_id") %>%
    dplyr::arrange(-sum_methods)
}

#' Get step_ids from recipe
#'
#' @param rec A recipe object.
#' @param type character vector indicating the type class. Options `c("all",
#'   "da", "prepro")`.
#'
#' @return character vector
#' @export
#' @autoglobal
#' @tests testthat
#' data(test_prep_rec)
#' rec <- test_prep_rec
#' expect_equal(
#'   steps_ids(rec), 
#'   c("subset_taxa__Jalebi",
#'     "filter_taxa__Palmier", 
#'     "maaslin__Coussin_de_Lyon", 
#'     "metagenomeseq__Nazook", 
#'     "deseq__Pan_dulce"    
#'    )
#' )
#' expect_equal(
#'   steps_ids(rec, "da"), 
#'   c("maaslin__Coussin_de_Lyon", "metagenomeseq__Nazook", "deseq__Pan_dulce")
#' )
#' expect_equal(
#'   steps_ids(rec, "prepro"), 
#'   c("subset_taxa__Jalebi", "filter_taxa__Palmier")
#' )
#' expect_error(steps_ids(rec, "das"))
#' expect_type(steps_ids(rec), "character")
#' 
#' @examples 
#' data(test_rec)
#' 
#' ## We can extract the step identifiers from a recipe with `step_ids`
#' ids <- steps_ids(test_rec)
#' ids
#' 
#' ## With the `type` parameter, extract the prepro and da steps separately.
#' da_ids <- steps_ids(test_rec, type = "da")
#' da_ids
#' 
#' prepro_ids <- steps_ids(test_rec, type = "prepro")
#' prepro_ids
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
    "da" = purrr::discard(
      out, stringr::str_detect(out, "subset|filter|rarefaction")
    ),
    "prepro" = purrr::keep(
      out, stringr::str_detect(out, "subset|filter|rarefaction")
      )
  )
}

#' @keywords internal
#' @autoglobal
tick <- function() {
  crayon::green(cli::symbol$tick)
}

#' @keywords internal
#' @autoglobal
cross <- function() {
  crayon::red(cli::symbol$cross)
}

#' @keywords internal
#' @autoglobal
info <- function() {
  crayon::blue(cli::symbol$info)
}

#' @keywords internal
#' @autoglobal
dot <- function() {
  crayon::blue(cli::symbol$circle_filled)
}

#' Export step parameters as json.
#'
#' @param rec A recipe object.
#' @param file_name The path and file name of the optout file.
#'
#' @importFrom glue double_quote
#' @return invisible
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' 
#' ## Create a recipe with steps
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
#'   step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") %>%
#'   step_maaslin()
#'  
#' ## Prep recipe   
#' rec <- prep(rec, parallel = TRUE)
#' 
#' ## Export to json file
#' export_steps(rec, tempfile(fileext = ".json"))
export_steps <- function(rec, file_name) {
  
  inp <- rec@steps
  if (methods::is(rec, "prep_recipe")) {
    inp <- c(rec@steps, rec@bakes)
  } 
  
  to_cat <-
    inp %>%
    purrr::map_chr(~ {
      params <-
        names(.x) %>%
        purrr::map_chr(function(.y) {
          msg <- .x[[.y]]
          if (is.character(.x[[.y]]) | is.factor(.x[[.y]])) {
            msg <- double_quote(.x[[.y]])
          }
          stringr::str_c(
            "   ", double_quote(.y), ": ", paste0(msg, collapse = ""), ","
          )
        }) %>%
        stringr::str_c(collapse = "\n")

      stringr::str_c("{\n", params, "\n}")
    })

  cat(to_cat, file = file_name)
}

#' Import steps from json file
#'
#' @param rec A recipe object.
#' @param file Path to the input file.
#' @param parallel if FALSE, no palatalization. if TRUE, parallel execution
#'   using future and furrr packages.
#' @param workers Number of workers for palatalization.
#'
#' @return recipe-class object
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Initialize the recipe with a phyloseq object
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#' rec
#'
#' ## Import steps
#' json_file <- system.file("extdata", "test.json", package = "dar")
#' rec <- import_steps(rec, json_file)
#' rec
#' 
#' ## If the json file contains 'bake', the recipe is automatically prepared. 
#' json_file <- system.file("extdata", "test_bake.json", package = "dar")
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species") %>% 
#'   import_steps(json_file)
#'   
#' rec
#' cool(rec)
import_steps <- function(rec, file, parallel = TRUE, workers = 8) {
  lines <-
    readr::read_lines(file) %>%
    purrr::discard(stringr::str_detect(., "[{]|[}]"))

  id_idx <-  which(stringr::str_detect(lines, "id"))
  for (i in seq_along(id_idx)) {
    low_idx <- 1
    if (i != 1) { low_idx <- id_idx[i - 1] + 1 }
    if (stringr::str_detect(lines[id_idx[i]], "bake_")) { next }
    fun_name <- lines[id_idx[i]] %>%
      stringr::str_remove_all(".*: |,|\\\"") %>%
      stringr::str_remove_all("__.*")

    lines[low_idx:id_idx[i]] %>%
      purrr::map_chr(function(.x) {
        params <- .x %>%
          stringr::str_squish() %>%
          stringr::str_split(pattern = ": ") %>%
          unlist()

        param <- params[[1]] %>%
          stringr::str_remove_all("\\\"")

        value <- params[[2]] %>%
          stringr::str_squish() %>%
          stringr::str_remove_all(",$")

        if (stringr::str_detect(value, "%in%")) {
          value <- stringr::str_remove_all(value, "\\\"") %>%
            stringr::str_replace_all("\\\\", "\\\"")

          value <- encodeString(value, quote = "'")
        }

        if (stringr::str_count(value) == 0) { value <- expression(NULL) }
        stringr::str_c(param, " = ", paste0(value, sep = ""))

      }) %>% stringr::str_c(collapse = ", ") %>%
      stringr::str_c("rec <<- step_", fun_name, "(rec, ", ., ")") %>%
      parse(text = .) %>%
      eval()
  }

  if (any(stringr::str_detect(lines, "bake__"))) {
    rlang::inform(c(
      "!" = "bakes found in imported recipe",
      i = glue::glue("running {crayon::bgMagenta('prep()')}")
    ))

    rec <- prep(rec, parallel = parallel, workers = workers)
    for (i in seq_along(id_idx)) {
      low_idx <- 1
      if (i != 1) { low_idx <- id_idx[i - 1] + 1 }
      if (!stringr::str_detect(lines[id_idx[i]], "bake_")) { next }
      fun_name <- lines[id_idx[i]] %>%
        stringr::str_remove_all(".*: |,|\\\"") %>%
        stringr::str_remove_all("__.*")

      lines[low_idx:id_idx[i]] %>%
        purrr::map_chr(function(.x) {
          params <- .x %>%
            stringr::str_squish() %>%
            stringr::str_split(pattern = ": ") %>%
            unlist()

          param <- params[[1]] %>%
            stringr::str_remove_all("\\\"")

          value <- params[[2]] %>%
            stringr::str_squish() %>%
            stringr::str_remove_all(",$")

          if (stringr::str_detect(value, "%in%")) {
            value <- stringr::str_remove_all(value, "\\\"") %>%
              stringr::str_replace_all("\\\\", "\\\"")

            value <- encodeString(value, quote = "'")
          }

          if (stringr::str_count(value) == 0) { value <- expression(NULL) }
          stringr::str_c(param, " = ", paste0(value, sep = ""))
        }) %>% stringr::str_c(collapse = ", ") %>%
        stringr::str_c("rec <<- ", fun_name, "(rec, ", ., ")") %>%
        parse(text = .) %>%
        eval()
    }
  }
  rec
}

#' Checks if recipe contains a rarefaction step
#'
#' @param rec A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#'
#' @export
#' @autoglobal
#' @return boolean
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' rec <-
#'   phyloseq::subset_samples(
#'     GlobalPatterns, SampleType %in% c("Soil", "Skin")
#'   ) %>%
#'   recipe(var_info  = "SampleType", tax_info = "Genus") %>%
#'   step_rarefaction()
#'
#' contains_rarefaction(rec)
contains_rarefaction <- function(rec) {
  steps_ids(rec) %>%
    stringr::str_detect("^rarefaction_") %>%
    any()
}
