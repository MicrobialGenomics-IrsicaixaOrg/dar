# GETTERS ----------------------------------------------------------------------

#' Returns var_info from Recipe-class object
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_var
#' @return Tibble containing `var_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' get_var(rec)
get_var <- function(rec) {
  check_any_recipe(rec)
  rec@var_info
}

#' Returns tax_info from Recipe-class object
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_tax
#' @return Tibble containing `tax_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' get_tax(rec)
get_tax <- function(rec) {
  check_any_recipe(rec)
  rec@tax_info
}

#' Returns phyloseq from Recipe-class object
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @aliases get_phy
#' @return Phyloseq class object
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' get_phy(rec)
get_phy <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq
}

# SETTERS ----------------------------------------------------------------------

#' Adds variable of interest to the Recipe
#'
#' @param rec A `Recipe` object (must not be prepped).
#' @param var_info A character string of column names corresponding to variables
#'   that will be used in any context.
#'
#' @aliases add_var
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' rec <- add_var(rec, var_info = "RiskGroup2")
#' rec
add_var <- function(rec, var_info) {
  check_recipe(rec) 
  rec@var_info <- tibble::tibble(vars = var_info)
  rec
}

#' Adds taxonomic level of interest in the Recipe.
#'
#' @param rec A `Recipe` object (must not be prepped).
#' @param tax_info A character string of taxonomic levels that will be used in
#'   any context.
#'
#' @aliases add_tax
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy)
#' rec <- add_tax(rec, tax_info = "Species")
#' rec
add_tax <- function(rec, tax_info) {
  check_recipe(rec) 
  rec@tax_info <- tibble::tibble(tax_lev = tax_info)
  rec
}

# PHYLOSEQ EXTRACTION ----------------------------------------------------------

#' Extracts tax_table from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' tax_table(rec)
tax_table <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq@tax_table %>%
    to_tibble("taxa_id") %>%
    dplyr::select(taxa_id, taxa = !!get_tax(rec)[[1]])
}

#' Extracts sample_data from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' sample_data(rec)
sample_data <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq %>%
    phyloseq::sample_data() %>%
    to_tibble("sample_id") %>%
    dplyr::select(sample_id, !!get_var(rec)[[1]])
}

#' Extracts otu_table from phyloseq inside a Recipe
#'
#' @param rec A `Recipe` or `PrepRecipe` object.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#' otu_table(rec)
otu_table <- function(rec) {
  check_any_recipe(rec)
  rec@phyloseq %>%
    phyloseq::otu_table() %>%
    to_tibble("taxa_id")
}


# EXECUTION (PREP) -------------------------------------------------------------

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_prep <- function(x, ...) { c("furrr", "future") }

#' Performs all the steps defined in a Recipe
#'
#' For a Recipe with at least one preprocessing or DA operation run the steps in
#' a convenient order.
#'
#' @param rec A `Recipe` object.
#' @param parallel if FALSE, no palatalization. if TRUE, parallel execution
#'   using future and furrr packages.
#' @param workers Number of workers for palatalization.
#' @param force Force the reexecution of all steps. This remove previous
#'   results.
#'
#' @aliases prep
#' @return A `PrepRecipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- 
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Class") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = function(x) sum(x > 0) >= (0.03 * length(x))) |>
#'   step_maaslin()
#'
#' da_results <- prep(rec)
prep <- function(rec,
                 parallel = TRUE,
                 workers = future::availableCores(constraints = "connections-16"),
                 force = FALSE) {
  
  if (inherits(rec, "PrepRecipe") && !force) {
    rlang::abort(c(
      "The input Recipe has already been prepped!",
      i = glue::glue("To force the rerun of all steps please run {crayon::bgMagenta('prep(rec, force = TRUE)')}")
    ))
  }
  
  check_any_recipe(rec)
  
  check <- utils::capture.output(required_deps(rec))
  if (length(check) > 0) {
    rlang::abort(c(
      "Not all necessary dependencies are installed.",
      i = glue::glue("Use {crayon::bgMagenta('required_deps(rec)')} to see how to install them.")
    ))
  }
  
  ## Phyloseq preprocessing steps
  calls_filter <- 
    rec@steps %>%
    purrr::keep(~ stringr::str_detect(.x[["id"]], "subset|filter")) %>% 
    purrr::map(step_to_call)
  
  calls_rarefaction <- 
    rec@steps %>%
    purrr::keep(~ stringr::str_detect(.x[["id"]], "run_rarefaction")) %>%
    purrr::map(step_to_call) 

  to_execute <- c(calls_filter, calls_rarefaction)
  
  for (.x in to_execute) { rec <- base::eval(.x) }
  
  ## DA steps
  da_steps <- 
    rec@steps %>% 
    purrr::discard(~ stringr::str_detect(.x[["id"]], "subset|filter|rarefaction")) %>%
    rarefy_msg()

  da_names <- purrr::map_chr(da_steps, ~ .x[["id"]])
  da_calls <- da_steps %>% purrr::map(step_to_call)
  
  if (parallel) {
    recipes_pkg_check(required_pkgs_prep(), "prep()")
    with(future::plan(future::multisession, workers = workers), local = TRUE)
    
    res <- 
      da_calls %>% 
      furrr::future_map(~ {
          rec <- rec
          base::eval(.x)
        },
        .options = furrr::furrr_options(
          seed = TRUE,
          packages = "dar",
          conditions = structure("condition", exclude = "message")
        )
      )
    
    names(res) <- da_names
    
  } else {
    res <- 
      da_calls %>%
      purrr::map(~ {
        rec <- rec
        base::eval(.x)
      })
      
    names(res) <- da_names
  }
  
  prep_recipe(rec, res, list())
}


# RESULTS EXTRACTION -----------------------------------------------------------

#' Returns data.frame with OTU intersection between methods
#'
#' @param rec A `PrepRecipe` object.
#' @param steps character vector with step_ids to take in account.
#' @param tidy Boolean indicating if result must be in tidy format.
#'
#' @aliases intersection_df
#' @return data.frame class object
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' df <- intersection_df(test_prep_rec)
#' head(df)
intersection_df <- function(rec, steps = steps_ids(rec, "da"), tidy = FALSE) {
  
  check_prep_recipe(rec)
  
  df <- 
    names(rec@results) %>%
    purrr::keep(. %in% steps) %>%
    purrr::set_names() %>%
    purrr::map_dfc( ~ {
      taxa <- rec@results[[.x]][[1]] %>% 
        dplyr::filter(signif == TRUE) %>% 
        dplyr::pull(taxa_id)
        
      rownames(rec@phyloseq@otu_table) %>%
        tibble::tibble(taxa_id = .) %>%
        dplyr::mutate(!!.x := dplyr::if_else(taxa_id %in% taxa, 1, 0)) %>%
        dplyr::select(!!.x)
    }) %>%
    dplyr::mutate(taxa_id = rownames(rec@phyloseq@otu_table), .before = 1) %>%
    as.data.frame()
  
  if (tidy) { df <- tidyr::pivot_longer(df, -taxa_id) }
  df
}

#' Overlap of significant OTUs between tested methods.
#'
#' @param rec A `PrepRecipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param type Indicates whether to use all taxa ("all") or only those that are
#'   differentially abundant in at least one method ("da"). Default as "all". 
#'
#' @aliases overlap_df
#' @return df
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' df <- overlap_df(test_prep_rec, steps_ids(test_prep_rec, "da"))
#' head(df)
overlap_df <- function(rec, steps = steps_ids(rec, "da"), type = "all") {
  
  check_prep_recipe(rec)
  
  df <- 
    intersection_df(rec) %>%  
    tibble::as_tibble() %>% 
    dplyr::select(dplyr::all_of(steps))
  
  if (type == "da") {
    df <- 
      dplyr::rowwise(df) %>% 
      dplyr::mutate(sum = sum(dplyr::across(dplyr::all_of(steps))), .before = 1) %>% 
      dplyr::filter(sum != 0) %>% 
      dplyr::select(dplyr::all_of(steps))
  }
    
  names(df) %>%
    purrr::map_dfr(~ {
      names(df) %>%
        purrr::map_dfc(function(.y) {
          res <- tibble::tibble(
            var_1 = dplyr::pull(df, .x),
            var_2 = dplyr::pull(df, .y),
            sum = var_1 - var_2
          )
          
          tibble::tibble(
            !!dplyr::sym(.y) := nrow(dplyr::filter(res, sum == 0)) / nrow(res)
          )
        })
    }) %>%
    data.frame(row.names = names(df))
}

#' Extract results from defined bake
#'
#' @param rec A `PrepRecipe` object.
#' @param bake Name or index of the bake to extract.
#'
#' @aliases cool
#' @return tbl_df
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#' rec <- bake(test_prep_rec)
#' cool(rec)
cool <- function(rec, bake = 1) {
  
  check_prep_recipe(rec)
  
  all_bakes <- rec@bakes
  all_names <- all_bakes %>% purrr::map_chr(~ as.character(.x[["id"]]))
  
  if (is.numeric(bake) && length(all_bakes) < bake) {
    rlang::abort(c(
      "Bake index is not defined in the PrepRecipe!",
      glue::glue("Run {crayon::bgMagenta('bake(PrepRecipe)')} and then try with {crayon::bgMagenta('cool()')}")
    ))
  }
  
  if (!is.numeric(bake) && !bake %in% all_names) {
    rlang::abort(c(
      "Bake name is not defined in the PrepRecipe!",
      glue::glue("Run {crayon::bgMagenta('bake(PrepRecipe)')} and then try with {crayon::bgMagenta('cool()')}")
    ))
  }
  
  if (!is.numeric(bake)) {
    bake <- which(all_names == bake)
  }
  
  to_execute <- 
    rec@bakes %>%
    .[bake] %>%
    purrr::map_chr(step_to_expr)
  
  eval(parse(text = to_execute))
}
