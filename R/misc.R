#' Make a random identification field for steps
#'
#' @param prefix A single character string
#' @return A character string with the prefix and random letters separated by
#'   and underscore.
#'
#' @export
#' @autoglobal
#' @tests
#' set.seed(123); first_id <- rand_id()
#' set.seed(123); expect_identical(rand_id(), first_id)
#' expect_match(first_id, "^step__"); expect_true(sub("^step__", "", first_id) %in% pastry_df$id)
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
#' @tests
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
#' dar:::get_comparisons("RiskGroup2", get_phy(test_rec))
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
#' @tests
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
#'   get_phy(test_rec) |>
#'   phyloseq::otu_table()
#'
#' dar:::to_tibble(otu_table)
to_tibble <- function(df, id_name = "otu_id") {
  df %>%
    data.frame(check.names = FALSE) %>%
    tibble::as_tibble(rownames = id_name)
}

#' Convert a step object into an unevaluated R call
#'
#' @description
#' This function takes a step from a `Recipe` and converts it into a native R
#' call (language object) using `rlang::call2`. This avoids the pitfalls of
#' pasting strings and parsing them, handling complex R objects like functions
#' or formulas naturally. The resulting call will have `rec` as its first
#' argument, followed by the parameters defined in the step.
#'
#' @param step A list or `step` class object containing the parameters for the
#'   operation, including an `id` string (e.g., `"method__randomstring"`).
#'
#' @return An unevaluated call (language object).
#'
#' @noRd
#' @keywords internal
#' @autoglobal
#' @tests
#' # 1. Test standard string and numeric parameters
#' step_standard <- step("maaslin",
#'   id = "maaslin__123",
#'   transform = "LOG",
#'   min_abundance = 0.1
#' )
#' call_standard <- step_to_call(step_standard)
#'
#' expect_true(is.call(call_standard))
#' expect_true(is.function(call_standard[[1]]))
#' expect_equal(call_standard[[1]], run_maaslin)
#' expect_equal(rlang::call_args(call_standard)$transform, "LOG")
#' expect_equal(rlang::call_args(call_standard)$min_abundance, 0.1)
#'
#' # 2. Test with a function parameter (The main reason for this refactor)
#' my_fun <- function(x) sum(x > 0) >= (0.03 * length(x))
#' step_func <- step("filter_taxa",
#'   id = "filter_taxa__abc",
#'   .f = my_fun
#' )
#' call_func <- step_to_call(step_func)
#'
#' expect_true(is.call(call_func))
#' expect_true(is.function(call_func[[1]]))
#' expect_equal(call_func[[1]], run_filter_taxa)
#' expect_true(is.function(rlang::call_args(call_func)$.f))
#'
#' # 3. Test with a formula parameter
#' step_formula <- step("deseq",
#'   id = "deseq__xyz",
#'   design = ~ RiskGroup2
#' )
#' call_formula <- step_to_call(step_formula)
#'
#' expect_true(is.call(call_formula))
#' expect_true(is.function(call_formula[[1]]))
#' expect_equal(call_formula[[1]], run_deseq)
#' expect_true(inherits(rlang::call_args(call_formula)$design, "formula"))
step_to_call <- function(step) {
  fn_obj <- step_runner(step)
  args <- step
  if (!"id" %in% names(formals(fn_obj))) {
    args[["id"]] <- NULL
  }
  rlang::call2(fn_obj, rec = quote(rec), !!!args)
}

#' Converts a step object into an executable R code string
#'
#' @description
#' This legacy function takes a step from a `Recipe` and converts its parameters
#' into a formatted text string that can be evaluated using `eval(parse(text = ...))`.
#' It handles the conversion of strings, NULLs, formulas, functions, and specific
#' cases like `weights` for `step_bake`.
#'
#' @param step A list or `step` class object containing the parameters for the
#'   operation, including an `id` string (e.g., `"method__randomstring"`).
#'
#' @return A character string containing the formatted R code (e.g.,
#'   `"rec %>% run_method(param = 'value')"`) ready to be parsed and evaluated.
#'
#' @noRd
#' @keywords internal
#' @autoglobal
step_to_expr <- function(step) {
  params <-
    step %>%
    purrr::discard(names(.) == "id" |
      (names(.) == "engine_args" & lengths(.) == 0L)) %>%
    purrr::map2_chr(names(.), ~ {

      # 1. NULL values
      if (is.null(.x)) {
        return(glue::glue("{.y} = NULL"))
      }

      # 2. Character values
      if (is.character(.x)) {
        .x <- stringr::str_c("'", .x, "'", collapse = ", ")
        return(glue::glue("{.y} = c({.x})"))
      }

      # 3. Formulas
      if (inherits(.x, "formula")) {
        return(paste0(.y, " = ", paste0(.x, collapse = "")))
      }

      # 4. Functions
      if (is.function(.x)) {
        func_str <- paste(deparse(.x), collapse = " ")
        return(glue::glue("{.y} = {func_str}"))
      }

      # 5. Named bake weights
      if (.y == "weights" && inherits(step, "step_bake") && !is.null(.x)) {
        text <- .x %>%
          purrr::map2_chr(names(.), ~ { paste0(.y, " = ", paste0(.x)) }) %>%
          stringr::str_c(collapse = ", ")
        return(glue::glue("{.y} = c({text})"))
      }

      # 6. Nested lists (for example engine_args)
      if (is.list(.x)) {
        return(glue::glue("{.y} = {step_value_expr(.x, .y)}"))
      }

      # 7. Numeric, logical and other scalar values
      glue::glue("{.y} = {.x}")
    }) %>%
    stringr::str_c(collapse = ", ")

  method <- step_method(step)

  # Build the final pipeline expression
  glue::glue("rec %>% run_{method}({params})")
}

#' Convert a persistable step value to reconstructible R code
#' @noRd
step_value_expr <- function(value, context = "step argument") {
  unsupported <- methods::is(value, "S4") || is.environment(value) ||
    typeof(value) %in% c("externalptr", "weakref") || inherits(value, "connection")
  if (unsupported) {
    cli::cli_abort(
      "Cannot export {.arg {context}} because it cannot be reconstructed faithfully.",
      class = "dar_error_unserializable_step"
    )
  }
  if (is.null(value)) {
    return("NULL")
  }
  if (is.function(value)) {
    captured <- codetools::findGlobals(value, merge = FALSE)$variables
    if (length(captured) > 0L) {
      cli::cli_abort(
        "Cannot export {.arg {context}} because the function captures external values: {.val {captured}}.",
        class = "dar_error_unserializable_step"
      )
    }
    return(paste(deparse(value), collapse = " ") |> stringr::str_squish())
  }
  if (inherits(value, "formula")) {
    return(paste(deparse(value), collapse = " ") |> stringr::str_squish())
  }
  if (is.list(value)) {
    if (length(value) == 0L) {
      return("list()")
    }
    value_names <- names(value)
    if (is.null(value_names) || anyNA(value_names) ||
        any(!nzchar(value_names)) || anyDuplicated(value_names)) {
      cli::cli_abort(
        "Cannot export {.arg {context}}: nested lists must have non-empty, unique names.",
        class = "dar_error_unserializable_step"
      )
    }
    entries <- purrr::map2_chr(value, value_names, function(item, name) {
      paste0(
        deparse(name), " = ",
        step_value_expr(item, paste0(context, "$", name))
      )
    })
    return(paste0("list(", paste(entries, collapse = ", "), ")"))
  }
  if (is.atomic(value)) {
    return(paste(utils::capture.output(dput(value)), collapse = " "))
  }
  cli::cli_abort(
    "Cannot export {.arg {context}} because its type is not supported.",
    class = "dar_error_unserializable_step"
  )
}

#' Finds common OTU between method results
#'
#' @param rec A Recipe object.
#' @param steps character vector with step ids to take in account
#'
#' @return tibble
#' @export
#' @autoglobal
#'
#' @tests
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
#' expect_equal(nrow(res_1), 52)
#' expect_equal(nrow(res_2), 27)
#' expect_equal(nrow(res_3), 41)
#' expect_s3_class(res_1, "tbl_df")
#'
#' @examples
#' data(test_prep_rec)
#'
#' ## From a PrepRecipe we can extract a tibble with all intersections
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
  df <- intersection_df(rec, steps) %>% tibble::as_tibble()
  if (!is.null(get_model(rec))) {
    return(
      df %>%
        tidyr::pivot_longer(
          dplyr::all_of(steps), names_to = "name", values_to = "value"
        ) %>%
        dplyr::filter(.data$value == 1) %>%
        dplyr::group_by(.data$taxa_id, .data$contrast_id, .data$effect) %>%
        dplyr::summarise(
          step_ids = paste(sort(unique(.data$name)), collapse = ", "),
          sum_methods = sum(.data$value), .groups = "drop"
        ) %>%
        dplyr::left_join(tax_table(rec), by = "taxa_id") %>%
        dplyr::arrange(dplyr::desc(.data$sum_methods))
    )
  }
  df %>%
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
#' @param rec A Recipe object.
#' @param type character vector indicating the type class. Options `c("all",
#'   "da", "prepro")`.
#' @param include_skipped For a prepared modeled recipe, include DA steps that
#'   were declared but skipped as statistically incompatible.
#'
#' @return character vector
#' @export
#' @autoglobal
#' @tests
#' data(test_prep_rec)
#' print(test_prep_rec) |> expect_snapshot()
#' rec <- test_prep_rec
#' expect_equal(
#'   steps_ids(rec),
#'   c("subset_taxa__Viennoiserie",
#'     "filter_taxa__Karakudamono",
#'     "maaslin__Welsh_cake",
#'     "deseq__Coussin_de_Lyon"
#'    )
#' )
#' expect_equal(
#'   steps_ids(rec, "da"),
#'   c("maaslin__Welsh_cake", "deseq__Coussin_de_Lyon")
#' )
#' expect_equal(
#'   steps_ids(rec, "prepro"),
#'   c("subset_taxa__Viennoiserie", "filter_taxa__Karakudamono")
#' )
#' expect_error(steps_ids(rec, "das"))
#' expect_type(steps_ids(rec), "character")
#'
#' @examples
#' data(test_rec)
#'
#' ## We can extract the step identifiers from a Recipe with `step_ids`
#' ids <- steps_ids(test_rec)
#' ids
#'
#' ## With the `type` parameter, extract the prepro and da steps separately.
#' da_ids <- steps_ids(test_rec, type = "da")
#' da_ids
#'
#' prepro_ids <- steps_ids(test_rec, type = "prepro")
#' prepro_ids
steps_ids <- function(rec, type = "all", include_skipped = FALSE) {
  check_any_recipe(rec)
  if (!type %in% c("all", "da", "prepro")) {
    cli::cli_abort(
      c(
        "x" = "Incorrect step type!",
        "i" = "Please use one of: {.val {c('all', 'da', 'prepro')}}."
      ),
      class = "dar_error_invalid_step_type"
    )

  }

  selected_steps <- switch(
    type,
    "all" = rec@steps,
    "da" = purrr::keep(rec@steps, is_da_step),
    "prepro" = purrr::keep(rec@steps, is_preprocessing_step)
  )
  out <- purrr::map_chr(selected_steps, ~ .x[["id"]])
  if (methods::is(rec, "PrepRecipe") && !include_skipped &&
      length(rec@execution) > 0L && type %in% c("all", "da")) {
    skipped <- rec@execution$skipped_steps$step_id %||% character()
    out <- setdiff(out, skipped)
  }
  out
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
#' @param rec A Recipe object.
#' @param file_name The path and file name of the optout file.
#'
#' @importFrom glue double_quote
#' @return invisible
#' @export
#' @autoglobal
export_steps <- function(rec, file_name) {
  check_any_recipe(rec)
  inp <- rec@steps
  if (methods::is(rec, "PrepRecipe")) {
    inp <- c(rec@steps, rec@bakes)
  }
  if (!is.null(get_model(rec))) {
    model <- get_model(rec)
    model$id <- "model__central"
    inp <- c(list(model), inp)
  }

  to_cat <-
    inp %>%
    purrr::map_chr(~ {
      params <-
        c(setdiff(names(.x), "id"), intersect("id", names(.x))) %>%
        purrr::map_chr(function(.y) {
          val <- .x[[.y]]

          if (identical(.y, "engine_args") && length(val) == 0L) {
            return(NA_character_)
          }

          # functions & formulas
          if (is.function(val) || inherits(val, "formula")) {
            val <-
              paste(deparse(val), collapse = " ") %>%
              stringr::str_squish() %>%
              glue::double_quote()

          # characters & factors
          } else if (is.character(val) || is.factor(val)) {
            val <- glue::double_quote(val)
            if (!is.null(names(.x[[.y]])) && all(nzchar(names(.x[[.y]])))) {
              val <- paste0(names(.x[[.y]]), " = ", val)
            }
            val <- val %>%
              stringr::str_c(collapse = ", ") %>%
              stringr::str_c("[c(", ., ")]")

          # NULL
          } else if (is.null(val)) {
            val <- "NULL"

          # recursively reconstructible named lists
          } else if (is.list(val)) {
            val <- paste0("[", step_value_expr(val, .y), "]")

          # reject values that the text format cannot preserve
          } else if (methods::is(val, "S4") || is.environment(val) ||
                     typeof(val) %in% c("externalptr", "weakref") ||
                     inherits(val, "connection")) {
            step_value_expr(val, .y)
          }

          stringr::str_c(
            "   ",
            glue::double_quote(.y), ": ", paste0(val, collapse = ""),
            ","
          )
        }) %>%
        stats::na.omit() %>%
        as.character() %>%
        stringr::str_c(collapse = "\n")

      stringr::str_c("{\n", params, "\n}")
    })

  writeLines(to_cat, file_name)
}

#' Import steps from json file
#'
#' @param rec A Recipe object.
#' @param file Path to the input file.
#' @param parallel if FALSE, no palatalization. if TRUE, parallel execution
#'   using future and furrr packages.
#' @param workers Number of workers for palatalization.
#'
#' @return recipe-class object
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' suppressWarnings(
#'   recipe(metaHIV_phy, "RiskGroup2", "Class") |>
#'     import_steps(
#'       system.file("extdata", "test_bake.json", package = "dar"),
#'       parallel = FALSE
#'     )
#' ) |>
#'  expect_snapshot()
#' @examples
#' data(metaHIV_phy)
#'
#' ## Initialize the Recipe with a phyloseq object
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' rec
#'
#' ## Import steps
#' json_file <- system.file("extdata", "test.json", package = "dar")
#' rec <- import_steps(rec, json_file)
#' rec
#'
#' ## If the json file contains 'bake', the Recipe is automatically prepared.
#' json_file <- system.file("extdata", "test_bake.json", package = "dar")
#' rec <-
#'   recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
#'   import_steps(json_file)
#'
#' rec
#' cool(rec)
import_steps <- function(
  rec,
  file,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16")
) {
  lines <-
    readr::read_lines(file) %>%
    purrr::discard(stringr::str_detect(., "[{]|[}]"))

  id_idx <- which(stringr::str_detect(lines, "id"))
  for (i in seq_along(id_idx)) {
    low_idx <- 1
    if (i != 1) {
      low_idx <- id_idx[i - 1] + 1
    }
    if (stringr::str_detect(lines[id_idx[i]], "bake_")) {
      next
    }
    fun_name <- lines[id_idx[i]] %>%
      stringr::str_remove_all(".*: |,|\\\"") %>%
      stringr::str_remove_all("__.*") %>%
      stringr::str_replace_all("\\[|\\]|c\\(", "")

    constructor <- if (identical(fun_name, "model")) {
      "add_model"
    } else {
      paste0("step_", fun_name)
    }
    instruction_lines <- lines[low_idx:id_idx[i]]
    if (identical(fun_name, "model")) {
      instruction_lines <- instruction_lines[-length(instruction_lines)]
    }

    extract_instructions(instruction_lines) %>%
      stringr::str_replace_all("\\[|\\]", "") %>%
      stringr::str_c(collapse = ", ") %>%
      stringr::str_c("rec <<- ", constructor, "(rec, ", ., ")") %>%
      parse(text = .) %>%
      eval()
  }

  if (any(stringr::str_detect(lines, "bake__"))) {
    cli::cli_inform(c(
      "!" = "bakes found in imported recipe",
      "i" = "running {.fun prep}"
    ))

    rec <- prep(rec, parallel = parallel, workers = workers)
    for (i in seq_along(id_idx)) {
      low_idx <- 1
      if (i != 1) {
        low_idx <- id_idx[i - 1] + 1
      }
      if (!stringr::str_detect(lines[id_idx[i]], "bake_")) {
        next
      }
      fun_name <- lines[id_idx[i]] %>%
        stringr::str_remove_all(".*: |,|\\\"") %>%
        stringr::str_remove_all("__.*") %>%
        stringr::str_replace_all("\\[|\\]|c\\(", "")

      extract_instructions(lines[low_idx:id_idx[i]]) %>%
        stringr::str_replace_all("\\[|\\]", "") %>%
        stringr::str_c(collapse = ", ") %>%
        stringr::str_c("rec <<- ", fun_name, "(rec, ", ., ")") %>%
        parse(text = .) %>%
        eval()
    }
  }
  rec
}

#' @noRd
#' @keywords internal
#' @autoglobal
extract_instructions <- function(lines) {
  lines %>%
    purrr::map_chr(function(.x) {
      params <- .x %>%
        stringr::str_squish() %>%
        stringr::str_split(pattern = ": ", n = 2) %>%
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
    })
}

#' Checks if Recipe contains a rarefaction step
#'
#' @param rec A Recipe object. The step will be added to the sequence of
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
#'   ) |>
#'   recipe() |>
#'   step_rarefaction()
#'
#' contains_rarefaction(rec)
contains_rarefaction <- function(rec) {
  check_any_recipe(rec)
  any(purrr::map_lgl(rec@steps, function(configured_step) {
    identical(step_method(configured_step), "rarefaction")
  }))
}


#' Perform Rarefaction on Phyloseq Object
#'
#' This function performs rarefaction on a phyloseq object if the `rarefy`
#' parameter is set to TRUE. Rarefaction is a process that randomly subsamples
#' the data to a specified depth. This is done to account for differences in
#' sequencing depth between samples. However, this process is not without
#' controversy. Rarefaction can lead to loss of information and can also lead to
#' false positives in differential abundance testing. For more information, see
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2
#'
#' @param phy A phyloseq object.
#' @param rarefy A logical value indicating whether to perform rarefaction. If
#'   'no_seed', rarefaction is performed without a set seed. If FALSE, no
#'   rarefaction is performed.
#'
#' @return A phyloseq object after rarefaction if `rarefy` is TRUE or "no_seed",
#'   otherwise the original phyloseq object is returned.
#' @keywords internal
#' @examples
#' data(metaHIV_phy)
#'
#' ## With seed
#' # phy_rarefied <- dar:::use_rarefy(metaHIV_phy, TRUE)
#'
#' ## Witout seed
#' # phy_rarefied <- dar:::use_rarefy(metaHIV_phy, "no_seed")
use_rarefy <- function(phy, rarefy) {
  if (isTRUE(rarefy)) {
    phy <- phyloseq::rarefy_even_depth(phy, rngseed = 1234, verbose = FALSE)
  }
  if (rarefy == "no_seed") {
    phy <- phyloseq::rarefy_even_depth(phy, verbose = FALSE)
  }
  phy
}

#' @noRd
#' @keywords internal
#' @autoglobal
rarefy_msg <- function(steps) {
  info_rarefy <- "Rarefaction is a process that randomly subsamples the data..."
  rarefy_vals <- purrr::map(steps, ~ .x[["rarefy"]])

  has_true <- any(purrr::map_lgl(rarefy_vals, ~ isTRUE(.x) || identical(.x, quote(T))))
  has_no_seed <- any(purrr::map_lgl(rarefy_vals, ~ identical(.x, "no_seed")))

  if (has_true) {
    cli::cli_inform(c(
      "!" = "Rarefaction is enabled (with set seed).",
      "i" = "This process is not without controversy. Use {.fun rarefaction_help} to get more info."
    ))
  }

  if (has_no_seed) {
    cli::cli_inform(c(
      "!" = "Rarefaction is enabled (without set seed).",
      "i" = "This process is not without controversy. Use {.fun rarefaction_help} to get more info."
    ))
  }

  steps
}


#' Information about the Rarefaction process
#'
#' Prints detailed information regarding the use of rarefaction in microbiome
#' analysis, its implications, and relevant literature.
#'
#' @return invisible
#' @export
#' @autoglobal
#' @examples
#' rarefaction_help()
rarefaction_help <- function() {
  cli::cli_rule(left = "{.strong Rarefaction in Microbiome Analysis}")
  cli::cli_text("")
  cli::cli_alert_info(c("{.strong What is Rarefaction?} It is a process that randomly subsamples the data to a specified depth to account for differences in sequencing depth between samples."))
  cli::cli_text("")
  cli::cli_alert_warning("{.strong The Controversy:}")
  cli::cli_bullets(c(
    "i" = "It can lead to an irreversible {.strong loss of information}.",
    "i" = "It can lead to {.strong false positives} in differential abundance testing."
  ))
  cli::cli_text("")
  cli::cli_alert_success("{.strong Reproducibility:} If performed with a {.emph set seed}, it ensures reproducible results, but removes randomness in the subsampling process.")
  cli::cli_alert_info("{.strong More info:} {.url https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0650-2}")
  cli::cli_text("")
  invisible(TRUE)
}
