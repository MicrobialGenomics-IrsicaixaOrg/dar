# CLASS RECIPE -----------------------------------------------------------------

## class unions ----

#' @rdname Recipe-class
methods::setClassUnion("phyloseq_or_null", c("phyloseq", "NULL"))

#' @rdname Recipe-class
methods::setClassUnion("list_or_NULL", c("list", "NULL"))

## class def ----

#' Recipe-class object
#'
#' A Recipe is a description of the steps to be applied to a data set in order
#' to prepare it for data analysis.
#'
#' @slot phyloseq Phyloseq-class object.
#' @slot steps List of step-class objects that will be used by DA.
#' @slot model Optional centralized statistical model specification created by
#'   [add_model()].
#'
#' @name Recipe-class
#' @rdname Recipe-class
#' @exportClass Recipe
#' @return Recipe-class object
#' @autoglobal
methods::setClass(
  Class = "Recipe",
  slots = c(
    phyloseq = "phyloseq_or_null",
    model = "list_or_NULL",
    steps = "list"
  ),
  prototype = list(
    phyloseq = NULL,
    model = NULL,
    steps = list()
  )
)

## constructor ----

#' Create a Recipe for preprocessing data
#'
#' A Recipe is a description of the steps to be applied to a data set in order
#' to prepare it for data analysis.
#'
#' @param microbiome_object Phyloseq-class object or
#'   TreeSummarizedExperiment-class object.
#' @param steps list with steps.
#'
#' @return An object of class `Recipe` containing the microbiome object, an
#'   optional centralized model and the configured processing or DA steps.
#'
#' @aliases Recipe
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' data(GlobalPatterns, package = "mia")
#'
#' # 1. Error: Invalid microbiome_object type
#' expect_error(
#'   recipe(data.frame(a = 1:5)),
#'   class = "dar_error_invalid_microbiome_object"
#' )
#'
#' # 2. Error: Invalid taxonomy rank names (fails with made-up names)
#' bad_phy <- metaHIV_phy
#' tax_tab <- phyloseq::tax_table(bad_phy)
#' colnames(tax_tab) <- paste0("BadRank", seq_len(ncol(tax_tab)))
#' phyloseq::tax_table(bad_phy) <- tax_tab
#'
#' expect_error(
#'   recipe(bad_phy),
#'   class = "dar_error_invalid_rank_names"
#' )
#'
#' # 3. Success: Valid taxonomy ranks in UPPERCASE (tests stringr normalization)
#' upper_phy <- metaHIV_phy
#' tax_tab_up <- phyloseq::tax_table(upper_phy)
#' colnames(tax_tab_up) <- toupper(colnames(tax_tab_up))
#' phyloseq::tax_table(upper_phy) <- tax_tab_up
#'
#' expect_s4_class(recipe(upper_phy), "Recipe")
#'
#' # 4. Success: Valid TreeSummarizedExperiment
#' expect_s4_class(recipe(GlobalPatterns), "Recipe")
#'
#' # 5. Success: Valid phyloseq
#' expect_s4_class(recipe(metaHIV_phy), "Recipe")
#'
#' @examples
#' data(metaHIV_phy)
#'
#' ## Configure, prepare and obtain a consensus
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Phylum") |>
#'   step_filter_by_prevalence(threshold = 0.3) |>
#'   step_deseq() |>
#'   step_maaslin()
#'
#' da_results <- prep(rec, parallel = FALSE) |>
#'   bake()
#' cool(da_results)
#'
#' ## A recipe without a model can be used for preprocessing
#' rec <- recipe(metaHIV_phy) |>
#'   step_filter_by_prevalence(threshold = 0.1)
recipe <- function(microbiome_object = NULL, steps = list()) {
  if (!is(microbiome_object, "phyloseq") &&
      !is(microbiome_object, "TreeSummarizedExperiment")) {
    cli::cli_abort(
      c(
        "x" = "{.arg microbiome_object} must be a {.cls phyloseq} or {.cls TreeSummarizedExperiment} object.",
        "i" = "You supplied an object of class: {.cls {class(microbiome_object)}}."
      ),
      class = "dar_error_invalid_microbiome_object"
    )
  }

  if (is(microbiome_object, "TreeSummarizedExperiment")) {
    microbiome_object <- mia::convertToPhyloseq(microbiome_object)
  }

  if (!is.null(microbiome_object@tax_table)) {
    tax_names <- phyloseq::rank_names(microbiome_object) %>% stringr::str_to_sentence()
    tax_tbl <- phyloseq::tax_table(microbiome_object)
    colnames(tax_tbl) <- tax_names
    phyloseq::tax_table(microbiome_object) <- tax_tbl
    tax_expected <- mia::getTaxonomyRanks() %>% stringr::str_to_sentence()
    if (!all(tax_names %in% tax_expected)) {
      cli::cli_abort(
        c(
          "x" = "{.arg rank} must be a value from {.fun taxonomyRanks}.",
          "i" = "Rename the columns from the {.code tax_table} slot of your input {.cls phyloseq} or {.cls TreeSummarizedExperiment} with standard names.",
          "i" = "Standard names: {.val {mia::getTaxonomyRanks()}}."
        ),
        class = "dar_error_invalid_rank_names"
      )
    }
  }

  methods::new(
    Class = "Recipe",
    phyloseq = microbiome_object,
    model = NULL,
    steps = steps
  )
}

## validity ----

#' Validate the structural invariants of a Recipe
#'
#' @param object A `Recipe` object.
#'
#' @return `TRUE` when the object is valid, otherwise a character vector with
#'   every detected problem.
#' @noRd
#' @keywords internal
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#'
#' complete_rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#' expect_true(methods::validObject(complete_rec, test = TRUE))
#' expect_true(methods::validObject(recipe(metaHIV_phy), test = TRUE))
#'
#' testthat::local_mocked_bindings(
#'   nsamples = function(x) 0,
#'   ntaxa = function(x) 0,
#'   .package = "phyloseq"
#' )
#' empty_problems <- recipe_validity_problems(complete_rec)
#' expect_true(any(grepl("at least one sample", empty_problems)))
#' expect_true(any(grepl("at least one taxon", empty_problems)))
recipe_validity_problems <- function(object) {
  problems <- character()
  phy <- object@phyloseq

  if (is.null(phy)) {
    problems <- c(problems, "`phyloseq` must not be NULL")
    return(problems)
  }

  if (phyloseq::nsamples(phy) < 1) {
    problems <- c(problems, "`phyloseq` must contain at least one sample")
  }
  if (phyloseq::ntaxa(phy) < 1) {
    problems <- c(problems, "`phyloseq` must contain at least one taxon")
  }

  problems <- c(
    problems,
    step_slot_validation_problems(object@steps, "steps", c("preprocessing", "da")),
    duplicate_step_id_problems(object@steps)
  )

  if (!is.null(object@model)) {
    problems <- c(problems, model_validity_problems(object))
  }

  if (length(problems) == 0) TRUE else problems
}

methods::setValidity(Class = "Recipe", method = recipe_validity_problems)

## printing ----

methods::setMethod("show", signature = "Recipe", definition = function(object) {
  cli::cat_rule(crayon::blue("DAR Recipe"))
  cat("Inputs:\n\n")

  ## Phyloseq
  phy <- get_phy(object)
  ntax <- phyloseq::ntaxa(phy)
  nsam <- phyloseq::nsamples(phy)
  cat(
    glue::glue(
      "     {info()} phyloseq object with {crayon::blue(ntax)} taxa and ",
      "{crayon::blue(nsam)} samples"
    ),
    "\n"
  )

  ## Variable

  targets <- recipe_targets(object)
  if (length(targets) == 0) {
    cat(
      glue::glue(
        "     {cross()} undefined analysis target. Use ",
        "{crayon::bgMagenta('add_model()')} to define the analysis!"
      ),
      "\n"
    )
  } else {
    var <- targets[[1]]
    var_vals <- analysis_sample_data(object) %>% dplyr::pull(.env$var)
    if (is.character(var_vals) | is.factor(var_vals)) {
      levs <- factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
      msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
    }
    if (is.numeric(var_vals)) { msg <- glue::glue("class: numeric") }
    cat(glue::glue("     {info()} variable of interes {crayon::blue(var)} ({msg})"), "\n")
  }

  ## Taxa
  tax_level <- recipe_tax_level(object)
  if (length(tax_level) == 0) {
    cat(
      glue::glue(
        "     {cross()} undefined taxonomic level. Use ",
        "{crayon::bgMagenta('add_model()')} to define the analysis!"
      ),
      "\n"
    )
  } else {
    cat(glue::glue("     {info()} taxonomic level {crayon::blue(tax_level)}"), "\n\n")
  }

  if (!is.null(get_model(object))) {
    model_text <- paste(deparse(get_model(object)$formula), collapse = " ")
    cat("Statistical model:\n\n")
    cat(glue::glue("     {info()} {crayon::blue(model_text)}"), "\n\n")
  }

  ## Steps
  if (length(object@steps) > 0) {
    cat("Preprocessing steps:\n\n")
    object@steps %>%
      purrr::walk(
        ~ {
          if (is_preprocessing_step(.x)) {
            id <-
              glue::glue("id = {.x[['id']]}") %>%
              crayon::silver()

            class(.x)[[1]]
            cat(c(glue::glue("     {dot()} {class(.x)[[1]]}() {id}"), "\n"))
          }
        }
      )
  }
  cat("\n")
  if (length(object@steps) > 0) {
    cat("DA steps:\n\n")
    object@steps %>%
      purrr::walk(
        ~ {
          if (is_da_step(.x)) {
            id <-
              glue::glue("id = {.x[['id']]}") %>%
              crayon::silver()

            class(.x)[[1]]
            cat(c(glue::glue("     {dot()} {class(.x)[[1]]}() {id}"), "\n"))
          }
        }
      )
  }
})


# CLASS PrepRecipe ------------------------------------------------------------

## class def ----

#' PrepRecipe-class object
#'
#' A PrepRecipe is Recipe with the results corresponding to the steps defined
#' in the Recipe.
#'
#' @slot results Contains the results of all defined analysis in the Recipe.
#' @slot bakes Contains the executed bakes.
#' @slot execution Manifest containing resolved contrasts, samples removed by
#'   the missing-value policy, and executed or skipped methods.
#'
#' @name PrepRecipe-class
#' @rdname PrepRecipe-class
#' @exportClass PrepRecipe
#' @return PrepRecipe-class object
#' @autoglobal
methods::setClass(
  Class = "PrepRecipe",
  contains = "Recipe",
  slots = c(results = "list", bakes = "list", execution = "list"),
  prototype = list(execution = list())
)

#' @noRd
prep_recipe_validity_problems <- function(object) {
  problems <- c(
    step_slot_validation_problems(object@bakes, "bakes", "bake"),
    duplicate_step_id_problems(c(object@steps, object@bakes), "`steps` and `bakes`")
  )
  if (length(problems) == 0L) TRUE else problems
}

methods::setValidity(Class = "PrepRecipe", method = prep_recipe_validity_problems)

## constructor ----

#' Create a PrepRecipe.
#'
#' A PrepRecipe is Recipe with the results corresponding to the steps defined
#' in the Recipe.
#'
#' @param rec A Recipe object.
#' @param results list with the results
#' @param bakes list with saved bakes
#'
#' @return An object of class `PrepRecipe`.
#' @autoglobal
#' @keywords internal
#' @aliases PrepRecipe
#' @tests
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species")
#'
#' prepped <- prep_recipe(rec, results = list(), bakes = list())
#' expect_s4_class(prepped, "PrepRecipe")
#' expect_true(methods::validObject(prepped, test = TRUE))
#'
#' rec@model$tax_level <- "Missing_rank"
#' expect_error(
#'   prep_recipe(rec, results = list(), bakes = list()),
#'   regexp = "Missing_rank"
#' )
prep_recipe <- function(rec, results, bakes, execution = list()) {
  methods::new(
    Class = "PrepRecipe",
    results = results,
    bakes = bakes,
    execution = execution,
    rec
  )
}

## printing ----

#' @param object A Recipe object.
#'
#' @rdname Recipe-class
#' @autoglobal
methods::setMethod(
  "show",
  signature = "PrepRecipe",
  definition = function(object) {
    cli::cat_rule(crayon::blue("DAR Results"))
    cat("Inputs:\n\n")

    ## Phyloseq
    phy <- get_phy(object)
    ntax <- phyloseq::ntaxa(phy)
    nsam <- phyloseq::nsamples(phy)
    cat(
      glue::glue(
        "     {info()} phyloseq object with {crayon::blue(ntax)} taxa and ",
        "{crayon::blue(nsam)} samples"
      ),
      "\n"
    )

    ## Variable
    targets <- recipe_targets(object)
    var <- if (length(targets) > 0L) targets else NULL

    if (is.null(var)) {
      cat(
        glue::glue(
          "     {cross()} undefined analysis target. Use ",
          "{crayon::bgMagenta('add_model()')} to define the analysis!"
        ),
        "\n"
      )
    } else {
      var <- targets[[1]]
      var_vals <- analysis_sample_data(object) %>% dplyr::pull(.env$var)
      if (is.character(var_vals) | is.factor(var_vals)) {
        levs <-
          factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
        msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
      }
      if (is.numeric(var_vals)) { msg <- glue::glue("class: numeric") }
      cat(
        glue::glue("     {info()} variable of interes {crayon::blue(var)} ({msg})"),
        "\n"
      )
    }

    ## Taxa
    tax_level <- recipe_tax_level(object)
    var <- if (length(tax_level) > 0L) tax_level else NULL

    if (is.null(var)) {
      cat(
        glue::glue(
          "     {cross()} undefined taxonomic level. Use ",
          "{crayon::bgMagenta('add_model()')} to define the analysis!"
        ),
        "\n"
      )
    } else {
      cat(glue::glue("     {info()} taxonomic level {crayon::blue(tax_level)}"), "\n\n")
    }

    ## Results
    cat("Results:\n\n")
    names(object@results) %>%
      purrr::walk(
        ~ {
          n_taxa <-
            tidy_results(object, steps = .x, significant_only = TRUE) %>%
            dplyr::pull(taxa_id) %>%
            unique() %>%
            length()

          n_taxa <- crayon::silver(glue::glue("diff_taxa = {n_taxa}"))

          cat(c(glue::glue("     {tick()} {.x} {n_taxa}"), "\n"))
        }
      )

    if (length(object@execution$skipped_steps) > 0L &&
        nrow(object@execution$skipped_steps) > 0L) {
      cat("\nSkipped methods:\n\n")
      purrr::pwalk(
        object@execution$skipped_steps,
        function(step_id, engine, compatible, reason, ...) {
          cat(glue::glue("     {cross()} {step_id}: {reason}"), "\n")
        }
      )
    }

    if (length(object@results) > 0) {
      n_overlap <-
        find_intersections(object) %>%
        dplyr::filter(sum_methods == length(steps_ids(object, type = "da"))) %>%
        nrow()

      cli::cat_line()
      cat(
        glue::glue(
          "     {info()} {n_overlap} taxon-contrast effects are present in all tested methods"
        ),
        "\n\n"
      )
    }

    ## Bakes
    if (length(object@bakes) > 0) {
      cat("Bakes:\n\n")

      object@bakes %>%
        purrr::iwalk(
          ~ {
            msg <-
              .x %>%
              purrr::map2_chr(
                names(.),
                ~ {
                  if (is.null(.x)) { .x <- "NULL" }
                  if (.y == "weights" & .x != "NULL") {
                    .x <-
                      glue::glue("{names(.x)} = {.x}") %>%
                      stringr::str_c(collapse = ", ") %>%
                      stringr::str_c("c(", ., ")")
                  }
                  glue::glue("{.y}: {.x}")
                }
              ) %>% stringr::str_c(collapse = ", ")

            cat(
              c(
                glue::glue(
                  "     {dot()} {crayon::blue(crayon::bold(paste0(.y, ' ->')))} ",
                  "{crayon::silver(msg)}"
                ),
                "\n"
              )
            )
          }
        )
    }
  }
)
