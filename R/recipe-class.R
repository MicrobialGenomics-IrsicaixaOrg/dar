# CLASS RECIPE -----------------------------------------------------------------

## class unions ----

#' @rdname Recipe-class
methods::setClassUnion("phyloseq_or_null", c("phyloseq", "NULL"))

#' @rdname Recipe-class
methods::setClassUnion("tibble_or_NULL", c("tbl_df", "NULL"))

#' @rdname Recipe-class
methods::setClassUnion("list_or_NULL", c("list", "NULL"))

## class def ----

#' Recipe-class object
#'
#' A Recipe is a description of the steps to be applied to a data set in order
#' to prepare it for data analysis.
#'
#' @slot phyloseq Phyloseq-class object.
#' @slot var_info A tibble that contains the current set of terms in the data
#'   set. This initially defaults to the same data contained in `var_info`.
#' @slot tax_info A tibble that contains the current set of taxonomic levels
#'   that will be used in the analysis.
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
    var_info = "tibble_or_NULL",
    tax_info = "tibble_or_NULL",
    model = "list_or_NULL",
    steps = "list"
  ),
  prototype = list(
    phyloseq = NULL,
    var_info = NULL,
    tax_info = NULL,
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
#' @param var_info A character string of column names corresponding to variables
#'   that will be used in any context.
#' @param tax_info A character string of taxonomic levels that will be used in
#'   any context.
#' @param steps list with steps.
#'
#' @return An object of class `Recipe` with sub-objects: \item{phyloseq}{object
#'   of class `phyloseq` with taxa abundance information.} \item{var_info}{A
#'   tibble that contains the current set of terms in the data set. This
#'   initially defaults to the same data contained in `var_info`.}
#'   \item{tax_info}{A tibble that contains the current set of taxonomic levels
#'   that will be used in the analysis.}
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
#' # 4. Error: Invalid var_info missing in metadata
#' expect_error( 
#'   recipe(metaHIV_phy, var_info = "error_var", tax_info = "Species"),
#'   class = "dar_error_missing_vars"
#' )
#' 
#' # 5. Error: Invalid tax_info missing in tax_table
#' expect_error( 
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "error_tax"),
#'   class = "dar_error_missing_tax"
#' )
#' 
#' # 6. Success: Valid TreeSummarizedExperiment
#' expect_s4_class(recipe(GlobalPatterns), "Recipe")
#' 
#' # 7. Success: Valid phyloseq
#' expect_s4_class(recipe(metaHIV_phy), "Recipe")
#' 
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") |>
#'   step_deseq() |>
#'   step_maaslin()
#'
#' ## Prep recipe
#' da_results <- prep(rec)
#'
#' ## Consensus strategy
#' n_methods <- 2
#' da_results <- bake(da_results, count_cutoff = n_methods)
#'
#' ## Results
#' cool(da_results)
#'
#' ## You can also crate a recipe without var and tax info
#' rec <- recipe(metaHIV_phy)
#'
#' rec
#'
#' ## And define them later
#' rec <- rec |>
#'   add_var("RiskGroup2") |>
#'   add_tax("Genus")
#'
#' rec
#'
#' ## When trying to add an identical step to an existing one, the system
#' ## returns an information message.
#' rec <- step_aldex(rec)
#' rec <- step_aldex(rec)
#'
#' ## The same with bake
#' da_results <- bake(da_results)
#' da_results <- bake(da_results)
recipe <- function(microbiome_object = NULL, 
                   var_info = NULL, 
                   tax_info = NULL, 
                   steps = list()) {
  
  if (!is(microbiome_object, "phyloseq") && 
      !is(microbiome_object, "TreeSummarizedExperiment")) {
    cli::cli_abort(c(
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
      cli::cli_abort(c(
        "x" = "{.arg rank} must be a value from {.fun taxonomyRanks}.",
        "i" = "Rename the columns from the {.code tax_table} slot of your input {.cls phyloseq} or {.cls TreeSummarizedExperiment} with standard names.",
        "i" = "Standard names: {.val {mia::getTaxonomyRanks()}}."
      ),
      class = "dar_error_invalid_rank_names"
      )
    }
  }
  
  if (!is.null(var_info)) {
    s_data <- as(phyloseq::sample_data(microbiome_object), "data.frame")
    if (!all(var_info %in% colnames(s_data))) {
      missing_vars <- setdiff(var_info, colnames(s_data))
      cli::cli_abort(c(
        "x" = "Column{?s} {.val {missing_vars}} {?does not/do not} exist in the {.arg microbiome_object} metadata.",
        "i" = "Available columns: {.val {colnames(s_data)}}."
      ),
      class = "dar_error_missing_vars"
      )
    }
  }

  if (!is.null(tax_info)) {
    tax_info <- stringr::str_to_sentence(tax_info)
    tax_ranks <- phyloseq::rank_names(microbiome_object)
    if (!all(tax_info %in% tax_ranks)) {
      missing_tax <- setdiff(tax_info, tax_ranks)
      cli::cli_abort(c(
        "x" = "Taxonomic rank{?s} {.val {missing_tax}} {?does not/do not} exist in the {.arg microbiome_object}.",
        "i" = "Available ranks: {.val {tax_ranks}}."
      ),
      class = "dar_error_missing_tax"
      )
    }
  }

  var_info <- tibble::tibble(vars = var_info)
  tax_info <- tibble::tibble(tax_lev = tax_info)

  methods::new(
    Class = "Recipe",
    phyloseq = microbiome_object,
    var_info = var_info,
    tax_info = tax_info,
    model = NULL,
    steps = steps
  )
}

## validity ----

#' Validate the structural invariants of a Recipe
#'
#' `var_info` and `tax_info` may be unset because users can add them after
#' creating a recipe. When present, however, their values must describe the
#' data stored in the recipe.
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
#' complete_rec <- recipe(
#'   metaHIV_phy,
#'   var_info = "RiskGroup2",
#'   tax_info = "Species"
#' )
#' incomplete_rec <- recipe(metaHIV_phy)
#'
#' expect_true(methods::validObject(complete_rec, test = TRUE))
#' expect_true(methods::validObject(incomplete_rec, test = TRUE))
#'
#' invalid_var <- complete_rec
#' invalid_var@var_info <- tibble::tibble(vars = "missing_variable")
#' expect_match(
#'   methods::validObject(invalid_var, test = TRUE),
#'   "missing_variable"
#' )
#'
#' invalid_tax <- complete_rec
#' invalid_tax@tax_info <- tibble::tibble(tax_lev = "Missing_rank")
#' expect_match(
#'   methods::validObject(invalid_tax, test = TRUE),
#'   "Missing_rank"
#' )
#'
#' malformed_var <- complete_rec
#' malformed_var@var_info <- tibble::tibble(wrong = "RiskGroup2")
#' expect_match(
#'   methods::validObject(malformed_var, test = TRUE),
#'   "var_info"
#' )
#'
#' malformed_tax <- complete_rec
#' malformed_tax@tax_info <- tibble::tibble(wrong = "Species")
#' expect_match(
#'   methods::validObject(malformed_tax, test = TRUE),
#'   "tax_info"
#' )
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

  var_info <- object@var_info
  if (!is.null(var_info) && ncol(var_info) > 0) {
    if (!identical(names(var_info), "vars")) {
      problems <- c(
        problems,
        "`var_info` must contain only a `vars` column"
      )
    } else if (!is.character(var_info$vars)) {
      problems <- c(problems, "`var_info$vars` must be a character vector")
    } else {
      sample_vars <- tryCatch(
        colnames(as(phyloseq::sample_data(phy), "data.frame")),
        error = function(cnd) character()
      )
      missing_vars <- setdiff(var_info$vars, sample_vars)
      if (length(missing_vars) > 0) {
        problems <- c(
          problems,
          glue::glue(
            "Variables not found in `sample_data`: {paste(missing_vars, collapse = ', ')}"
          )
        )
      }
    }
  }

  tax_info <- object@tax_info
  if (!is.null(tax_info) && ncol(tax_info) > 0) {
    if (!identical(names(tax_info), "tax_lev")) {
      problems <- c(
        problems,
        "`tax_info` must contain only a `tax_lev` column"
      )
    } else if (!is.character(tax_info$tax_lev)) {
      problems <- c(problems, "`tax_info$tax_lev` must be a character vector")
    } else {
      missing_tax <- setdiff(tax_info$tax_lev, phyloseq::rank_names(phy))
      if (length(missing_tax) > 0) {
        problems <- c(
          problems,
          glue::glue(
            "Taxonomic ranks not found in `tax_table`: {paste(missing_tax, collapse = ', ')}"
          )
        )
      }
    }
  }

  if (!is.null(object@model)) {
    problems <- c(problems, model_validity_problems(object))
  }

  if (length(problems) == 0) TRUE else problems
}

methods::setValidity(
  Class = "Recipe",
  method = recipe_validity_problems
)

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

  if (length(get_var(object)) == 0) {
    cat(
      glue::glue(
        "     {cross()} undefined variable of interest. Use ",
        "{crayon::bgMagenta('add_var()')} to add it to Recipe!"
      ),
      "\n"
    )
  } else {
    var <- get_var(object) %>% dplyr::pull(1)
    var_vals <- sample_data(object) %>% dplyr::pull(.env$var)
    if (is.character(var_vals) | is.factor(var_vals)) {
      levs <- factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
      msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
    }
    if (is.numeric(var_vals)) {
      msg <- glue::glue("class: numeric")
    }
    cat(glue::glue(
      "     {info()} variable of interes {crayon::blue(var)} ({msg})"
      ),
      "\n"
    )
  }

  ## Taxa
  if (length(get_tax(object)) == 0) {
    cat(
      glue::glue(
        "     {cross()} undefined taxonomic level. Use ",
        "{crayon::bgMagenta('add_tax()')} to add it to Recipe!"
      ),
      "\n"
    )
  } else {
    cat(glue::glue(
      "     {info()} taxonomic level {crayon::blue(get_tax(object))}"),
      "\n\n"
    )
  }

  if (!is.null(get_model(object))) {
    model_text <- paste(deparse(get_model(object)$formula), collapse = " ")
    cat("Statistical model:\n\n")
    cat(glue::glue("     {info()} {crayon::blue(model_text)}"), "\n\n")
  }

  ## Steps
  if (length(object@steps) > 0) {
    cat("Preporcessing steps:\n\n")
    object@steps %>%
      purrr::walk(~ {
        if (stringr::str_detect(.x[['id']], "subset|filter|rarefaction")) {
          id <-
            glue::glue("id = {.x[['id']]}") %>%
            crayon::silver()

          class(.x)[[1]]
          cat(c(glue::glue("     {dot()} {class(.x)[[1]]}() {id}"), "\n"))
        }
      })
  }
  cat("\n")
  if (length(object@steps) > 0) {
    cat("DA steps:\n\n")
    object@steps %>%
      purrr::walk(~ {
        if (!stringr::str_detect(.x[['id']], "subset|filter|rarefaction")) {
          id <-
            glue::glue("id = {.x[['id']]}") %>%
            crayon::silver()

          class(.x)[[1]]
          cat(c(glue::glue("     {dot()} {class(.x)[[1]]}() {id}"), "\n"))
        }
      })
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
#' rec <- recipe(metaHIV_phy, "RiskGroup2", "Species")
#'
#' prepped <- prep_recipe(rec, results = list(), bakes = list())
#' expect_s4_class(prepped, "PrepRecipe")
#' expect_true(methods::validObject(prepped, test = TRUE))
#'
#' rec@tax_info <- tibble::tibble(tax_lev = "Missing_rank")
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
    var <- NULL
    if (nrow(get_var(object)) > 0) {
      var <- get_var(object) %>% dplyr::pull(vars)
    }
    
    if (is.null(var)) {
      cat(
        glue::glue(
          "     {cross()} undefined variable of interest. Use ", 
          "{crayon::bgMagenta('add_var()')} to add it to Recipe!"
        ),
        "\n"
      )
    } else {
      var <- get_var(object) %>% dplyr::pull(1)
      var_vals <- sample_data(object) %>% dplyr::pull(.env$var)
      if (is.character(var_vals) | is.factor(var_vals)) {
        levs <-
          factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
        msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
      }
      if (is.numeric(var_vals)) {
        msg <- glue::glue("class: numeric")
      }
      cat(glue::glue(
        "     {info()} variable of interes {crayon::blue(var)} ({msg})"),
        "\n")
    }
    
    ## Taxa
    var <- NULL
    if (nrow(get_tax(object)) > 0) {
      var <- get_tax(object) %>% dplyr::pull(tax_lev)
    }
    
    if (is.null(var)) {
      cat(
        glue::glue(
          "     {cross()} undefined taxonomic level. Use ", 
          "{crayon::bgMagenta('add_tax()')} to add it to Recipe!"
        ),
        "\n"
      )
    } else {
      cat(glue::glue(
        "     {info()} taxonomic level {crayon::blue(get_tax(object))}"
      ),
      "\n\n")
    }
    
    ## Results
    cat("Results:\n\n")
    names(object@results) %>%
      purrr::discard(stringr::str_detect(., "step_subster|step_filter")) %>%
      purrr::walk( ~ {
        n_taxa <-
          object@results[[.x]][[1]] %>%
          dplyr::filter(signif == TRUE) %>% 
          dplyr::pull(taxa_id) %>%
          unique() %>%
          length()
        
        n_taxa <- crayon::silver(glue::glue("diff_taxa = {n_taxa}"))
        
        cat(c(glue::glue("     {tick()} {.x} {n_taxa}"), "\n"))
      })

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
        dplyr::filter(
          sum_methods == length(steps_ids(object, type = "da"))
        ) %>%
        nrow()
      
      cli::cat_line()
      overlap_unit <- if (is.null(get_model(object))) "taxa" else
        "taxon-contrast effects"
      cat(glue::glue(
        "     {info()} {n_overlap} {overlap_unit} are present in all tested methods"),
        "\n\n"
      )
    }
    
    ## Bakes
    if (length(object@bakes) > 0) {
      cat("Bakes:\n\n")
      
      object@bakes %>%
        purrr::iwalk(~ {
          msg <-
            .x %>%
            purrr::map2_chr(names(.), ~ {
              if (is.null(.x)) {
                .x <- "NULL"
              }
              if (.y == "weights" & .x != "NULL") {
                .x <- 
                  glue::glue("{names(.x)} = {.x}") %>%
                  stringr::str_c(collapse = ", ") %>%
                  stringr::str_c("c(", ., ")")
              }
              glue::glue("{.y}: {.x}")
            }) %>% stringr::str_c(collapse = ", ")
          
          cat(c(
            glue::glue(
              "     {dot()} {crayon::blue(crayon::bold(paste0(.y, ' ->')))} ", 
              "{crayon::silver(msg)}"
            ),
            "\n"
          ))
        })
    }
  }
)
