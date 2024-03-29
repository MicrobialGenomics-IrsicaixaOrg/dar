# CLASS RECIPE -----------------------------------------------------------------

## class unions ----

#' @rdname Recipe-class
methods::setClassUnion("phyloseq_or_null", c("phyloseq", "NULL"))

#' @rdname Recipe-class
methods::setClassUnion("tibble_or_NULL", c("tbl_df", "NULL"))

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
    steps = "list"
  ),
  prototype = list(
    phyloseq = NULL,
    var_info = NULL,
    tax_info = NULL,
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
#' colnames(metaHIV_phy@tax_table) <-
#'   c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Sp")
#'   
#' expect_error( 
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species") 
#' )
#' 
#' data(GlobalPatterns, package = "mia")
#' expect_s4_class(recipe(GlobalPatterns), "Recipe")
#' 
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))") |>
#'   step_metagenomeseq(rm_zeros = 0.01) |>
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
#' rec <- step_ancom(rec)
#' rec <- step_ancom(rec)
#'
#' ## The same with bake
#' da_results <- bake(da_results)
#' da_results <- bake(da_results)
recipe <- function(microbiome_object = NULL, 
                   var_info = NULL, 
                   tax_info = NULL, 
                   steps = list()) {

  var_info <- tibble::tibble(vars = var_info)
  tax_info <- tibble::tibble(tax_lev = tax_info)

  expected <- c("Kingdom", "Phylum", "Class", "Order", 
                "Family", "Genus", "Species")
  
  if (!is(microbiome_object, "phyloseq") && 
      !is(microbiome_object, "TreeSummarizedExperiment")) {
    rlang::abort(
      c(
        "The input object must be a phyloseq or TreeSummarizedExperiment object."
      )
    )
  }
  
  if (is(microbiome_object, "TreeSummarizedExperiment")) {
    microbiome_object <- 
      mia::makePhyloseqFromTreeSummarizedExperiment(microbiome_object)
  }
  
  if (!all(colnames(microbiome_object@tax_table) %in% expected)) {
    rlang::abort(
      c(
        glue::glue(
          "'rank' must be a value from {crayon::bgMagenta('taxonomyRanks()')}."
        ),
        glue::glue(
          "Rename the columns from the tax_table slot of your input phyloseq ",
          "with standard names.",
        ),
        glue::glue(
          "Standard names: {stringr::str_c(expected, collapse = ', ')}."
        )
      )
    )
  }
  
  methods::new(
    Class = "Recipe",
    phyloseq = microbiome_object,
    var_info = var_info,
    tax_info = tax_info,
    steps = steps
  )
}

## validity ----

methods::setValidity(
  Class = "Recipe",
  method = function(object) {
    TRUE
  }
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

# METHODS RECIPE ---------------------------------------------------------------

## get_var ----

#' Returns var_info from Recipe-class object
#'
#' @param rec A `Recipe` object
#'
#' @aliases get_var
#' @return Tibble containing `var_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract variable of interest
#' get_var(rec)
methods::setGeneric("get_var", function(rec) standardGeneric("get_var"))

#' @rdname get_var
#' @export
methods::setMethod(
  f = "get_var",
  signature = "Recipe",
  definition = function(rec) { rec@var_info }
)

## get_tax ----

#' Returns tax_info from Recipe-class object
#'
#' @param rec A `Recipe` object
#'
#' @aliases get_tax
#' @return Tibble containing `tax_info`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract taxonomic level
#' get_tax(rec)
methods::setGeneric("get_tax", function(rec) standardGeneric("get_tax"))

#' @rdname get_tax
#' @export
#' @autoglobal
methods::setMethod(
  f = "get_tax",
  signature = "Recipe",
  definition = function(rec) { rec@tax_info }
)

## get_phy ----

#' Returns phyloseq from Recipe-class object
#'
#' @param rec A `Recipe` object
#'
#' @aliases get_phy
#' @return Phyloseq class object
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract phyloseq object
#' get_phy(rec)
methods::setGeneric("get_phy", function(rec) standardGeneric("get_phy"))

#' @rdname get_phy
#' @export
#' @autoglobal
methods::setMethod(
  f = "get_phy",
  signature = "Recipe",
  definition = function(rec) { rec@phyloseq }
)

## add_var ----

#' Adds variable of interest to the Recipe
#'
#' @param rec A `Recipe` object.
#' @param var_info A character string of column names corresponding to variables
#'   that will be used in any context.
#'
#' @aliases add_var
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <- recipe(metaHIV_phy)
#'
#' ## add var info
#' rec <- add_var(rec, var_info = "RiskGroup2")
#' rec
#'
#' ## add var info to a prep-recipe returns an error
#' data(test_prep_rec)
#' err <- testthat::expect_error(
#'   add_var(test_prep_rec, var_info = "RiskGroup2")
#' )
#' 
#' err
methods::setGeneric("add_var", function(rec, var_info)
  standardGeneric("add_var"))

#' @rdname add_var
#' @export
#' @autoglobal
methods::setMethod(
  f = "add_var",
  signature = "Recipe",
  definition = function(rec, var_info) {
    rec@var_info <- tibble::tibble(vars = var_info)
    rec
  }
)

## add_tax ----

#' Adds taxonomic level of interest in the Recipe.
#'
#' @param rec A `Recipe` object.
#' @param tax_info A character string of taxonomic levels that will be used in
#'   any context.
#'
#' @aliases add_tax
#' @return A `Recipe` object.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy)
#'
#' ## add var info
#' rec <- add_tax(rec, tax_info = "Species")
#' rec
#'
#' ## add tax info to a prep-Recipe returns an error
#' data(test_prep_rec)
#' err <- testthat::expect_error(
#'   add_tax(test_prep_rec, tax_info = "Species")
#' )
#' 
#' err
methods::setGeneric("add_tax", function(rec, tax_info)
  standardGeneric("add_tax"))

#' @rdname add_tax
#' @export
#' @autoglobal
methods::setMethod(
  f = "add_tax",
  signature = "Recipe",
  definition = function(rec, tax_info) {
    rec@tax_info <- tibble::tibble(vars = tax_info)
    rec
  }
)



## Phyloseq slots as tibble from Recipe ----

#' Extracts tax_table from phyloseq inside a Recipe
#'
#' @param rec A Recipe or Recipe step.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract tax_table from phyloseq object
#' tax_table(rec)
methods::setGeneric("tax_table", function(rec) standardGeneric("tax_table"))

#' @rdname tax_table
#' @export
#' @autoglobal
methods::setMethod(
  f = "tax_table",
  signature = "Recipe",
  definition = function(rec) {
    rec@phyloseq@tax_table %>%
      to_tibble("taxa_id") %>%
      dplyr::select(taxa_id, taxa = !!get_tax(rec)[[1]])
  }
)

#' Extracts sample_data from phyloseq inside a Recipe
#'
#' @param rec A Recipe or Recipe step.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract sample_data from phyloseq object
#' sample_data(rec)
methods::setGeneric("sample_data", function(rec) standardGeneric("sample_data"))

#' @rdname sample_data
#' @export
#' @autoglobal
methods::setMethod(
  f = "sample_data",
  signature = "Recipe",
  definition = function(rec) {
    rec@phyloseq %>%
      phyloseq::sample_data() %>%
      to_tibble("sample_id") %>%
      dplyr::select(sample_id, !!get_var(rec)[[1]])
  }
)

#' Extracts otu_table from phyloseq inside a Recipe
#'
#' @param rec A Recipe or Recipe step.
#'
#' @return A tibble
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#'
#' ## Define recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Species")
#'
#' ## Extract otu_table from phyloseq object
#' otu_table(rec)
methods::setGeneric("otu_table", function(rec) standardGeneric("otu_table"))

#' @rdname otu_table
#' @export
#' @autoglobal
methods::setMethod(
  f = "otu_table",
  signature = "Recipe",
  definition = function(rec) {
    rec@phyloseq %>%
      phyloseq::otu_table() %>%
      to_tibble("taxa_id")
  }
)

# CLASS PrepRecipe ------------------------------------------------------------

## class def ----

#' PrepRecipe-class object
#'
#' A PrepRecipe is Recipe with the results corresponding to the steps defined
#' in the Recipe.
#'
#' @slot results Contains the results of all defined analysis in the Recipe.
#' @slot bakes Contains the executed bakes.
#'
#' @name PrepRecipe-class
#' @rdname PrepRecipe-class
#' @exportClass PrepRecipe
#' @return PrepRecipe-class object
#' @autoglobal
methods::setClass(
  Class = "PrepRecipe",
  contains = "Recipe",
  slots = c(results = "list", bakes = "list")
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
prep_recipe <- function(rec, results, bakes) {
  methods::new(
    Class = "PrepRecipe",
    results = results,
    bakes = bakes,
    rec
  )
}

## validity ----

methods::setValidity(
  Class = "Recipe",
  method = function(object) {
    TRUE
  }
)

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
    
    if (length(object@results) > 0) {
      n_overlap <-
        find_intersections(object) %>%
        dplyr::filter(
          sum_methods == length(steps_ids(object, type = "da"))
        ) %>%
        nrow()
      
      cli::cat_line()
      cat(glue::glue(
        "     {info()} {n_overlap} taxa are present in all tested methods"),
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


# METHODS PREP_RECIPE ----------------------------------------------------------

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_prep <- function(x, ...) {
  c("furrr", "future")
}


## add_var ----

#' @rdname add_var
#' @export
#' @autoglobal
methods::setMethod(
  f = "add_var",
  signature = "PrepRecipe",
  definition = function(rec, var_info) {
    rlang::abort("var_info can only be added to a non-PrepRecipe")
  }
)

## add_tax ----

#' @rdname add_tax
#' @export
#' @autoglobal
methods::setMethod(
  f = "add_tax",
  signature = "PrepRecipe",
  definition = function(rec, tax_info) {
    rlang::abort("tax_info can only be added to a non-PrepRecipe")
  }
)

## prep ----

#' Performs all the steps defined in a Recipe
#'
#' For a Recipe with at least one preprocessing or DA operation run the steps in
#' a convenient order.
#'
#' @param rec A `Recipe` object. and furrr packages.
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
#'
#' ## Define Recipe
#' rec <-
#'   recipe(metaHIV_phy, var_info = "RiskGroup2", tax_info = "Class") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.03 * length(x))") |>
#'   step_maaslin()
#'
#' ## Prep Recipe
#' da_results <- prep(rec)
#'
#' ## If you try
#'
#' ## Consensus strategy
#' n_methods <- 2
#' da_results <- bake(da_results, count_cutoff = n_methods)
#' da_results
#'
#' ## If you try to run prep on an object of class PrepRecipe it returns an 
#' ## error.
#' err <- testthat::expect_error(prep(da_results))
#' err
#' 
#' ## You can force the overwrite with:
#' prep(rec, force = TRUE)
#'
#' ## This function can operate in parallel thanks to future and furrr packages.
#' prep(rec, parallel = TRUE, workers = 2)
methods::setGeneric(
  name = "prep",
  def = function(rec,
                 parallel = TRUE,
                 workers = 4,
                 force = FALSE) {
    standardGeneric("prep")
  }
)

#' @rdname prep
#' @export
methods::setMethod(
  f = "prep",
  signature = "Recipe",
  definition = function(rec, parallel, workers, force) {
    if ("results" %in% methods::slotNames(rec) & !force) {
      rlang::abort(c(
        "The input Recipe has already been prep!",
        i = glue::glue(
          "To force the rerun of all steps plese run ", 
          "{crayon::bgMagenta('prep(rec, force = T)')}"
        )
      ))
    }
    
    check <- utils::capture.output(required_deps(rec))
    if (length(check) > 0) {
      rlang::abort(c(
        "Not all necessary dependencies are installed.",
        i = glue::glue(
          "Use {crayon::bgMagenta('required_deps(rec)')} to see how to ", 
          "install them."
        )
      ))
    }
    
    ## Phyloseq preporcessing steps
    to_execute <-
      rec@steps %>%
      purrr::map_chr(step_to_expr) %>%
      purrr::keep(stringr::str_detect(., "run_subset|run_filter"))
    
    to_execute <-
      rec@steps %>%
      purrr::map_chr(step_to_expr) %>%
      purrr::keep(stringr::str_detect(., "run_rarefaction")) %>%
      c(to_execute, .)
    
    for (.x in to_execute) {
      rec <- eval(parse(text = .x))
    }
    
    ## DA steps
    names <-
      purrr::map_chr(rec@steps, ~ .x[["id"]]) %>%
      purrr::discard(stringr::str_detect(., "subset|filter|rarefaction"))
    
    if (parallel) {
      recipes_pkg_check(required_pkgs_prep(), "prep()")
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential))
      
      res <-
        rec@steps %>%
        purrr::map_chr(step_to_expr) %>%
        rarefy_msg() %>%
        purrr::discard(
          stringr::str_detect(., "run_subset|run_filter|run_rarefaction")
        ) %>%
        furrr::future_map( ~ {
          rec <- rec
          eval(parse(text = .x))
        },
        .options = furrr::furrr_options(
          seed = TRUE,
          packages = "dar",
          conditions = structure("condition", exclude = "message")
        ))
      
      names(res) <- names
    }
    
    if (!parallel) {
      res <-
        rec@steps %>%
        purrr::map(step_to_expr) %>%
        rarefy_msg() %>%
        purrr::discard(
          stringr::str_detect(., "run_subset|run_filter|run_rarefaction")
        ) %>%
        purrr::map(~ {
          rec <- rec
          eval(parse(text = .x))
        })
      
      names(res) <- names
    }
    
    prep_recipe(rec, res, list())
  }
)

## Find intersections ----

#' Returns data.frame with OTU intersection between methods
#'
#' @param rec A `Recipe` object.
#' @param steps character vector with step_ids to take in account.
#' @param tidy Boolan indicating if result must be in tidy format.
#'
#' @aliases intersection_df
#' @return data.frame class object
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#'
#' df <- intersection_df(test_prep_rec)
#' head(df)
#'
#' ## intersection_df function needs a prep-Recipe. If you pass a a non-prep
#' ## recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(intersection_df(test_rec))
#' err
methods::setGeneric(
  name = "intersection_df",
  def = function(rec, steps = steps_ids(rec, "da"), tidy = FALSE) {
    standardGeneric("intersection_df")
  }
)

#' @rdname intersection_df
#' @export
#' @autoglobal
methods::setMethod(
  f = "intersection_df",
  signature = "Recipe",
  definition = function(rec, steps, tidy) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then retry last command."
      )
    ))
  }
)

#' @rdname intersection_df
#' @export
#' @autoglobal
methods::setMethod(
  f = "intersection_df",
  signature = "PrepRecipe",
  definition = function(rec, steps, tidy) {
    df <- names(rec@results) %>%
      purrr::keep(. %in% steps) %>%
      purrr::set_names() %>%
      purrr::map_dfc( ~ {
        taxa <- 
          rec@results[[.x]][[1]] %>% 
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
)

## Overlap df ----

#' Overlap of significant OTUs between tested methods.
#'
#' @param rec A `Recipe` object.
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
#'
#' ## Running the function returns a UpSet plot ordered by frequency.
#' df <- overlap_df(test_prep_rec, steps_ids(test_prep_rec, "da"))
#' head(df)
#'
#' ## If you want to exclude a method for the plot, you can remove it with the
#' ## step parameter. In the following example we eliminate from the graph the
#' ## results of maaslin
#' overlap_df(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#'
#' ## overlap_df function needs a prep-Recipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(overlap_df(test_rec))
#' err
methods::setGeneric(
  name = "overlap_df",
  def = function(rec, steps = steps_ids(rec, "da"), type = "all") {
    standardGeneric("overlap_df")
  }
)

#' @rdname overlap_df
#' @export
#' @autoglobal
methods::setMethod(
  f = "overlap_df",
  signature = "Recipe",
  definition = function(rec, steps, type) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('overlap_df()')}"
      )
    ))
  }
)

#' @rdname overlap_df
#' @export
#' @autoglobal
methods::setMethod(
  f = "overlap_df",
  signature = "PrepRecipe",
  definition = function(rec, steps, type) {
    df <- 
      intersection_df(rec) %>%  
      tibble::as_tibble() %>% 
      dplyr::select(dplyr::all_of(steps))
    
    if (type == "da") {
      df <- 
        dplyr::rowwise(df) %>% 
        dplyr::mutate(
          sum = sum(dplyr::across(dplyr::all_of(steps))),
          .before = 1
        ) %>% 
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
)

## Cool: Extract results from defined  <-  ----

#' Extract results from defined bake
#'
#' @param rec A `Recipe` object.
#' @param bake Name or index of the bake to extract.
#'
#' @aliases cool
#' @return tbl_df
#' @export
#' @autoglobal
#' @examples
#' data(test_prep_rec)
#'
#' ## First we need to add bakes (extraction strategies) to the PrepRecipe.
#' rec <- bake(test_prep_rec)
#'
#' ## Finally we can extract the results with the cool method
#' cool(rec)
#'
#' ## By default cool extracts the results of the first bake. If we have more
#' ## bakes we can extract the one that you want with the bake parameter.
#' rec <- bake(rec, count_cutoff = 1)
#' cool(rec, 2)
#'
#' ## bake and cool methods needs a prep-Recipe. If you pass a non-PrepRecipe
#' ## the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(cool(test_rec))
#' err
methods::setGeneric(
  name = "cool",
  def = function(rec, bake = 1) {
    standardGeneric("cool")
  }
)

#' @rdname cool
#' @export
#' @autoglobal
methods::setMethod(
  f = "cool",
  signature = "Recipe",
  definition = function(rec, bake) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('cool()')}"
      )
    ))
  }
)

#' @rdname cool
#' @export
#' @autoglobal
methods::setMethod(
  f = "cool",
  signature = "PrepRecipe",
  definition = function(rec, bake) {
    all_bakes <- rec@bakes
    all_names <- all_bakes %>% purrr::map_chr(~ as.character(.x[["id"]]))
    
    if (is.numeric(bake) & length(all_bakes) < bake) {
      rlang::abort(c(
        "Bake index is not defined in the PrepRecipe!",
        glue::glue(
          "Run {crayon::bgMagenta('bake(PrepRecipe)')} and then try with ", 
          "{crayon::bgMagenta('cool()')}"
        )
      ))
    }
    if (!is.numeric(bake) & !bake %in% all_names) {
      rlang::abort(c(
        "Bake name is not defined in the PrepRecipe!",
        glue::glue(
          "Run {crayon::bgMagenta('bake(PrepRecipe)')} and then try with ", 
          "{crayon::bgMagenta('cool()')}"
        )
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
)
