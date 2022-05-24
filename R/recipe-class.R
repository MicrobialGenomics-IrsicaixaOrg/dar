# CLASS ---------------------------------------------------------------------------

## class unions ----

#' @rdname recipe-class
methods::setClassUnion("phyloseq_or_null", c("phyloseq", "NULL"))

#' @rdname recipe-class
methods::setClassUnion("tibble_or_NULL", c("tbl_df", "NULL"))

## class def ----

#' recipe-class object
#'
#' A recipe is a description of the steps to be applied to a data set in order to prepare
#' it for data analysis.
#'
#' @slot phyloseq Phyloseq-class object
#' @slot var_info A tibble that contains the current set of terms in the data set. This
#'   initially defaults to the same data contained in `var_info`.
#' @slot tax_info A tibble that contains the current set of taxonomic levels that will be
#'   used in the analysis.
#' @slot steps List of step-class objects that will be used by DA.
#'
#' @name recipe-class
#' @rdname recipe-class
#' @exportClass recipe
methods::setClass(
  Class = "recipe",
  slots = c(
    phyloseq = "phyloseq_or_null",
    var_info = "tibble_or_NULL",
    tax_info = "tibble_or_NULL",
    steps = "list",
    results = "list"
  ),
  prototype = list(
    phyloseq = NULL,
    var_info = NULL,
    tax_info = NULL,
    steps = list(),
    results = list()
  )
)

## constructor ----

#' Create a recipe for preprocessing data
#'
#' A recipe is a description of the steps to be applied to a data set in order to prepare
#' it for data analysis.
#'
#' @param phyloseq Phyloseq-class object.
#' @param var_info A character string of column names corresponding to variables that will
#'   be used in any context.
#' @param tax_info A character string of taxonomic levels that will be used in any
#'   context.
#'
#' @return An object of class `recipe` with sub-objects:
#'   \item{phyloseq}{object of class `phyloseq` with taxa abundance information.}
#'   \item{var_info}{A tibble that contains the current set of terms in the data set.
#'   This initially defaults to the same data contained in `var_info`.}
#'   \item{tax_info}{A tibble that contains the current set of taxonomic levels that will
#'   be used in the analysis.}
#'
#' @aliases recipe
#' @export
recipe <- function(phyloseq = NULL, var_info = NULL, tax_info = NULL) {

  var_info <- tibble::tibble(vars = var_info)
  tax_info <- tibble::tibble(tax_lev = tax_info)

  methods::new(
    Class = "recipe",
    phyloseq = phyloseq,
    var_info = var_info,
    tax_info = tax_info,
    steps = list(),
    results = list()
  )
}

## validity ----

methods::setValidity(
  Class = "recipe",
  method = function(object) {
    TRUE
  }
)

## printing ----

# methods::setMethod("show", signature = "recipe", definition = function(rec) {})


# METHODS ---------------------------------------------------------------------------

## get_var ----

#' Returns var_info from recipe-class object
#'
#' @param rec A `recipe` object
#'
#' @aliases get_var
#' @return Tibble containing `var_info`.
#' @export
methods::setGeneric("get_var", function(rec) standardGeneric("get_var"))

#' @rdname get_var
#' @export
methods::setMethod(
  f = "get_var",
  signature = "recipe",
  definition = function(rec) { rec@var_info }
)

## get_tax ----

#' Returns tax_info from recipe-class object
#'
#' @param rec A `recipe` object
#'
#' @aliases get_tax
#' @return Tibble containing `tax_info`.
#' @export
methods::setGeneric("get_tax", function(rec) standardGeneric("get_tax"))

#' @rdname get_tax
#' @export
methods::setMethod(
  f = "get_tax",
  signature = "recipe",
  definition = function(rec) { rec@tax_info }
)

## get_phy ----

#' Returns phyloseq from recipe-class object
#'
#' @param rec A `recipe` object
#'
#' @aliases get_phy
#' @return Phyloseq class object
#' @export
methods::setGeneric("get_phy", function(rec) standardGeneric("get_phy"))

#' @rdname get_phy
#' @export
methods::setMethod(
  f = "get_phy",
  signature = "recipe",
  definition = function(rec) { rec@phyloseq }
)

## package deps ----

#' Methods for tracking which additional packages are needed for steps.
#'
#' @param rec A recipe or recipe step.
#'
#' @return A character vector
#' @export
methods::setGeneric("required_pkgs_recipe", function(rec) standardGeneric("required_pkgs_recipe"))

#' @rdname required_pkgs_recipe
#' @export
methods::setMethod(
  f = "required_pkgs_recipe",
  signature = "recipe",
  definition = function(rec) {
    res <- purrr::map(rec@steps, required_pkgs)
    res <- unique(unlist(res))
    res <- res[length(res) != 0]
    res
  }
)

## Phyloseq slots as tibble from recipe ----

#' Extracts tax_table from phyloseq inside a recipe
#'
#' @param rec A recipe or recipe step.
#'
#' @return A tibble
#' @export
methods::setGeneric("tax_table", function(rec) standardGeneric("tax_table"))

#' @rdname tax_table
#' @export
methods::setMethod(
  f = "tax_table",
  signature = "recipe",
  definition = function(rec) {
    rec@phyloseq@tax_table %>%
      to_tibble("taxa_id") %>%
      dplyr::select(taxa_id, taxa = !!get_tax(rec)[[1]])
  }
)

#' Extracts sample_data from phyloseq inside a recipe
#'
#' @param rec A recipe or recipe step.
#'
#' @return A tibble
#' @export
methods::setGeneric("sample_data", function(rec) standardGeneric("sample_data"))

#' @rdname sample_data
#' @export
methods::setMethod(
  f = "sample_data",
  signature = "recipe",
  definition = function(rec) {
    rec@phyloseq %>%
      phyloseq::sample_data() %>%
      to_tibble("sample_id") %>%
      dplyr::select(sample_id, !!get_var(rec)[[1]])
  }
)

#' Extracts otu_table from phyloseq inside a recipe
#'
#' @param rec A recipe or recipe step.
#'
#' @return A tibble
#' @export
methods::setGeneric("otu_table", function(rec) standardGeneric("otu_table"))

#' @rdname otu_table
#' @export
methods::setMethod(
  f = "otu_table",
  signature = "recipe",
  definition = function(rec) {
    rec@phyloseq %>%
      phyloseq::otu_table() %>%
      to_tibble("taxa_id")
  }
)


#' @noRd
#' @keywords internal
required_pkgs_backe <- function(x, ...) {  c("furrr", "future") }


#' Estimate a preprocessing recipe
#'
#' For a recipe with at least one preprocessing operation, estimate the required
#'   parameters from a training set that can be later applied to other data
#'   sets.
#'
#' @param rec
#' @param parallel
#' @param workers
#'
#' @export
prep <- function(rec, parallel = TRUE, workers = 8) {

  ## Phyloseq preporcessing steps
  to_execute <-
    rec@steps %>%
    purrr::map_chr(step_to_expr) %>%
    purrr::keep(stringr::str_detect(., "run_subset|run_filter"))

  for (.x in to_execute) {
    rec <- eval(parse(text = .x))
  }

  ## DA steps
  names <-
    purrr::map_chr(rec@steps, ~ .x[["id"]]) %>%
    purrr::discard(stringr::str_detect(., "subset|filter"))

  if (parallel) {
    recipes_pkg_check(required_pkgs_backe())
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential))

    res <-
      rec@steps %>%
      purrr::map_chr(step_to_expr) %>%
      purrr::discard(stringr::str_detect(., "run_subset|run_filter")) %>%
      furrr::future_map(~ {
        rec <- rec
        eval(parse(text = .x))
      }, .options = furrr::furrr_options(seed = TRUE))

    names(res) <- names
  }

  if (!parallel) {
    res <-
      rec@steps %>%
      purrr::map(step_to_expr) %>%
      purrr::discard(stringr::str_detect(., "run_subset|run_filter")) %>%
      purrr::map(~ {
        rec <- rec
        eval(parse(text = .x))
      })

    names(res) <- names
  }

  rec@results <- res

  ## Otu Overlaps
  rec@results$intersections <- find_intersections(rec)

  rec
}
