# CLASS RECIPE -------------------------------------------------------------------------

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
#' @param steps list with steps.
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
recipe <- function(phyloseq = NULL, var_info = NULL, tax_info = NULL, steps = list()) {

  var_info <- tibble::tibble(vars = var_info)
  tax_info <- tibble::tibble(tax_lev = tax_info)

  methods::new(
    Class = "recipe",
    phyloseq = phyloseq,
    var_info = var_info,
    tax_info = tax_info,
    steps = steps
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

methods::setMethod("show", signature = "recipe", definition = function(object) {
  cli::cat_rule(crayon::blue("DAR Recipe"))
  cat("Inputs:\n\n")

  ## Phyloseq
  phy <- get_phy(object)
  ntax <- phyloseq::ntaxa(phy)
  nsam <- phyloseq::nsamples(phy)
  cat(glue::glue("     {info()} phyloseq object with {crayon::blue(ntax)} taxa and {crayon::blue(nsam)} samples"), "\n")

  ## Variable
  var <- get_var(object) %>% dplyr::pull(1)
  var_vals <- sample_data(object) %>% dplyr::pull(.env$var)
  if (is.character(var_vals) | is.factor(var_vals)) {
    levs <- factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
    msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
  }
  if (is.numeric(var_vals)) {
    msg <- glue::glue("class: numeric")
  }
  cat(glue::glue("     {info()} variable of interes {crayon::blue(var)} ({msg})"), "\n")

  ## Taxa
  cat(glue::glue("     {info()} taxonomic level {crayon::blue(get_tax(object))}"), "\n\n")

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

# METHODS RECIPE -----------------------------------------------------------------------

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

# CLASS PREP_RECIPE -------------------------------------------------------------------------


## class def ----

#' prep_recipe-class object
#'
#' A prep_recipe is recipe with the results corresponding to the steps defined in the
#' recipe.
#'
#' @slot results Contains the results of all defined analysis in the recipe.
#' @slot bakes Contains the executed bakes.
#'
#' @name prep_recipe-class
#' @rdname prep_recipe-class
#' @exportClass prep_recipe
methods::setClass(
  Class = "prep_recipe",
  contains = "recipe",
  slots = c(results = "list", bakes = "list")
)

## constructor ----

#' Create a recipe prep_recipe.
#'
#' A prep_recipe is recipe with the results corresponding to the steps defined in the
#' recipe.
#'
#' @param rec A recipe object.
#' @param results list with the results
#' @param bakes list with saved bakes
#'
#' @return An object of class `prep_recipe`.
#'
#' @aliases prep_recipe
#' @export
prep_recipe <- function(rec, results, bakes) {
  methods::new(
    Class = "prep_recipe",
    results = results,
    bakes = bakes,
    rec
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

methods::setMethod("show", signature = "prep_recipe", definition = function(object) {
  cli::cat_rule(crayon::blue("DAR Results"))
  cat("Inputs:\n\n")

  ## Phyloseq
  phy <- get_phy(object)
  ntax <- phyloseq::ntaxa(phy)
  nsam <- phyloseq::nsamples(phy)
  cat(glue::glue("     {info()} phyloseq object with {crayon::blue(ntax)} taxa and {crayon::blue(nsam)} samples"), "\n")

  ## Variable
  var <- get_var(object) %>% dplyr::pull(1)
  var_vals <- sample_data(object) %>% dplyr::pull(.env$var)
  if (is.character(var_vals) | is.factor(var_vals)) {
    levs <- factor(var_vals) %>% levels() %>% stringr::str_c(collapse = ", ")
    msg <- glue::glue("class: {class(var_vals)}, levels: {levs}")
  }
  if (is.numeric(var_vals)) {
    msg <- glue::glue("class: numeric")
  }
  cat(glue::glue("     {info()} variable of interes {crayon::blue(var)} ({msg})"), "\n")

  ## Taxa
  cat(glue::glue("     {info()} taxonomic level {crayon::blue(get_tax(object))}"), "\n\n")


  ## Results
  cat("Results:\n\n")
  names(object@results) %>%
    purrr::discard(stringr::str_detect(., "step_subster|step_filter")) %>%
    purrr::walk(~ {
      n_taxa <-
        object@results[[.x]][[1]] %>%
        dplyr::pull(.data$taxa_id) %>%
        unique() %>%
        length()

      n_taxa <- crayon::silver(glue::glue("diff_taxa = {n_taxa}"))

      cat(c(glue::glue("     {tick()} {.x} {n_taxa}"), "\n"))
    })

  n_overlap <-
    find_intersections(object) %>%
    dplyr::filter(.data$sum_methods == length(steps_ids(object, type = "da"))) %>%
    nrow()

  cli::cat_line()
  cat(glue::glue("     {info()} {n_overlap} taxa are present in all tested methods"), "\n")
})


# METHODS PREP_RECIPE -----------------------------------------------------------------------

#' @noRd
#' @keywords internal
required_pkgs_prep <- function(x, ...) {  c("furrr", "future") }

## prep ----

#' Performs all the steps defined in a recipe
#'
#' For a recipe with at least one preprocessing or DA operation run the steps in a
#' convenient order.
#'
#' @param rec A `recipe` object.
#' @param parallel if FALSE, no palatalization. if TRUE, parallel execution using future
#'   and furrr packages.
#' @param workers Number of workers for palatalization.
#' @param force Force the reexecution of all steps. This remove previous results.
#'
#' @aliases prep
#' @return A `prep_recipe` object.
#' @export
methods::setGeneric(
  name = "prep",
  def = function(rec, parallel = TRUE, workers = 8, force = FALSE) {
    standardGeneric("prep")
  }
)

#' @rdname prep
#' @export
methods::setMethod(
  f = "prep",
  signature = "recipe",
  definition = function(rec, parallel, workers, force) {

    if ("results" %in% methods::slotNames(rec) & !force) {
      rlang::abort(c(
        "The input recipe has already been prep!",
        i = glue::glue(
          "To force the rerun of all steps plese run {crayon::bgMagenta('prep(rec, force = T)')}"
        )
      ))
    }

    check <- utils::capture.output(required_deps(rec))
    if (length(check) > 0) {
      rlang::abort(c(
        "Not all necessary dependencies are installed.",
        i = glue::glue(
          "Use {crayon::bgMagenta('required_deps(rec)')} to see how to install them."
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
        purrr::discard(stringr::str_detect(., "run_subset|run_filter|run_rarefaction")) %>%
        furrr::future_map( ~ {
          rec <- rec
          eval(parse(text = .x))
        }, .options = furrr::furrr_options(seed = TRUE))

      names(res) <- names
    }

    if (!parallel) {
      res <-
        rec@steps %>%
        purrr::map(step_to_expr) %>%
        purrr::discard(stringr::str_detect(., "run_subset|run_filter|run_rarefaction")) %>%
        purrr::map( ~ {
          rec <- rec
          eval(parse(text = .x))
        })

      names(res) <- names
    }

    prep_recipe(rec, res, list())
  }
)

## bake ----

#' Extract results from recipe
#'
#' For a prep recipe with extracts the results.
#'
#' @param rec A `recipe` object.
#' @param count_cutoff Indicates the minimum number of methods in which an OTU must be
#'   present (Default: NULL). If count_cutoff is NULL count_cutoff is equal to
#'   `length(steps_ids(rec, "da")) - length(exclude)`
#' @param weights Named vector with the ponderation value for each method.
#' @param exclude Method ids to exclude.
#'
#' @aliases bake
#' @return Tibble containing differentially expressed taxa.
#' @export
methods::setGeneric(
  name = "bake",
  def = function(rec = rec,
                 count_cutoff = NULL,
                 weights = NULL,
                 exclude = NULL) {
    standardGeneric("bake")
  }
)

#' @rdname bake
#' @export
methods::setMethod(
  f = "bake",
  signature = "prep_recipe",
  definition = function(rec, count_cutoff, weights, exclude) {

    ids <- steps_ids(rec, type = "da") %>% .[!. %in% exclude]
    if (is.null(count_cutoff)) {
      count_cutoff <- length(steps_ids(rec, "da")) - length(exclude)
    }

    if (is.null(weights)) {
      ids %>%
        purrr::set_names() %>%
        purrr::map_chr(~ 1) %>%
        as.numeric()

      weights <- rep(1, length(ids))
      names(weights) <- ids
    }

    not_weights <- ids[!ids %in% names(weights)]
    if (length(not_weights) > 0) {
      not_weights <- not_weights %>%  stringr::str_c(collapse = ", ")
      rlang::abort(c(
        glue::glue("Some non-excluded methods are not present in the weights vector ({crayon::yellow(not_weights)})."),
        "Please include a value for this/these step/s in the vector of weights.",
        "Alternatively, explicitly exclude the method/s via the exclude parameter."
      ))
    }

    # weighted matrix
    res <-
      ids %>%
      purrr::map2_dfc(seq_along(.), ~ {
        intersection_df(rec, ids) %>%
          tibble::as_tibble() %>%
          dplyr::select(!!.x) %>%
          dplyr::mutate(!!(.x) := (!!dplyr::sym(.x) * weights[.y]))
      }) %>%
      dplyr::bind_cols(intersection_df(rec, ids)[1], .) %>%
      tidyr::pivot_longer(cols = -1) %>%
      dplyr::group_by(taxa_id) %>%
      dplyr::summarise(
        step_ids = purrr::map_chr(.data$name, ~ .x) %>% stringr::str_c(collapse = ", "),
        sum_methods = sum(.data$value)
      ) %>%
      dplyr::right_join(tax_table(rec), ., by = "taxa_id") %>%
      dplyr::arrange(-.data$sum_methods) %>%
      dplyr::filter(.data$sum_methods >= count_cutoff) %>%
      dplyr::select(.data$taxa_id, .data$taxa)

    not_incl <- NULL
    if (!is.null(exclude)) {
      coll <- stringr::str_c(exclude, collapse = ", ")
      not_incl <- glue::glue("Results from {crayon::blue(coll)} are excluded")
    }

    rlang::inform(
      c(
        i = glue::glue("Bake for {crayon::blue('count_cutoff =')} {crayon::blue(count_cutoff)}"),
        i = not_incl,
        ""
      )
    )

    res
  }
)

#' @rdname bake
#' @export
methods::setMethod(
  f = "bake",
  signature = "recipe",
  definition = function(rec, count_cutoff, weights, exclude) {
    rlang::abort(c(
      "This function needs a prep recipe!",
      glue::glue("Run {crayon::bgMagenta('prep(rec)')} and then try with {crayon::bgMagenta('bake()')}")
    ))
  }
)

## Find intersections ----

#' Returns data.frame with OTU intersection between methods
#'
#' @param rec A `recipe` object.
#' @param steps character vector with step_ids to take in account.
#'
#' @aliases intersection_df
#' @return data.frame class object
#' @export
methods::setGeneric(
  name = "intersection_df",
  def = function(rec, steps = steps_ids(rec, "da")) { standardGeneric("intersection_df") }
)

#' @rdname intersection_df
#' @export
methods::setMethod(
  f = "intersection_df",
  signature = "recipe",
  definition = function(rec, steps) {
    rlang::abort(c(
      "This function needs a prep recipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then retry last command."
      )
    ))
  }
)

#' @rdname intersection_df
#' @export
methods::setMethod(
  f = "intersection_df",
  signature = "prep_recipe",
  definition = function(rec, steps) {
    names(rec@results) %>%
      purrr::keep(. %in% steps) %>%
      purrr::set_names() %>%
      purrr::map_dfc(~ {
        taxa <- rec@results[[.x]][[1]][["taxa_id"]]
        rownames(rec@phyloseq@otu_table) %>%
          tibble::tibble(taxa_id = .) %>%
          dplyr::mutate(!!.x := dplyr::if_else(taxa_id %in% taxa, 1, 0)) %>%
          dplyr::select(!!.x)
      }) %>%
      dplyr::mutate(taxa_id = rownames(rec@phyloseq@otu_table), .before = 1) %>%
      as.data.frame()
  }
)

## Intersection plot----

#' Plot results using UpSet plot
#'
#' @param rec A `recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param ordered_by How the intersections in the matrix should be ordered by. Options
#'   include frequency (entered as "freq"), degree, or both in any order.
#'
#' @aliases intersection_plt
#' @return UpSet plot
#' @export
methods::setGeneric(
  name = "intersection_plt",
  def = function(rec, steps = steps_ids(rec, "da"), ordered_by = c("freq", "degree")) {
    standardGeneric("intersection_plt")
  }
)

#' @rdname intersection_plt
#' @export
methods::setMethod(
  f = "intersection_plt",
  signature = "recipe",
  definition = function(rec, steps) {
    rlang::abort(c(
      "This function needs a prep recipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with {crayon::bgMagenta('intersection_plt()')}"
      )
    ))
  }
)

#' @rdname intersection_plt
#' @export
methods::setMethod(
  f = "intersection_plt",
  signature = "prep_recipe",
  definition = function(rec, steps, ordered_by) {
    UpSetR::upset(
      data = intersection_df(rec, steps),
      nsets = length(rec@results),
      sets.bar.color = "#56B4E9",
      order.by = ordered_by
    )
  }
)
