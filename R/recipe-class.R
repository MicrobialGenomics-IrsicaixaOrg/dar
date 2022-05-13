methods::setClassUnion("phyloseq_or_null", c("phyloseq", "NULL"))
methods::setClassUnion("char_or_NULL", c("character", "NULL"))

#' Title
#'
#' @slot phyloseq
#' @slot var_info
#' @slot tax_info
#' @slot steps
#'
#' @name recipe-class
#' @rdname recipe-class
#' @exportClass recipe
methods::setClass(
  Class = "recipe",
  slots = c(
    phyloseq = "phyloseq_or_null",
    var_info = "char_or_NULL",
    tax_info = "char_or_NULL",
    steps = "list"
  ),
  prototype = list(
    phyloseq = NULL,
    var_info = NULL,
    tax_info = NULL,
    steps = list()
  )
)

recipe <- function(phyloseq = NULL, var_info = NULL, tax_info = NULL, steps = list()) {
  methods::new(
    Class = "recipe",
    phyloseq = phyloseq,
    var_info = var_info,
    tax_info = tax_info,
    steps = steps
  )
}

methods::setValidity(
  Class = "recipe",
  method = function(rec) {
    TRUE
  }
)


methods::setGeneric("get_var", function(rec) standardGeneric("get_var"))
methods::setMethod(
  f = "get_var",
  signature = "recipe",
  definition = function(rec) { rec@var_info }
)

methods::setGeneric("get_tax", function(rec) standardGeneric("get_tax"))
methods::setMethod(
  f = "get_tax",
  signature = "recipe",
  definition = function(rec) { rec@tax_info }
)


