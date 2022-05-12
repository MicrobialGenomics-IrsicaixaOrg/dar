methods::setClassUnion("null_phyloseq", c("NULL", "phyloseq"))
methods::setClassUnion("null_recipes", c("NULL", "list"))
methods::setClassUnion("null_da_results", c("NULL", "list"))

#' Title
#'
#' @slot phyloseq null_phyloseq.
#' @slot recipes null_recipes.
#' @slot da_results null_da_results.
#'
#' @name DarSet-class
#' @rdname DarSet-class
#' @exportClass DarSet
methods::setClass(
  Class = "DarSet",
  slots = c(
    phyloseq = "null_phyloseq",
    recipes = "null_recipe",
    da_results = "null_da_results"
  ),
  prototype = list(
    phyloseq = NULL,
    recipes = NULL,
    da_results = NULL
  )
)

DarSet <- function(phyloseq = NULL, recipes = NULL, da_results = NULL) {
  methods::new(
    Class = "DarSet",
    phyloseq = phyloseq,
    recipes = recipes,
    da_results = da_results
  )
}

methods::setValidity(
  Class = "DarSet",
  method = function(object) {
    TRUE
  }
)
