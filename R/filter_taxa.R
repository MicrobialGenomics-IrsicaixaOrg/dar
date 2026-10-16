#' Filter taxa based on across-sample OTU abundance criteria
#'
#' Add a preprocessing step that applies one safely serializable predicate to
#' each taxon's abundance vector across samples. Filtering is performed during
#' [prep()]; taxa for which the predicate returns `FALSE` are removed.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param .f A function with exactly one argument named `x`, without a default,
#'   or a character string containing that function expression. It must return
#'   a single logical value for the abundance vector `x`. Only the restricted
#'   expression syntax described below is accepted; lists of functions are not
#'   supported.
#' @param id A character string that is unique to this step to identify it.
#'
#' @details The predicate body must be a single expression, without braces,
#'   assignments, external variables or user-defined helper calls. Write
#'   thresholds as literal numbers inside the expression, even when using the
#'   function directly in R. These restrictions keep recipes safely exportable
#'   and importable; they are checked when the step is added, not just on export.
#'
#'   Accepted expressions use `x`, finite numeric or non-missing logical
#'   literals, parentheses, arithmetic operators (`+`, `-`, `*`, `/`, `^`,
#'   `%%`, `%/%`), comparisons (`>`, `>=`, `<`, `<=`, `==`, `!=`), logical
#'   operators (`&`, `|`, `!`), and calls to `sum`, `length`, `mean`, `median`,
#'   `min`, `max`, `var`, `sd`, `any`, `all` or `is.na`. The only supported
#'   named argument is `na.rm` in `sum`, `mean`, `median`, `min`, `max`, `var`,
#'   `sd`, `any` and `all`. Namespace-qualified calls are not accepted.
#'
#' @include recipe-class.R
#' @family filter phy steps
#' @aliases step_filter_taxa
#' @return A `Recipe` with the filtering step added. The constructor does not
#'   filter the data or return a logical vector. During [prep()], the step uses
#'   [phyloseq::filter_taxa()] with `prune = TRUE`; there is no `prune` argument
#'   in this wrapper.
#' @export
#' @autoglobal
#' @tests 
#' data(test_prep_rec)
#' expect_error(step_filter_taxa(test_prep_rec))
#' @examples
#' data(metaHIV_phy)
#'
#' ## Keep taxa present in at least 3% of samples
#' recipe(metaHIV_phy) |>
#'   step_filter_taxa(.f = function(x) sum(x > 0) >= 0.03 * length(x)) |>
#'   prep()
#'
#' ## The same predicate can be supplied as a string
#' recipe(metaHIV_phy) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= 0.03 * length(x)")
step_filter_taxa <- function(rec, .f, id = rand_id("filter_taxa")) {
  check_recipe(rec)
  recipes_pkg_check(required_pkgs_filter_taxa(), "step_filter_taxa()")
  .f <- normalize_filter_predicate(.f)

  add_step(rec, step(subclass = "filter_taxa", .f = .f, id = id))
}

#' @noRd
#' @autoglobal
#' @keywords internal
run_filter_taxa <- function(rec, .f, id) {
  rec@phyloseq <- phyloseq::filter_taxa(get_phy(rec), .f, prune = TRUE)

  rec
}

#' @noRd
#' @keywords internal
required_pkgs_filter_taxa <- function(x, ...) {  c("bioc::phyloseq") }
