#' @include recipe-class.R
methods::setGeneric("step_ancom", function(rec, out_cut = 0, zero_cut = 1, lib_cut = 0, neg_lb = FALSE, id = rand_id("ancom"))
  standardGeneric("step_ancom"))
methods::setMethod(
  f = "step_ancom",
  signature = c(rec = "recipe"),
  definition = function(rec,
                        out_cut,
                        zero_cut,
                        lib_cut,
                        neg_lb,
                        id) {

  add_step(
      rec,
      step_ancom_new(
        out_cut = out_cut,
        zero_cut = zero_cut,
        lib_cut = lib_cut,
        neg_lb = neg_lb,
        id = id
      )
    )

  }
)

step_ancom_new <- function(out_cut, zero_cut, lib_cut, neg_lb, id) {
  step(
    subclass = "ancom",
    out_cut = out_cut,
    zero_cut = zero_cut,
    lib_cut = lib_cut,
    neg_lb = neg_lb,
    id = id
  )
}
