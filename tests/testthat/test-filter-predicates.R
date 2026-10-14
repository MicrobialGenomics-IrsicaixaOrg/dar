make_filter_predicate_phy <- function() {
  counts <- rbind(
    zero = rep(0, 6),
    constant = rep(2, 6),
    variable = 0:5,
    outlier = c(0, 0, 0, 0, 0, 10)
  )
  colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))
  taxonomy <- matrix(
    paste0("Species_", rownames(counts)), ncol = 1,
    dimnames = list(rownames(counts), "Species")
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy)
  )
}

filter_stats_predicates <- list(
  median = function(x) median(x) > 1,
  var = function(x) var(x) > 1,
  sd = function(x) sd(x) > 1
)

test_that("stats predicates filter taxa in prep as functions and strings", {
  phy <- make_filter_predicate_phy()
  expected_ids <- list(
    median = c("constant", "variable"),
    var = c("variable", "outlier"),
    sd = c("variable", "outlier")
  )
  for (name in names(filter_stats_predicates)) {
    predicate <- filter_stats_predicates[[name]]
    expected <- phyloseq::filter_taxa(phy, predicate, prune = TRUE)
    for (as_string in c(FALSE, TRUE)) {
      argument <- if (as_string) {
        paste(deparse(predicate), collapse = " ")
      } else {
        predicate
      }
      rec <- step_filter_taxa(recipe(phy), argument, id = "filter__stats")
      prepared <- prep(rec, parallel = FALSE)
      expect_s4_class(prepared, "PrepRecipe")
      expect_identical(
        phyloseq::taxa_names(get_phy(prepared)), expected_ids[[name]]
      )
      expect_identical(get_phy(prepared), expected)
      expect_identical(get_phy(rec), phy)
    }
  }
})

test_that("stats predicates execute after JSON decoding and recipe import", {
  phy <- make_filter_predicate_phy()
  for (predicate in filter_stats_predicates) {
    for (as_string in c(FALSE, TRUE)) {
      argument <- if (as_string) {
        paste(deparse(predicate), collapse = " ")
      } else {
        predicate
      }
      rec <- step_filter_taxa(recipe(phy), argument, id = "filter__stats")
      path <- tempfile(fileext = ".json")
      export_steps(rec, path)
      document <- jsonlite::fromJSON(path, simplifyVector = FALSE)
      decoded <- dar:::decode_filter_predicate(document$steps[[1]]$arguments$.f)
      expect_identical(decoded(0:5), predicate(0:5))

      original <- recipe(phy)
      imported <- import_steps(original, path)
      expect_identical(
        dar:::recipe_serialization_document(imported),
        dar:::recipe_serialization_document(rec)
      )
      expect_identical(
        get_phy(prep(imported, parallel = FALSE)),
        get_phy(prep(rec, parallel = FALSE))
      )
      expect_length(original@steps, 0L)
      expect_identical(get_phy(original), phy)
    }
  }
})

test_that("stats predicates preserve na.rm through serialization", {
  predicates <- list(
    function(x) median(x, na.rm = TRUE) > 1,
    function(x) var(x, na.rm = TRUE) >= 1,
    function(x) sd(x, na.rm = TRUE) >= 1
  )
  values <- c(1, 2, 3, NA_real_)
  for (predicate in predicates) {
    normalized <- dar:::normalize_filter_predicate(predicate)
    encoded <- dar:::encode_filter_predicate(normalized)
    decoded <- dar:::decode_filter_predicate(encoded)
    expect_identical(normalized(values), predicate(values))
    expect_identical(decoded(values), predicate(values))
    expect_true(decoded(values))
    expect_identical(dar:::encode_filter_predicate(decoded), encoded)
  }
})

test_that("predicate environments expose only fixed trusted stats bindings", {
  source_environment <- new.env(parent = globalenv())
  shadow_called <- FALSE
  shadow <- function(...) {
    shadow_called <<- TRUE
    stop("A captured function must not execute")
  }
  for (name in names(filter_stats_predicates)) {
    source_environment[[name]] <- shadow
  }

  for (predicate in filter_stats_predicates) {
    environment(predicate) <- source_environment
    normalized <- dar:::normalize_filter_predicate(predicate)
    decoded <- dar:::decode_filter_predicate(
      dar:::encode_filter_predicate(normalized)
    )
    for (safe in list(normalized, decoded)) {
      expect_true(safe(0:5))
      env <- environment(safe)
      expect_identical(parent.env(env), baseenv())
      expect_identical(ls(env, all.names = TRUE), c("median", "sd", "var"))
      expect_true(environmentIsLocked(env))
      for (name in names(filter_stats_predicates)) {
        expect_identical(env[[name]], getExportedValue("stats", name))
        expect_true(bindingIsLocked(name, env))
      }
      expect_error(assign("median", shadow, envir = env))
      expect_error(assign("sum", shadow, envir = env))
    }
  }
  expect_false(shadow_called)
})

test_that("safe stats bindings do not expand accepted predicate syntax", {
  rec <- recipe(make_filter_predicate_phy())
  unsafe <- c(
    "function(x) median(x) > cutoff",
    "function(x) stats::median(x) > 1",
    "function(x) { median(x) > 1 }",
    "function(x) cutoff <- median(x)",
    "function(x) helper(x)",
    "function(x) eval(x)",
    "function(x) get(x)",
    "function(x) source(x)",
    "function(x) system(x)",
    "function(x) quantile(x)",
    "function(x) rnorm(x)",
    "function(x) median(x, trim = TRUE) > 1",
    "function(y) median(y) > 1",
    "function(x = 1) median(x) > 1",
    "function(x, y) median(x) > 1",
    "function(x) median(x[1]) > 1"
  )
  for (predicate in unsafe) {
    expect_error(
      step_filter_taxa(rec, predicate),
      class = "dar_error_invalid_filter_predicate"
    )
  }
  expect_error(
    step_filter_taxa(rec, list(function(x) median(x) > 1)),
    class = "dar_error_invalid_filter_predicate"
  )
  expect_length(rec@steps, 0L)
})

test_that("tampered stats predicate imports remain transactional", {
  phy <- make_filter_predicate_phy()
  configured <- recipe(phy) |>
    step_filter_taxa(function(x) median(x) > 1, id = "filter__stats")
  original <- recipe(phy) |>
    step_filter_by_prevalence(threshold = 0.1, id = "prevalence__original")
  before <- dar:::recipe_serialization_document(original)
  path <- tempfile(fileext = ".json")
  export_steps(configured, path)
  valid <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  for (operator in c("system", "eval", "get", "::", "quantile", "rnorm")) {
    document <- valid
    document$steps[[1]]$arguments$.f$body$arguments[[1]]$value$operator <-
      operator
    jsonlite::write_json(
      document, path, auto_unbox = TRUE, null = "null", na = "null", digits = NA
    )
    expect_error(
      import_steps(original, path), class = "dar_error_unsafe_recipe_value"
    )
    expect_identical(dar:::recipe_serialization_document(original), before)
    expect_identical(get_phy(original), phy)
  }
})
