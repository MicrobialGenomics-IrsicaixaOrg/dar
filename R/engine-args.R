# ADVANCED ENGINE ARGUMENTS ---------------------------------------------------

#' @noRd
engine_args_registry <- function() {
  list(
    deseq = list(
      public = c(
        "test", "fitType", "betaPrior", "type", "max_significance",
        "log2FC", "rarefy", "id", "engine_args"
      ),
      stages = list(
        size_factors = list(
          package = "DESeq2", fun = "estimateSizeFactors",
          protected = c("object", "geoMeans")
        ),
        fit = list(
          package = "DESeq2", fun = "DESeq",
          protected = c("object", "full", "reduced", "modelMatrixType")
        ),
        results = list(
          package = "DESeq2", fun = "results",
          protected = c(
            "object", "contrast", "name", "listValues", "format",
            "saveCols", "test", "addMLE", "tidy"
          )
        ),
        shrink = list(
          package = "DESeq2", fun = "lfcShrink",
          protected = c(
            "dds", "coef", "contrast", "res", "type", "returnList",
            "format", "saveCols"
          )
        )
      )
    ),
    aldex = list(
      public = c(
        "max_significance", "mc.samples", "denom", "rarefy", "id",
        "engine_args"
      ),
      stages = list(
        clr = list(
          package = "ALDEx2", fun = "aldex.clr",
          protected = c("reads", "conds")
        ),
        fit = list(
          package = "ALDEx2", fun = "aldex.glm", protected = "clr"
        ),
        effect = list(
          package = "ALDEx2", fun = "aldex.glm.effect", protected = "clr"
        )
      )
    ),
    ancom = list(
      public = c(
        "fix_formula", "rand_formula", "p_adj_method", "prv_cut",
        "lib_cut", "s0_perc", "group", "struc_zero", "neg_lb", "alpha",
        "n_cl", "verbose", "global", "pairwise", "dunnet", "trend",
        "rarefy", "id", "engine_args"
      ),
      stages = list(
        fit = list(
          package = "ANCOMBC", fun = "ancombc2",
          protected = c(
            "data", "taxa_are_rows", "assay.type", "assay_name", "rank",
            "tax_level", "aggregate_data", "meta_data", "fix_formula",
            "rand_formula"
          )
        )
      )
    ),
    corncob = list(
      public = c(
        "phi.formula", "formula_null", "phi.formula_null", "link",
        "phi.link", "test", "boot", "B", "filter_discriminant",
        "fdr_cutoff", "fdr", "log2FC", "rarefy", "id", "engine_args"
      ),
      stages = list(
        fit = list(
          package = "corncob", fun = "bbdml",
          protected = c("formula", "phi.formula", "data", "link", "phi.link")
        )
      )
    ),
    linda = list(
      public = c(
        "prev_filter", "mean_abund_filter", "max_abund_filter",
        "winsorize", "outlier_pct", "adaptive", "zero_handling",
        "pseudo_count", "corr_cut", "p_adj_method", "alpha", "n_cpus",
        "rarefy", "id", "engine_args"
      ),
      stages = list(
        fit = list(
          package = "MicrobiomeStat", fun = "linda",
          protected = c(
            "feature.dat", "meta.dat", "formula", "feature.dat.type",
            "prev.filter", "mean.abund.filter", "max.abund.filter",
            "is.winsor", "outlier.pct", "adaptive", "zero.handling",
            "pseudo.cnt", "corr.cut", "p.adj.method", "alpha", "n.cores"
          )
        )
      )
    ),
    maaslin = list(
      public = c(
        "min_abundance", "min_prevalence", "min_variance", "normalization",
        "transform", "max_significance", "random_effects", "correction",
        "standardize", "reference", "median_comparison_abundance",
        "rarefy", "id", "engine_args"
      ),
      stages = list(
        fit = list(
          package = "maaslin3", fun = "maaslin3",
          protected = c(
            "input_data", "input_metadata", "output", "formula",
            "fixed_effects", "reference", "random_effects", "evaluate_only"
          )
        ),
        contrast = list(
          package = "maaslin3", fun = "maaslin_contrast_test",
          protected = c(
            "maaslin3_fit", "contrast_mat", "rhs", "evaluate_only"
          )
        )
      )
    ),
    wilcox = list(
      public = c(
        "norm_method", "max_significance", "p_adj_method", "rarefy", "id",
        "engine_args"
      ),
      stages = list(
        test = list(
          package = "stats", fun = "wilcox.test", protected = c("x", "y")
        )
      )
    ),
    lefse = list(
      public = c(
        "kruskal.threshold", "wilcox.threshold", "lda.threshold",
        "subclassCol", "assay", "trim.names", "rarefy", "id",
        "engine_args"
      ),
      stages = list(
        transform = list(
          package = "lefser", fun = "relativeAb", protected = "se"
        ),
        fit = list(
          package = "lefser", fun = "lefser",
          protected = c("relab", "classCol")
        )
      )
    )
  )
}

#' @noRd
engine_stage_spec <- function(engine, stage) {
  registry <- engine_args_registry()
  if (!engine %in% names(registry) ||
      !stage %in% names(registry[[engine]]$stages)) {
    cli::cli_abort(
      "Unknown advanced argument stage {.val {stage}} for engine {.val {engine}}.",
      class = "dar_error_invalid_engine_args"
    )
  }
  registry[[engine]]$stages[[stage]]
}

#' @noRd
engine_stage_function <- function(engine, stage) {
  spec <- engine_stage_spec(engine, stage)
  getExportedValue(spec$package, spec$fun)
}

#' Validate and normalize advanced native-engine arguments
#' @noRd
normalize_engine_args <- function(engine, engine_args = list()) {
  registry <- engine_args_registry()
  if (!engine %in% names(registry)) {
    cli::cli_abort(
      "Unknown DA engine {.val {engine}}.",
      class = "dar_error_invalid_engine_args"
    )
  }
  if (!is.list(engine_args)) {
    cli::cli_abort(
      "{.arg engine_args} for {.val {engine}} must be a named list.",
      class = "dar_error_invalid_engine_args"
    )
  }
  if (length(engine_args) == 0L) {
    return(list())
  }
  stage_names <- names(engine_args)
  if (is.null(stage_names) || anyNA(stage_names) || any(!nzchar(stage_names)) ||
      anyDuplicated(stage_names)) {
    cli::cli_abort(
      "{.arg engine_args} stages must have non-empty, unique names.",
      class = "dar_error_invalid_engine_args"
    )
  }
  invalid_stages <- setdiff(stage_names, names(registry[[engine]]$stages))
  if (length(invalid_stages) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Invalid advanced argument stages for {.val {engine}}: {.val {invalid_stages}}.",
        "i" = "Available stages: {.val {names(registry[[engine]]$stages)}}."
      ),
      class = "dar_error_invalid_engine_args"
    )
  }

  normalized <- purrr::imap(engine_args, function(arguments, stage) {
    if (!is.list(arguments)) {
      cli::cli_abort(
        "Stage {.val {stage}} for {.val {engine}} must contain a named list.",
        class = "dar_error_invalid_engine_args"
      )
    }
    if (length(arguments) == 0L) {
      return(NULL)
    }
    argument_names <- names(arguments)
    if (is.null(argument_names) || anyNA(argument_names) ||
        any(!nzchar(argument_names)) || anyDuplicated(argument_names)) {
      cli::cli_abort(
        "Arguments in stage {.val {stage}} for {.val {engine}} must have non-empty, unique names.",
        class = "dar_error_invalid_engine_args"
      )
    }
    spec <- engine_stage_spec(engine, stage)
    protected <- unique(c(registry[[engine]]$public, spec$protected))
    collisions <- intersect(argument_names, protected)
    if (length(collisions) > 0L) {
      cli::cli_abort(
        "Stage {.val {stage}} for {.val {engine}} cannot override dar-managed arguments: {.arg {collisions}}.",
        class = "dar_error_invalid_engine_args"
      )
    }

    native <- engine_stage_function(engine, stage)
    native_formals <- names(formals(native))
    if (!"..." %in% native_formals) {
      invalid <- setdiff(argument_names, native_formals)
      if (length(invalid) > 0L) {
        cli::cli_abort(
          "Unknown native arguments for {.val {engine}} stage {.val {stage}}: {.arg {invalid}}.",
          class = "dar_error_invalid_engine_args"
        )
      }
    }
    arguments
  })
  purrr::compact(normalized)
}

#' @noRd
engine_args_present <- function(step) {
  length(step$engine_args %||% list()) > 0L
}

#' @noRd
validate_engine_args_model <- function(rec) {
  if (!is.null(get_model(rec))) {
    return(invisible(TRUE))
  }
  advanced <- purrr::keep(rec@steps, function(item) {
    is_da_step(item) && engine_args_present(item)
  })
  if (length(advanced) > 0L) {
    ids <- purrr::map_chr(advanced, ~ .x$id)
    cli::cli_abort(
      c(
        "x" = "Advanced native-engine arguments require a centralized model.",
        "i" = "Legacy steps with {.arg engine_args}: {.val {ids}}."
      ),
      class = "dar_error_engine_args_requires_model"
    )
  }
  invisible(TRUE)
}

#' @noRd
check_engine_args_execution <- function(rec, engine, engine_args, step_id) {
  engine_args <- normalize_engine_args(engine, engine_args)
  if (is.null(get_model(rec)) && length(engine_args) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Advanced native-engine arguments require a centralized model.",
        "i" = "Legacy step {.val {step_id}} configured {.arg engine_args}."
      ),
      class = "dar_error_engine_args_requires_model"
    )
  }
  engine_args
}

#' Execute one registered native-engine stage
#' @noRd
exec_engine_stage <- function(engine, stage, engine_args, fixed,
                              defaults = list()) {
  extra <- engine_args[[stage]] %||% list()
  collisions <- intersect(names(extra), names(fixed))
  if (length(collisions) > 0L) {
    cli::cli_abort(
      "Stage {.val {stage}} for {.val {engine}} cannot override fixed arguments: {.arg {collisions}}.",
      class = "dar_error_invalid_engine_args"
    )
  }
  defaults <- defaults[setdiff(names(defaults), names(extra))]
  arguments <- c(fixed, defaults, extra)
  rlang::exec(engine_stage_function(engine, stage), !!!arguments)
}

#' @noRd
check_unused_engine_args <- function(engine, engine_args, used_stages) {
  unused <- setdiff(names(engine_args), used_stages)
  if (length(unused) > 0L) {
    cli::cli_abort(
      "Advanced argument stages configured but not executed for {.val {engine}}: {.val {unused}}.",
      class = "dar_error_invalid_engine_args"
    )
  }
  invisible(TRUE)
}
