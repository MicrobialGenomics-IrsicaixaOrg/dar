# SAFE RECIPE SERIALIZATION -------------------------------------------------

recipe_schema_format <- "dar-recipe"
recipe_schema_version <- 1L

#' @noRd
abort_recipe_schema <- function(message, class, record = NULL, field = NULL,
                                parent = NULL) {
  details <- character()
  if (!is.null(record)) {
    details <- c(details, paste0("Record: `", record, "`."))
  }
  if (!is.null(field)) {
    details <- c(details, paste0("Field: `", field, "`."))
  }
  cli::cli_abort(
    c("x" = message, stats::setNames(details, rep("i", length(details)))),
    class = c(class, "dar_error_invalid_recipe_schema"),
    record = record,
    field = field,
    parent = parent
  )
}

#' @noRd
validate_schema_fields <- function(x, allowed, required = allowed,
                                   record = NULL, field = NULL) {
  if (!is.list(x) || (length(x) > 0L && is.null(names(x)))) {
    abort_recipe_schema(
      "Expected a JSON object.", "dar_error_invalid_recipe_record",
      record, field
    )
  }
  object_names <- names(x) %||% character()
  if (anyNA(object_names) || any(!nzchar(object_names)) ||
      anyDuplicated(object_names)) {
    abort_recipe_schema(
      "JSON object fields must be non-empty and unique.",
      "dar_error_invalid_recipe_record", record, field
    )
  }
  unknown <- setdiff(object_names, allowed)
  missing <- setdiff(required, object_names)
  if (length(unknown) > 0L) {
    abort_recipe_schema(
      paste0("Unknown field(s): ", paste(unknown, collapse = ", "), "."),
      "dar_error_invalid_recipe_record", record, field
    )
  }
  if (length(missing) > 0L) {
    abort_recipe_schema(
      paste0("Missing field(s): ", paste(missing, collapse = ", "), "."),
      "dar_error_invalid_recipe_record", record, field
    )
  }
  invisible(x)
}

#' @noRd
validate_unique_json_fields <- function(x, path = "root") {
  if (!is.list(x)) {
    return(invisible(NULL))
  }
  object_names <- names(x)
  if (!is.null(object_names) && anyDuplicated(object_names)) {
    abort_recipe_schema(
      "Duplicate JSON fields are not allowed.",
      "dar_error_invalid_recipe_record", field = path
    )
  }
  for (i in seq_along(x)) {
    child <- if (is.null(object_names)) {
      paste0(path, "[[", i, "]]")
    } else {
      paste0(path, "$", object_names[[i]])
    }
    validate_unique_json_fields(x[[i]], child)
  }
  invisible(NULL)
}

#' @noRd
is_json_scalar <- function(x, type) {
  is.atomic(x) && typeof(x) == type && length(x) == 1L && !is.na(x)
}

#' @noRd
safe_atomic_values <- function(x) {
  missing <- is.na(x)
  safe <- x
  if (is.character(safe)) safe[missing] <- ""
  if (is.logical(safe)) safe[missing] <- FALSE
  if (is.integer(safe)) safe[missing] <- 0L
  if (is.double(safe)) safe[missing] <- 0
  list(
    values = unname(as.list(safe)),
    missing = unname(as.list(missing)),
    names = if (is.null(names(x))) NULL else unname(as.list(names(x)))
  )
}

#' @noRd
encode_recipe_value <- function(value, context = "value") {
  if (is.null(value)) {
    return(list(type = "null"))
  }
  if (is.function(value)) {
    cli::cli_abort(
      "Cannot export {.arg {context}} because executable functions are not supported.",
      class = "dar_error_unserializable_step"
    )
  }
  if (inherits(value, "formula")) {
    validate_safe_formula(value, context)
    return(list(
      type = "formula",
      text = paste(deparse(value), collapse = " ") |> stringr::str_squish()
    ))
  }
  if (is.factor(value)) {
    encoded <- safe_atomic_values(as.character(value))
    return(c(
      list(type = "factor"), encoded,
      list(
        levels = unname(as.list(levels(value))),
        ordered = is.ordered(value)
      )
    ))
  }
  if (is.atomic(value) && typeof(value) %in% c(
    "character", "logical", "integer", "double"
  )) {
    if (is.double(value) && any(!is.finite(value[!is.na(value)]))) {
      cli::cli_abort(
        "Cannot export {.arg {context}} because non-finite numbers are not supported.",
        class = "dar_error_unserializable_step"
      )
    }
    encoded <- safe_atomic_values(value)
    return(c(list(type = typeof(value)), encoded))
  }
  if (is.list(value) && !is.data.frame(value) && !methods::is(value, "S4")) {
    if (length(value) == 0L) {
      return(list(type = "list", entries = list()))
    }
    value_names <- names(value)
    if (is.null(value_names) || anyNA(value_names) ||
        any(!nzchar(value_names)) || anyDuplicated(value_names)) {
      cli::cli_abort(
        "Cannot export {.arg {context}}: nested lists must have non-empty, unique names.",
        class = "dar_error_unserializable_step"
      )
    }
    entries <- purrr::map2(value, value_names, function(item, name) {
      list(name = name, value = encode_recipe_value(item, paste0(context, "$", name)))
    })
    return(list(type = "list", entries = unname(entries)))
  }
  cli::cli_abort(
    "Cannot export {.arg {context}} because its type is not supported.",
    class = "dar_error_unserializable_step"
  )
}

#' @noRd
decode_atomic_value <- function(value, type, context) {
  validate_schema_fields(
    value, c("type", "values", "missing", "names"),
    record = context
  )
  if (!is.list(value$values) || !is.list(value$missing) ||
      length(value$values) != length(value$missing)) {
    abort_recipe_schema(
      "Atomic values and their missingness mask must be arrays of equal length.",
      "dar_error_invalid_recipe_record", context
    )
  }
  missing <- vapply(value$missing, function(x) {
    if (!is_json_scalar(x, "logical")) {
      abort_recipe_schema(
        "The missingness mask must contain only Boolean values.",
        "dar_error_invalid_recipe_record", context, "missing"
      )
    }
    x
  }, logical(1))
  scalar <- switch(
    type,
    character = function(x) {
      if (!is_json_scalar(x, "character")) stop("expected a string")
      x
    },
    logical = function(x) {
      if (!is_json_scalar(x, "logical")) stop("expected a Boolean")
      x
    },
    integer = function(x) {
      if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
          !is.finite(x) || x != as.integer(x)) stop("expected an integer")
      as.integer(x)
    },
    double = function(x) {
      if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
        stop("expected a finite number")
      }
      as.double(x)
    }
  )
  decoded <- tryCatch(
    vapply(value$values, scalar, switch(
      type, character = character(1), logical = logical(1),
      integer = integer(1), double = double(1)
    )),
    error = function(cnd) {
      abort_recipe_schema(
        conditionMessage(cnd), "dar_error_invalid_recipe_record",
        context, "values", parent = cnd
      )
    }
  )
  decoded[missing] <- switch(
    type, character = NA_character_, logical = NA, integer = NA_integer_,
    double = NA_real_
  )
  if (!is.null(value$names)) {
    if (!is.list(value$names) || length(value$names) != length(decoded)) {
      abort_recipe_schema(
        "Names must be null or an array matching the value length.",
        "dar_error_invalid_recipe_record", context, "names"
      )
    }
    decoded_names <- vapply(value$names, function(x) {
      if (!is_json_scalar(x, "character")) {
        abort_recipe_schema(
          "Value names must be strings.", "dar_error_invalid_recipe_record",
          context, "names"
        )
      }
      x
    }, character(1))
    if (any(!nzchar(decoded_names)) || anyDuplicated(decoded_names)) {
      abort_recipe_schema(
        "Value names must be non-empty and unique.",
        "dar_error_invalid_recipe_record", context, "names"
      )
    }
    names(decoded) <- decoded_names
  }
  decoded
}

#' @noRd
decode_recipe_value <- function(value, context = "value") {
  if (!is.list(value) || !is_json_scalar(value$type, "character")) {
    abort_recipe_schema(
      "Encoded values must be objects with a scalar `type` field.",
      "dar_error_invalid_recipe_record", context
    )
  }
  type <- value$type
  if (identical(type, "null")) {
    validate_schema_fields(value, "type", record = context)
    return(NULL)
  }
  if (type %in% c("character", "logical", "integer", "double")) {
    return(decode_atomic_value(value, type, context))
  }
  if (identical(type, "factor")) {
    validate_schema_fields(
      value,
      c("type", "values", "missing", "names", "levels", "ordered"),
      record = context
    )
    atomic <- value[c("type", "values", "missing", "names")]
    atomic$type <- "character"
    decoded <- decode_atomic_value(atomic, "character", context)
    if (!is.list(value$levels) || !is_json_scalar(value$ordered, "logical")) {
      abort_recipe_schema(
        "Factor levels must be an array and `ordered` must be Boolean.",
        "dar_error_invalid_recipe_record", context
      )
    }
    levels <- vapply(value$levels, function(x) {
      if (!is_json_scalar(x, "character")) {
        abort_recipe_schema(
          "Factor levels must be strings.", "dar_error_invalid_recipe_record",
          context, "levels"
        )
      }
      x
    }, character(1))
    return(factor(decoded, levels = levels, ordered = value$ordered))
  }
  if (identical(type, "formula")) {
    validate_schema_fields(value, c("type", "text"), record = context)
    if (!is_json_scalar(value$text, "character")) {
      abort_recipe_schema(
        "Formula text must be one string.", "dar_error_invalid_recipe_record",
        context, "text"
      )
    }
    return(decode_safe_formula(value$text, context))
  }
  if (identical(type, "list")) {
    validate_schema_fields(value, c("type", "entries"), record = context)
    if (!is.list(value$entries)) {
      abort_recipe_schema(
        "List entries must be an array.", "dar_error_invalid_recipe_record",
        context, "entries"
      )
    }
    entries <- purrr::map(seq_along(value$entries), function(index) {
      entry <- value$entries[[index]]
      entry_context <- paste0(context, "$entries[[", index, "]]")
      validate_schema_fields(entry, c("name", "value"), record = entry_context)
      if (!is_json_scalar(entry$name, "character") || !nzchar(entry$name)) {
        abort_recipe_schema(
          "List entry names must be non-empty strings.",
          "dar_error_invalid_recipe_record", entry_context, "name"
        )
      }
      list(
        name = entry$name,
        value = decode_recipe_value(entry$value, paste0(context, "$", entry$name))
      )
    })
    entry_names <- vapply(entries, `[[`, character(1), "name")
    if (anyDuplicated(entry_names)) {
      abort_recipe_schema(
        "Nested list names must be unique.",
        "dar_error_invalid_recipe_record", context, "entries"
      )
    }
    return(stats::setNames(lapply(entries, `[[`, "value"), entry_names))
  }
  abort_recipe_schema(
    paste0("Unsupported encoded value type: `", type, "`."),
    "dar_error_unsafe_recipe_value", context, "type"
  )
}

#' @noRd
safe_formula_operators <- c("~", "+", "-", "*", ":", "^", "|", "(")

#' @noRd
validate_safe_formula_node <- function(node, context) {
  if (is.name(node)) return(invisible(NULL))
  if (is.numeric(node) && length(node) == 1L && is.finite(node)) {
    return(invisible(NULL))
  }
  if (!is.call(node)) {
    abort_recipe_schema(
      "Formula contains an unsupported literal or expression.",
      "dar_error_unsafe_recipe_value", context, "formula"
    )
  }
  operator <- as.character(node[[1L]])
  if (length(operator) != 1L || !operator %in% safe_formula_operators) {
    abort_recipe_schema(
      paste0("Formula call `", paste(operator, collapse = ""), "` is not allowed."),
      "dar_error_unsafe_recipe_value", context, "formula"
    )
  }
  lapply(as.list(node)[-1L], validate_safe_formula_node, context = context)
  invisible(NULL)
}

#' @noRd
validate_safe_formula <- function(formula, context = "formula") {
  if (!inherits(formula, "formula")) {
    abort_recipe_schema(
      "Expected an R formula.", "dar_error_unsafe_recipe_value", context
    )
  }
  validate_safe_formula_node(formula, context)
  invisible(formula)
}

#' @noRd
decode_safe_formula <- function(text, context = "formula") {
  expression <- tryCatch(
    str2lang(text),
    error = function(cnd) {
      abort_recipe_schema(
        "Formula text could not be parsed.", "dar_error_unsafe_recipe_value",
        context, "formula", parent = cnd
      )
    }
  )
  formula <- tryCatch(
    stats::as.formula(expression, env = baseenv()),
    error = function(cnd) {
      abort_recipe_schema(
        "Formula text does not describe a formula.",
        "dar_error_unsafe_recipe_value", context, "formula", parent = cnd
      )
    }
  )
  validate_safe_formula(formula, context)
  formula
}

safe_predicate_calls <- c(
  "(", "+", "-", "*", "/", "^", "%%", "%/%", ">", ">=", "<", "<=",
  "==", "!=", "&", "|", "!", "sum", "length", "mean", "median",
  "min", "max", "var", "sd", "any", "all", "is.na"
)
safe_predicate_aggregates <- c(
  "sum", "mean", "median", "min", "max", "var", "sd", "any", "all"
)

#' @noRd
encode_predicate_node <- function(node, context = ".f") {
  if (is.name(node)) {
    name <- as.character(node)
    if (!identical(name, "x")) {
      cli::cli_abort(
        "Cannot export {.arg {context}} because it captures or references {.val {name}}.",
        class = "dar_error_unserializable_step"
      )
    }
    return(list(node = "symbol", name = "x"))
  }
  if (is.logical(node) && length(node) == 1L && !is.na(node)) {
    return(list(node = "logical", value = node))
  }
  if (is.numeric(node) && length(node) == 1L && is.finite(node)) {
    return(list(node = "number", value = as.double(node)))
  }
  if (!is.call(node)) {
    cli::cli_abort(
      "Cannot export {.arg {context}} because its predicate contains an unsupported expression.",
      class = "dar_error_unserializable_step"
    )
  }
  operator <- as.character(node[[1L]])
  if (length(operator) != 1L || !operator %in% safe_predicate_calls) {
    cli::cli_abort(
      "Cannot export {.arg {context}} because call {.fun {operator}} is not allowed.",
      class = "dar_error_unserializable_step"
    )
  }
  call_args <- as.list(node)[-1L]
  arg_names <- names(call_args) %||% rep("", length(call_args))
  if (any(nzchar(arg_names) &
          !(operator %in% safe_predicate_aggregates & arg_names == "na.rm"))) {
    cli::cli_abort(
      "Cannot export {.arg {context}} because the predicate contains an unsupported named argument.",
      class = "dar_error_unserializable_step"
    )
  }
  arguments <- Map(function(argument, name) {
    list(
      name = if (nzchar(name)) name else NULL,
      value = encode_predicate_node(argument, context)
    )
  }, call_args, arg_names)
  list(node = "call", operator = operator, arguments = unname(arguments))
}

#' @noRd
normalize_filter_predicate <- function(predicate, context = ".f") {
  if (is.character(predicate) && length(predicate) == 1L && !is.na(predicate)) {
    predicate <- tryCatch(
      str2lang(predicate),
      error = function(cnd) {
        cli::cli_abort(
          "{.arg {context}} must be a safely serializable one-expression function.",
          class = "dar_error_invalid_filter_predicate", parent = cnd
        )
      }
    )
    if (!is.call(predicate) || !identical(predicate[[1L]], as.name("function"))) {
      cli::cli_abort(
        "{.arg {context}} must be a function or a string containing one function expression.",
        class = "dar_error_invalid_filter_predicate"
      )
    }
    predicate <- rlang::new_function(
      predicate[[2L]], predicate[[3L]], env = baseenv()
    )
  }
  if (!is.function(predicate)) {
    cli::cli_abort(
      "{.arg {context}} must be a function or a string containing one function expression.",
      class = "dar_error_invalid_filter_predicate"
    )
  }
  predicate_formals <- formals(predicate)
  if (!identical(names(predicate_formals), "x") ||
      !identical(predicate_formals[[1L]], quote(expr = ))) {
    cli::cli_abort(
      "{.arg {context}} must have exactly one argument named {.code x} without a default.",
      class = "dar_error_invalid_filter_predicate"
    )
  }
  body <- body(predicate)
  tryCatch(
    encode_predicate_node(body, context),
    dar_error_unserializable_step = function(cnd) {
      cli::cli_abort(
        "The filter predicate is not safely serializable.",
        class = "dar_error_invalid_filter_predicate",
        parent = cnd
      )
    }
  )
  rlang::new_function(formals(predicate), body, env = baseenv())
}

#' @noRd
encode_filter_predicate <- function(predicate, context = ".f") {
  predicate <- normalize_filter_predicate(predicate, context)
  list(type = "predicate", argument = "x", body = encode_predicate_node(body(predicate), context))
}

#' @noRd
decode_predicate_node <- function(node, context = ".f") {
  if (!is.list(node) || !is_json_scalar(node$node, "character")) {
    abort_recipe_schema(
      "Predicate nodes must be typed JSON objects.",
      "dar_error_unsafe_recipe_value", context, "predicate"
    )
  }
  if (identical(node$node, "symbol")) {
    validate_schema_fields(node, c("node", "name"), record = context)
    if (!is_json_scalar(node$name, "character") || !identical(node$name, "x")) {
      abort_recipe_schema(
        "Predicate symbols are restricted to `x`.",
        "dar_error_unsafe_recipe_value", context, "predicate"
      )
    }
    return(as.name("x"))
  }
  if (identical(node$node, "logical")) {
    validate_schema_fields(node, c("node", "value"), record = context)
    if (!is_json_scalar(node$value, "logical")) {
      abort_recipe_schema(
        "Predicate logical literals must be finite Boolean values.",
        "dar_error_unsafe_recipe_value", context, "predicate"
      )
    }
    return(node$value)
  }
  if (identical(node$node, "number")) {
    validate_schema_fields(node, c("node", "value"), record = context)
    if (!is.numeric(node$value) || length(node$value) != 1L ||
        is.na(node$value) || !is.finite(node$value)) {
      abort_recipe_schema(
        "Predicate numeric literals must be finite numbers.",
        "dar_error_unsafe_recipe_value", context, "predicate"
      )
    }
    return(as.double(node$value))
  }
  if (!identical(node$node, "call")) {
    abort_recipe_schema(
      paste0("Unknown predicate node: `", node$node, "`."),
      "dar_error_unsafe_recipe_value", context, "predicate"
    )
  }
  validate_schema_fields(
    node, c("node", "operator", "arguments"), record = context
  )
  if (!is_json_scalar(node$operator, "character") ||
      !node$operator %in% safe_predicate_calls || !is.list(node$arguments)) {
    abort_recipe_schema(
      "Predicate call is not allowlisted.", "dar_error_unsafe_recipe_value",
      context, "predicate"
    )
  }
  args <- purrr::map(seq_along(node$arguments), function(index) {
    argument <- node$arguments[[index]]
    arg_context <- paste0(context, "$arguments[[", index, "]]")
    validate_schema_fields(argument, c("name", "value"), record = arg_context)
    name <- argument$name
    if (!is.null(name) && (!is_json_scalar(name, "character") ||
        !identical(name, "na.rm") ||
        !node$operator %in% safe_predicate_aggregates)) {
      abort_recipe_schema(
        "Predicate named argument is not allowlisted.",
        "dar_error_unsafe_recipe_value", arg_context, "name"
      )
    }
    list(name = name, value = decode_predicate_node(argument$value, arg_context))
  })
  values <- lapply(args, `[[`, "value")
  arg_names <- vapply(args, function(x) x$name %||% "", character(1))
  names(values) <- arg_names
  as.call(c(list(as.name(node$operator)), values))
}

#' @noRd
decode_filter_predicate <- function(value, context = ".f") {
  validate_schema_fields(
    value, c("type", "argument", "body"), record = context
  )
  if (!is_json_scalar(value$type, "character") ||
      !identical(value$type, "predicate") ||
      !is_json_scalar(value$argument, "character") ||
      !identical(value$argument, "x")) {
    abort_recipe_schema(
      "Filter predicates must declare the single argument `x`.",
      "dar_error_unsafe_recipe_value", context, "predicate"
    )
  }
  body <- decode_predicate_node(value$body, context)
  rlang::new_function(formals(function(x) NULL), body, env = baseenv())
}

#' @noRd
recipe_serialization_registry <- function() {
  registry <- purrr::imap(recipe_step_registry, function(spec, concrete_class) {
    constructor <- if (identical(spec$role, "bake")) {
      "bake"
    } else {
      paste0("step_", spec$method)
    }
    list(
      constructor = constructor,
      role = spec$role,
      concrete_class = concrete_class
    )
  })
  stats::setNames(registry, vapply(registry, `[[`, character(1), "constructor"))
}

#' @noRd
serialize_step_record <- function(configured_step) {
  descriptor <- step_descriptor(configured_step)
  constructor <- if (identical(descriptor$role, "bake")) {
    "bake"
  } else {
    paste0("step_", descriptor$method)
  }
  argument_names <- setdiff(names(configured_step), "id")
  arguments <- stats::setNames(lapply(argument_names, function(name) {
    value <- configured_step[[name]]
    if (identical(constructor, "step_filter_taxa") && identical(name, ".f")) {
      encode_filter_predicate(value, paste0(configured_step$id, "$.f"))
    } else {
      encode_recipe_value(value, paste0(configured_step$id, "$", name))
    }
  }), argument_names)
  list(constructor = constructor, id = configured_step$id, arguments = arguments)
}

#' @noRd
serialize_model_record <- function(model) {
  if (is.null(model)) return(NULL)
  stats::setNames(lapply(names(model), function(name) {
    encode_recipe_value(model[[name]], paste0("model$", name))
  }), names(model))
}

#' @noRd
recipe_serialization_document <- function(rec) {
  list(
    format = recipe_schema_format,
    version = recipe_schema_version,
    model = serialize_model_record(get_model(rec)),
    steps = unname(lapply(rec@steps, serialize_step_record)),
    bakes = if (methods::is(rec, "PrepRecipe")) {
      unname(lapply(rec@bakes, serialize_step_record))
    } else {
      list()
    }
  )
}

#' Export a recipe analysis configuration as versioned JSON
#'
#' Writes the centralized model, ordered preprocessing and DA steps, and any
#' bake configurations using the non-executable `dar-recipe` JSON schema.
#'
#' @param rec A [Recipe()] or [PrepRecipe()] object.
#' @param file_name Output path.
#'
#' @return `file_name`, invisibly.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' rec <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
#'   step_filter_by_prevalence(threshold = 0.1)
#' path <- tempfile(fileext = ".json")
#' export_steps(rec, path)
export_steps <- function(rec, file_name) {
  check_any_recipe(rec)
  if (!is.character(file_name) || length(file_name) != 1L ||
      is.na(file_name) || !nzchar(file_name)) {
    cli::cli_abort(
      "{.arg file_name} must be one non-empty path.",
      class = "dar_error_invalid_recipe_file"
    )
  }
  document <- recipe_serialization_document(rec)
  jsonlite::write_json(
    document, path = file_name, auto_unbox = TRUE, pretty = TRUE,
    null = "null", na = "null", digits = NA
  )
  invisible(file_name)
}

#' @noRd
looks_like_legacy_recipe_file <- function(text) {
  stringr::str_detect(text, "\\[c\\(") ||
  stringr::str_detect(text, "\\}\\s*\\{") ||
    stringr::str_detect(text, '"id"\\s*:\\s*\\[')
}

#' @noRd
read_recipe_document <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) ||
      !nzchar(file) || !file.exists(file)) {
    abort_recipe_schema(
      "{.arg file} must identify one existing readable file.",
      "dar_error_invalid_recipe_file", field = "file"
    )
  }
  text <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (!jsonlite::validate(text)) {
    if (looks_like_legacy_recipe_file(text)) {
      cli::cli_abort(
        c(
          "x" = "This file uses the executable pre-v1 dar recipe format and cannot be imported safely.",
          "i" = "Only migrate files you trust in an isolated environment with the older dar release, then export the resulting recipe with the current version."
        ),
        class = c("dar_error_legacy_recipe_schema", "dar_error_invalid_recipe_schema")
      )
    }
    abort_recipe_schema(
      "The recipe file is not valid JSON.", "dar_error_invalid_recipe_file",
      field = "file"
    )
  }
  document <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = FALSE),
    error = function(cnd) {
      abort_recipe_schema(
        "The recipe JSON could not be decoded.",
        "dar_error_invalid_recipe_file", field = "file", parent = cnd
      )
    }
  )
  validate_unique_json_fields(document)
  document
}

#' @noRd
constructor_arguments <- function(constructor) {
  fn <- get(constructor, envir = asNamespace("dar"))
  args <- formals(fn)
  args[setdiff(names(args), c("rec", "id", "..."))]
}

#' @noRd
decode_step_record <- function(record, section, index) {
  context <- paste0(section, "[[", index, "]]")
  validate_schema_fields(
    record, c("constructor", "id", "arguments"), record = context
  )
  if (!is_json_scalar(record$constructor, "character")) {
    abort_recipe_schema(
      "Constructor must be one string.", "dar_error_invalid_recipe_record",
      context, "constructor"
    )
  }
  registry <- recipe_serialization_registry()
  spec <- registry[[record$constructor]]
  if (is.null(spec)) {
    abort_recipe_schema(
      paste0("Constructor `", record$constructor, "` is not supported."),
      "dar_error_invalid_recipe_record", context, "constructor"
    )
  }
  expected_role <- if (identical(section, "bakes")) "bake" else c("preprocessing", "da")
  if (!spec$role %in% expected_role) {
    abort_recipe_schema(
      paste0("Constructor `", record$constructor, "` is not valid in `", section, "`."),
      "dar_error_invalid_recipe_record", context, "constructor"
    )
  }
  if (!is_json_scalar(record$id, "character") || !nzchar(record$id)) {
    abort_recipe_schema(
      "Step ID must be one non-empty string.",
      "dar_error_invalid_recipe_record", context, "id"
    )
  }
  if (!is.list(record$arguments) ||
      (length(record$arguments) > 0L && is.null(names(record$arguments)))) {
    abort_recipe_schema(
      "Step arguments must be a JSON object.",
      "dar_error_invalid_recipe_record", context, "arguments"
    )
  }
  formals <- constructor_arguments(record$constructor)
  allowed <- names(formals)
  required <- allowed[vapply(formals, identical, logical(1), quote(expr = ))]
  validate_schema_fields(
    record$arguments, allowed, required,
    record = context, field = "arguments"
  )
  arguments <- stats::setNames(lapply(names(record$arguments), function(name) {
    value_context <- paste0(context, "$arguments$", name)
    if (identical(record$constructor, "step_filter_taxa") && identical(name, ".f")) {
      decode_filter_predicate(record$arguments[[name]], value_context)
    } else {
      decode_recipe_value(record$arguments[[name]], value_context)
    }
  }), names(record$arguments))
  list(
    constructor = record$constructor,
    id = record$id,
    arguments = arguments,
    role = spec$role
  )
}

#' @noRd
decode_model_record <- function(model) {
  if (is.null(model)) return(NULL)
  fields <- c(
    "formula", "targets", "tax_level", "time", "reference", "contrasts",
    "na_action", "max_time_levels"
  )
  validate_schema_fields(model, fields, record = "model")
  decoded <- stats::setNames(lapply(fields, function(name) {
    decode_recipe_value(model[[name]], paste0("model$", name))
  }), fields)
  if (!inherits(decoded$formula, "formula")) {
    abort_recipe_schema(
      "The central model formula must use the formula value type.",
      "dar_error_invalid_recipe_record", "model", "formula"
    )
  }
  decoded
}

#' @noRd
decode_recipe_document <- function(document) {
  validate_schema_fields(
    document, c("format", "version", "model", "steps", "bakes"),
    record = "root"
  )
  if (!is_json_scalar(document$format, "character") ||
      !identical(document$format, recipe_schema_format)) {
    abort_recipe_schema(
      "Unsupported recipe serialization format.",
      "dar_error_unsupported_recipe_schema", "root", "format"
    )
  }
  if (!is.numeric(document$version) || length(document$version) != 1L ||
      is.na(document$version) || document$version != recipe_schema_version) {
    abort_recipe_schema(
      paste0("Unsupported recipe schema version; expected ", recipe_schema_version, "."),
      "dar_error_unsupported_recipe_schema", "root", "version"
    )
  }
  if (!is.list(document$steps) || !is.list(document$bakes)) {
    abort_recipe_schema(
      "`steps` and `bakes` must be JSON arrays.",
      "dar_error_invalid_recipe_record", "root"
    )
  }
  steps <- lapply(seq_along(document$steps), function(index) {
    decode_step_record(document$steps[[index]], "steps", index)
  })
  bakes <- lapply(seq_along(document$bakes), function(index) {
    decode_step_record(document$bakes[[index]], "bakes", index)
  })
  ids <- vapply(c(steps, bakes), `[[`, character(1), "id")
  duplicate_ids <- unique(ids[duplicated(ids)])
  if (length(duplicate_ids) > 0L) {
    abort_recipe_schema(
      paste0("Step IDs must be unique; duplicated: ", paste(duplicate_ids, collapse = ", "), "."),
      "dar_error_invalid_recipe_record", "root", "id"
    )
  }
  list(model = decode_model_record(document$model), steps = steps, bakes = bakes)
}

#' @noRd
apply_serialized_step <- function(rec, spec) {
  constructor <- get(spec$constructor, envir = asNamespace("dar"))
  do.call(
    constructor,
    c(list(rec = rec), spec$arguments, list(id = spec$id)),
    envir = parent.frame()
  )
}

#' @noRd
assert_prepared_import_matches <- function(rec, document) {
  current <- recipe_serialization_document(rec)
  canonical_json <- function(x) {
    as.character(jsonlite::toJSON(
      x, auto_unbox = TRUE, null = "null", na = "null", digits = NA
    ))
  }
  model_matches <- identical(
    canonical_json(current$model), canonical_json(document$model)
  )
  steps_match <- identical(
    canonical_json(current$steps), canonical_json(document$steps)
  )
  if (!model_matches || !steps_match) {
    abort_recipe_schema(
      "The serialized model and steps do not match the supplied `PrepRecipe`.",
      "dar_error_recipe_import_mismatch", "root"
    )
  }
  invisible(NULL)
}

#' Import a recipe configuration from versioned JSON
#'
#' Recipe files are treated as untrusted input. The complete document is
#' decoded and validated before supported constructors are applied. Importing
#' bake configurations into an unprepared recipe never runs analysis unless
#' `prepare = TRUE` is supplied explicitly.
#'
#' @param rec A [Recipe()] or [PrepRecipe()] object.
#' @param file Input JSON path.
#' @param prepare Whether an unprepared recipe may run [prep()] when the file
#'   contains bake configurations. Defaults to `FALSE`.
#' @param parallel Passed to [prep()] only when `prepare = TRUE`.
#' @param workers Passed to [prep()] only when `prepare = TRUE`.
#'
#' @return The configured `Recipe` or `PrepRecipe`.
#' @export
#' @autoglobal
#' @examples
#' data(metaHIV_phy)
#' configured <- recipe(metaHIV_phy) |>
#'   add_model(~ RiskGroup2, targets = "RiskGroup2", tax_level = "Species") |>
#'   step_filter_by_prevalence(threshold = 0.1)
#' path <- tempfile(fileext = ".json")
#' export_steps(configured, path)
#' restored <- import_steps(recipe(metaHIV_phy), path)
import_steps <- function(
  rec,
  file,
  prepare = FALSE,
  parallel = TRUE,
  workers = future::availableCores(constraints = "connections-16")
) {
  check_any_recipe(rec)
  if (!is.logical(prepare) || length(prepare) != 1L || is.na(prepare)) {
    cli::cli_abort(
      "{.arg prepare} must be one non-missing Boolean value.",
      class = "dar_error_invalid_recipe_import"
    )
  }
  document <- read_recipe_document(file)
  decoded <- decode_recipe_document(document)
  is_prepared <- methods::is(rec, "PrepRecipe")
  if (!is_prepared && length(decoded$bakes) > 0L && !prepare) {
    cli::cli_abort(
      c(
        "x" = "The recipe file contains bake configurations, but {.arg prepare} is {.code FALSE}.",
        "i" = "Set {.code prepare = TRUE} explicitly to run analysis, or import a file without bakes."
      ),
      class = "dar_error_import_requires_prepare"
    )
  }
  candidate <- rec
  if (is_prepared) {
    assert_prepared_import_matches(candidate, document)
  } else {
    if (!is.null(decoded$model)) {
      candidate <- do.call(add_model, c(list(rec = candidate), decoded$model))
    }
    for (step_spec in decoded$steps) {
      candidate <- apply_serialized_step(candidate, step_spec)
    }
    if (length(decoded$bakes) > 0L) {
      candidate <- prep(candidate, parallel = parallel, workers = workers)
    }
  }
  for (bake_spec in decoded$bakes) {
    candidate <- apply_serialized_step(candidate, bake_spec)
  }
  validate_recipe_object(candidate)
  candidate
}
