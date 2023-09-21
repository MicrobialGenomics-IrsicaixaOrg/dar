#' Loads Phyloseq data
#'
#' This function returns a validated Phyloseq object by loading it directly from
#' a file with the .rds extension. Alternatively, this function can also take
#' three text files as input that will be used to construct and validate the
#' Phyloseq object: - Counts matrix with the otu_id in the first column. -
#' Taxonomic annotation matrix with the otu_id in the first column. - Metadata
#' annotation with sample_id in the first column.
#'
#' @param data_path List of length 1 or 3, with the paths of the input files.
#'
#' @return a phylseq-class object
#' @export
#'
#' @tests testthat
#' 
#' ## read phyloseq from rds file
#' good_file <- system.file("extdata", "metaHIV_phy.rds", package = "dar")
#' bad_file <- system.file("extdata", "invented_file.rds", package = "dar")
#' expect_s4_class(read_data(good_file), "phyloseq")
#' expect_error(read_data(bad_file))
#' 
#' ## read phyloseq from multiple files
#' good_files <- c(
#'   system.file("extdata", "metaHIV_counts.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_taxas.txt", package = "dar")
#' )
#' expect_s4_class(read_data(good_files), "phyloseq")
#' 
#' only_two_files <- c(
#'   system.file("extdata", "metaHIV_counts.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_metadata.txt", package = "dar")
#' )
#' expect_error(read_data(only_otu_phy))
#' 
#' duplicated_files <- c(
#'   system.file("extdata", "metaHIV_counts.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
#' system.file("extdata", "metaHIV_metadata.txt", package = "dar")
#' )
#' expect_error(read_data(duplicated_files))
#' 
#' @examples
#' # From a phyloseq object saved with .rds extension.
#' system.file("extdata", "metaHIV_phy.rds", package = "dar") %>%
#'   read_data()
#'
#' # From the three components of a phyloseq object saved as a plain text.
#' data_path <- c(
#'   system.file("extdata", "metaHIV_counts.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_metadata.txt", package = "dar"),
#'   system.file("extdata", "metaHIV_taxas.txt", package = "dar")
#' )
#'
#' read_data(data_path)
read_data <- function(data_path) {
  if (!length(data_path) %in% c(1, 3)) {
    text <-
      "Input data must be length 1 (phyloseq objects) or length 3 (text files)"
    rlang::abort(c(text))
  }
  
  if (length(data_path) == 1) {
    phy <-
      purrr::pluck(data_path, 1) %>%
      read_phyloseq() %>%
      validate_phyloseq()
  }
  
  if (length(data_path) == 3) {
    files <-
      data_path %>%
      purrr::set_names() %>%
      purrr::map(read_file) %>%
      unique()
    
    if (length(files) < 3) {
      rlang::abort(c(x = "The three input files must be unique"))
    }
    
    counts_df <-
      files %>%
      purrr::map_dfr(~ {
        class <- try(validate_otu(.x), silent = TRUE)
        if (!inherits(class, "try-error")) {
          .x
        }
      })
    
    taxa_df <-
      files %>%
      purrr::discard(~ waldo::compare(.x, counts_df) %>% length() == 0) %>%
      purrr::map_dfr(~ {
        if (all(counts_df[[1]] %in% .x[[1]])) {
          .x
        }
      }) %>%
      validate_tax_table()
    
    metadata_df <-
      files %>%
      purrr::discard(
        ~ waldo::compare(.x, counts_df) %>% length() == 0 |
          waldo::compare(.x, taxa_df) %>% length() == 0
      ) %>%
      purrr::pluck(1) %>%
      validate_sample_data()
    
    phy <-
      phyloseq::phyloseq(
        data.frame(counts_df, row.names = 1) %>%
          phyloseq::otu_table(taxa_are_rows = TRUE),
        data.frame(taxa_df, row.names = 1) %>%
          as.matrix() %>%
          phyloseq::tax_table(),
        data.frame(metadata_df, row.names = 1) %>%
          phyloseq::sample_data()
      ) %>%
      validate_phyloseq()
  }
  phy
}

#' @rdname read_data
#' @keywords internal
#' @tests testthat
#' df <- read_file(system.file("extdata", "metaHIV_counts.txt", package = "dar"))
#' 
#' expect_invisible(validate_otu(df))
#' expect_error(validate_otu(dplyr::select(df, -1)))
#' expect_error(validate_otu(dplyr::mutate(
#'   df, otu_id = ifelse(otu_id == "Otu_1", NA_character_, otu_id)
#' )))
#' expect_error(validate_otu(dplyr::mutate(
#'   df, otu_id = ifelse(otu_id == "Otu_1", "Otu_2", otu_id)
#' )))
#' expect_error(validate_otu(dplyr::mutate(df, dplyr::across(2, as.character))))
#' expect_error(validate_otu(dplyr::mutate(df, dplyr::across(
#'  2, ~ ifelse(.x == 0, NA_integer_, .x)
#' ))))
#' expect_error(validate_otu(dplyr::mutate(df, dplyr::across(2, ~ .x * -1))))
validate_otu <- function(otu) {
  ids <- otu %>% dplyr::pull(1)
  
  check_1 <- NULL
  if (!is.character(ids)) {
    text <- stringr::str_c(
      "The first column must correspond to the otu identifiers. Therefore it ",
      "must be type character.",
      collapse = ""
    )
    check_1 <- c("x" = text)
  }
  
  check_2 <- NULL
  if (!all(not_na(ids))) {
    check_2 <- c("x" = "Otu identifiers cannot be NA.")
  }
  
  check_3 <- NULL
  if (length(unique(ids)) != length(ids)) {
    check_3 <- c("x" = "Otu identifiers bust be unique.")
  }
  
  check_4 <- NULL
  all_numeric <-
    dplyr::select(otu, -1) %>%
    purrr::map_lgl(is.numeric) %>%
    all()
  
  if (!all_numeric) {
    check_4 <-
      c("x" = "Apart from the first column, the others must be numeric.")
  }
  
  check_5 <- NULL
  no_nas <-
    dplyr::select(otu, -1) %>%
    purrr::map_lgl( ~ {
      all(not_na(.x))
    }) %>%
    all()
  
  if (!no_nas) {
    check_5 <- c("x" = "Count values cannot be NA.")
  }
  
  check_6 <- NULL
  all_positive <-
    dplyr::select(otu, -1) %>%
    purrr::map_lgl(~ {
      all(stats::na.omit(.x) >= 0)
    }) %>%
    all()
  
  if (!all_positive) {
    check_6 <- c("x" = "Count values must be positive.")
  }
  
  check <- c(check_1, check_2, check_3, check_4, check_5, check_6)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while validating otu counts table.", check))
  }
  
  invisible(TRUE)
}

#' @rdname read_data
#' @keywords internal
#' @tests testthat
#' df <- read_file(system.file("extdata", "metaHIV_metadata.txt", package = "dar"))
#' 
#' expect_s3_class(validate_sample_data(df), c("tbl_df", "tbl", "data.frame"))
#' expect_equal(dim(df), c(156, 4))
#' expect_error(validate_sample_data(dplyr::select(df, -1)))
#' expect_error(validate_sample_data(dplyr::mutate(
#'   df, sample_id = ifelse(sample_id == "Sample_169", NA_character_, sample_id)
#' )))
#' expect_error(validate_sample_data(dplyr::mutate(
#'   df, sample_id = ifelse(sample_id == "Sample_169", "Sample_162", sample_id)
#' )))
#' expect_error(validate_sample_data(dplyr::mutate(df, sample_id = 1:nrow(df))))
validate_sample_data <- function(sample_data) {
  ids <- sample_data %>% dplyr::pull(1)
  
  check_1 <- NULL
  if (!is.character(ids)) {
    text <- stringr::str_c(
      "The first column must correspond to the samples identifiers. Therefore ",
      "it must be type character.",
      collapse = ""
    )
    check_1 <- c("x" = text)
  }
  
  check_2 <- NULL
  if (!all(not_na(ids))) {
    check_2 <- c("x" = "Samples identifiers cannot be NA.")
  }
  
  check_3 <- NULL
  if (length(unique(ids)) != length(ids)) {
    check_3 <- c("x" = "Samples identifiers must be unique.")
  }
  
  check <- c(check_1, check_2, check_3)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while validating sample data.", check))
  }
  
  sample_data
}

#' @rdname read_data
#' @keywords internal
#' @tests testthat
#' df <- read_file(system.file("extdata", "metaHIV_taxas.txt", package = "dar"))
#' 
#' expect_s3_class(validate_tax_table(df), c("tbl_df", "tbl", "data.frame"))
#' expect_equal(dim(df), c(451, 8))
#' expect_error(validate_tax_table(dplyr::select(df, -1)))
#' expect_error(validate_tax_table(dplyr::mutate(
#'   df, otu_id = ifelse(otu_id == "Otu_1", NA_character_, otu_id)
#' )))
#' expect_error(validate_tax_table(dplyr::mutate(
#'   df, otu_id = ifelse(otu_id == "Otu_1", "Otu_2", otu_id)
#' )))
#' expect_error(validate_tax_table(dplyr::mutate(df, dplyr::across(2, as.factor))))
#' expect_error(validate_tax_table(dplyr::mutate(df, otu_id = 1:nrow(df))))
validate_tax_table <- function(tax_table) {
  ids <- tax_table %>% dplyr::pull(1)
  
  check_1 <- NULL
  if (!is.character(ids)) {
    text <- stringr::str_c(
      "The first column must correspond to the samples identifiers. Therefore ",
      "it must be type character.",
      collapse = ""
    )
    check_1 <- c("x" = text)
  }
  
  check_2 <- NULL
  if (!all(not_na(ids))) {
    check_2 <- c("x" = "Otu identifiers cannot be NA.")
  }
  
  check_3 <- NULL
  if (length(unique(ids)) != length(ids)) {
    check_3 <- c("x" = "Otu identifiers bust be unique.")
  }
  
  check_4 <- NULL
  all_char <-
    dplyr::select(tax_table, -1) %>% purrr::map_lgl(is.character) %>% all()
  if (!all_char) {
    check_4 <- c("x" = "All annotation columns must be character.")
  }
  
  check <- c(check_1, check_2, check_3, check_4)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while validating taxonomic anntation table.", check))
  }
  tax_table
}

#' @rdname read_data
#' @keywords internal
#' @tests testthat
#' ## phyloseq validation
#' phy <-
#'   system.file("extdata", "metaHIV_phy.rds", package = "dar") %>%
#'   read_phyloseq()
#' 
#' no_tax_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::sample_data(phy))
#' no_sam_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy), phyloseq::tax_table(phy))
#' only_otu_phy <- phyloseq::phyloseq(phyloseq::otu_table(phy))
#' 
#' expect_error(validate_phyloseq(
#'   read_file(system.file("extdata", "metaHIV_taxas.txt", package = "dar")), 
#'   slots = ""
#' ))
#' expect_s4_class(validate_phyloseq(phy), "phyloseq")
#' expect_s4_class(validate_phyloseq(no_sam_phy, "tax_table"), "phyloseq")
#' expect_s4_class(validate_phyloseq(no_tax_phy, "sample_data"), "phyloseq")
#' expect_s4_class(validate_phyloseq(only_otu_phy, ""), "otu_table")
#' expect_error(validate_phyloseq(no_tax_phy))
#' expect_error(validate_phyloseq(no_sam_phy))
#' expect_error(validate_phyloseq(only_otu_phy))
#' expect_error(validate_phyloseq(no_tax_phy, "tax_table"))
#' expect_error(validate_phyloseq(no_sam_phy, "sample_data"))
validate_phyloseq <- function(phy, slots = c("sample_data", "tax_table")) {
  check_1 <- NULL
  if (!isS4(phy) & "phyloseq" %not_in% methods::is(phy)) {
    check_1 <- c("x" = "The file is not a Phyloseq object.")
  }
  
  check_2 <- NULL
  if ("sample_data" %in% slots) {
    test <- try(phyloseq::sample_data(phy), silent = TRUE)
    if (inherits(test, "try-error")) {
      check_2 <-
        c("x" = "The phyloseq object does not contain the sample_data slot.")
    } else {
      phyloseq::sample_data(phy) %>%
        tibble::as_tibble(rownames = "sample_id") %>%
        validate_sample_data()
    }
  }
    
  check_3 <- NULL
  if ("tax_table" %in% slots) {
    test <- try(phyloseq::tax_table(phy), silent = TRUE)
    if (inherits(test, "try-error")) {
      check_3 <-
        c("x" = "The phyloseq object does not contain the tax_table slot.")
    } else {
      phyloseq::tax_table(phy) %>%
        as.data.frame() %>%
        tibble::as_tibble(rownames = "otu_id") %>%
        validate_tax_table()
    }
  }
    
  check <- c(check_1, check_2, check_3)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while validating phyloseq object.", check))
  }
  phy
}

#' @rdname read_data
#' @keywords internal
#' @examples
#' phy <-
#'   system.file("extdata", "metaHIV_phy.rds", package = "dar") %>%
#'   read_phyloseq()
read_phyloseq <- function(file_path) {
  check_1 <- NULL
  if (!file.exists(file_path)) {
    check_1 <-
      c("x" = glue::glue("The file { crayon::red(file_path) } not exists."))
  }
  
  check_2 <- NULL
  if (!stringr::str_ends(stringr::str_to_lower(file_path), pattern = ".rds")) {
    check_2 <- c("x" = glue::glue(
      "The file must end with { crayon::red('.rds') } extensions.")
    )
  }
  
  check <- c(check_1, check_2)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while reading phyloseq object.", check))
  }
  
  phy <- readRDS(file_path)
  
  validate_phyloseq(phy)
  
  phyloseq::filter_taxa(phy, function(x)
    sum(x) != 0, TRUE)
}

#' @rdname read_data
#' @keywords internal
#' @tests testthat
#' ## File no exists and bad file extension
#' expect_error(read_file("asdfas.tx"))
read_file <- function(file_path, ext = c(".txt|.csv|.tsv")) {
  check_1 <- NULL
  if (!file.exists(file_path)) {
    check_1 <-
      c("x" = glue::glue("The file { crayon::red(file_path) } not exists."))
  }
  
  check_2 <- NULL
  if (!stringr::str_ends(file_path, pattern = ext)) {
    check_2 <-
      c("x" = glue::glue(
        "The file must end with", 
        " {crayon::red(c('.txt, .csv or .tsv'))} extensions."
      ))
  }
  
  check <- c(check_1, check_2)
  if (!is.null(check)) {
    rlang::abort(c("Prblem while reading text file.", check))
  }
  
  data.table::fread(file_path) %>% tibble::as_tibble()
}
