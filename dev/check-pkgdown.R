if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop(
    "Package 'pkgdown' is required. Install it with install.packages('pkgdown').",
    call. = FALSE
  )
}

if (!rmarkdown::pandoc_available()) {
  stop(
    paste(
      "Pandoc is required to validate pkgdown.",
      "Install Pandoc or set RSTUDIO_PANDOC to its binary directory."
    ),
    call. = FALSE
  )
}

pkgdown::check_pkgdown(".")
