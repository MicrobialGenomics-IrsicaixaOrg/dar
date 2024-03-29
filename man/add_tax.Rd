% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recipe-class.R
\name{add_tax}
\alias{add_tax}
\alias{add_tax,Recipe-method}
\alias{add_tax,PrepRecipe-method}
\title{Adds taxonomic level of interest in the Recipe.}
\usage{
add_tax(rec, tax_info)

\S4method{add_tax}{Recipe}(rec, tax_info)

\S4method{add_tax}{PrepRecipe}(rec, tax_info)
}
\arguments{
\item{rec}{A \code{Recipe} object.}

\item{tax_info}{A character string of taxonomic levels that will be used in
any context.}
}
\value{
A \code{Recipe} object.
}
\description{
Adds taxonomic level of interest in the Recipe.
}
\examples{
data(metaHIV_phy)

## Define recipe
rec <-
  recipe(metaHIV_phy)

## add var info
rec <- add_tax(rec, tax_info = "Species")
rec

## add tax info to a prep-Recipe returns an error
data(test_prep_rec)
err <- testthat::expect_error(
  add_tax(test_prep_rec, tax_info = "Species")
)

err
}
