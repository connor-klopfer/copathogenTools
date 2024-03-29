% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolkit.R
\name{binarize_qpcr}
\alias{binarize_qpcr}
\title{Convert qPCR data to binary representation}
\usage{
binarize_qpcr(x)
}
\arguments{
\item{x}{the current column passed by \code{dplyr::transmute_all}. I beleive this
is simply a vector of numbers.}
}
\value{
the transformed vector to replace the passed \code{x} in the dataframe.
}
\description{
Function used by \code{transmute_all} to convert the qPCR data to binary,
representing ones or zeros. Each of the columns have different limits,
so it looks to see what values are below the max limit, since all
negative values would have the max number, either 35 or 50. This
function will not work if everybody is positive for a function. In
which case it will be off by one.
}
\note{
TODO: While this should work in 99.9% of cases, it has the potential
to break. I would like to use dplyr because it's fast and succinct, so
I need to find a way to pass a vector of variables to tansmute_all()
to allow for dynamic limits depending on the columns used.

The PROVIDE protocol give a value of -9 to missing values, if this changes, this code
will need to be updated.
}
