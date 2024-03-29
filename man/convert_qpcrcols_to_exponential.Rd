% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolkit.R
\name{convert_qpcrcols_to_exponential}
\alias{convert_qpcrcols_to_exponential}
\title{Converts the qPCR cycle numbers to an exponential value. Accoring to the
Platts mills paper, every change in cyle corresponds with a doubling in
in the relative amount of amplicon. To do this with the qCPR data, we'll
set it to 2 to the paower of the difference minus 1}
\usage{
convert_qpcrcols_to_exponential(df, qpcr_cols)
}
\description{
df: the dataframe to convert
returns: The dataframe with all the qpcr columns converted to the continuous column.
     values are set to the gemotric mean.
}
