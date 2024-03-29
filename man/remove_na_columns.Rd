% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolkit.R
\name{remove_na_columns}
\alias{remove_na_columns}
\title{When generalizing the column names, there's columns returned with NA
as the column name. These are not part of the 'desired' pathogens.
Remove these columns from the dataset.}
\usage{
remove_na_columns(df)
}
\description{
df: the dataframe to reduce
returns: the dataframe with all the NA columns removed.
}
