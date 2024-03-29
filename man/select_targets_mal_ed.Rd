% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maled_tools.R
\name{select_targets_mal_ed}
\alias{select_targets_mal_ed}
\title{Select target pathogens}
\usage{
select_targets_mal_ed(df, parent_dir = "../")
}
\arguments{
\item{df}{the dataframe with the larger set of variables.}

\item{parent_dir}{String for the parent directory: defaults to \code{"../data"} which allows the function to be
used in notebooks and scripts. This means the file is in a directory called \code{data}, which is a seperate
directory from where the r notebook is kept. If you are storing in a spperate parent directory, then
enter that in place of the default. If the file is in the same directory, enter \code{""}. If the file is in a
parent directory, enter \code{".."}}
}
\value{
The dataframe with only target pathogens.
}
\description{
Select the target pathogens from the pathogen list provided by Ross
and Dorothy. Selects only the target pathogens provided in the file "Pathogen Lists.xlsx"
provided in the Copathogen folder on the H://
}
\note{
This function relies on a CSV file 'pathogen_list_flagged.csv' which contaings names for the general terms,
provide, and mal-ed terms. Theres also a column indicating whether the variable is a pathogen or not. See
\code{\link{get_target_pathogens}}
}
