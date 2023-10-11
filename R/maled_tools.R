#' This is the toolkit for the mal-ed dataset.
#' USed for tools that are specific to the mal-ed set.

append_country_vars <- function(maled, parent_dir = "../data/"){
  countries <- get_country_ids(parent_dir)
  maled$participant <- as.character(maled$participant)
  countries$participant <- as.character(countries$participant)
  final <- dplyr::left_join(maled, countries, by = c("participant"))
  return(final)
}


clean_maled <- function(raw_df){
  #' Clean the MAL-ED dataset
  #'
  #' @description Parent function for importing the mal-ed. Imports that dataframe ant
  #' cleans the dataset for only the relvent data.
  #'
  #' @param raw_df: the dataframe directly from the
  #'
  #' @return: The mal-ed dataset. The

  # Change the names to a single term.
  names(raw_df) <- fix_names(raw_df)

  reduced <- reduce_mal_ed(raw_df)

  names(reduced) <- fix_names(reduced)

  final <- eliminate_empty_columns(reduced, ignore_columns = c(1:6, 80, 81))

  return(final)

}


convert_maled_binary <- function(d_f, ignore_cols, method = "simple"){
  #' Convert MAL-ED to binary matrix
  #'
  #' @description Converts the Mal-ed qpcr Data depending if they are below the threhsold value.
  #' Reight now defaults calls binarize_qpcr() using a dplyr transmute_all function.
  #'
  #' @param df: the dataframe of the MAL-ED dataset
  #' @param ignore_cols: columns to ignore for binarisation, these are ID columns such as ID and
  #' age that you do not want to mask.
  #'
  #' @return the MAL-ED dataframe as a binary representation.
  #'
  #' TODO: This should be in the mal-ed only toolbox, since it is only
  #' used for the mal-ed dataset

  ignore_cols <- setdiff(names(d_f), find_names_matching(names(d_f), "ct_value"))


  ignored_cols <- d_f %>% select(ignore_cols)
  # ignored_cols <- d_f[ignore_cols]
  # df[ignore_cols] <- NULL

  qpcr_cols <- d_f %>% select(-ignore_cols)

  if(method == "simple"){
    return(cbind(ignored_cols,
                 qpcr_cols %>% transmute_all(binarize_qpcr) %>% transmute_all(binarize_summations)))
  }else if(method == "simple_sum"){
    return(cbind(ignored_cols,
                 qpcr_cols %>% transmute_all(binarize_qpcr)))
  }
}


get_country_ids <- function(parent_dir){
  #' Get the Country Variable for the maled Dataset
  #'
  #'@description The Country variable for the mal-ed dataset is not aligned with the qPCR data
  #'and is missing in some cases. This extracts
  countries <- read.csv(paste(parent_dir, "country_vars.csv", sep = ""))
  names(countries) <- c('participant', 'country')
  return(countries)

}


reduce_mad_ed <- function(df){
  #' Remove uneeded columns in MAL-ED
  #'
  #' @description Reduce the mal_ed data frame to the Ct values, the ID columns,
  #' and remove any rows that are not designated as 'Diarrhea' Or asymptomatic.
  #' Uses Left Join to only obervations,  that have a mtaching obervation id.
  #'
  #' @param df: the MAL-ED dataframe with all columns.
  #'
  #' @return MAL-ED dataframe with non-id and non pathogen columns removed.
  final <- left_join(
    data.frame(cbind(df[, 1:6], # Append the identification columns.
                     df %>%
                       select(one_of(find_intersection_of_terms(., c("_ct_", "_value_")))))), # All the qcpr data
    cbind(df %>%
            select(one_of(find_intersection_of_terms(., c("Stool", "gut")))), # Gut function (Sym, Asym)
          df %>%
            select(one_of(find_intersection_of_terms(., c("Observation", "Id")))), # Observation ID
          df %>%
            select(one_of(find_names_matching(names(.), "Age_\\(days\\)")))) %>% # Age in days
      filter(!is.na(.[[find_intersection_of_terms(df, c("Stool", "gut"))]])), # Remove NA gut functions variables.
    by = find_intersection_of_terms(df, c("Observation", "Id")))
  return(final)
}


select_targets_mal_ed <- function(d_f, parent_dir = "../", vers = "60m"){
  #' @title Select target pathogens
  #'
  #' @description Select the target pathogens from the pathogen list provided by Ross
  #' and Dorothy. Selects only the target pathogens provided in the file "Pathogen Lists.xlsx"
  #' provided in the Copathogen folder on the H://
  #'
  #' @param df the dataframe with the larger set of variables.
  #'
  #'@param parent_dir String for the parent directory: defaults to \code{"../data"} which allows the function to be
  #'used in notebooks and scripts. This means the file is in a directory called \code{data}, which is a seperate
  #'directory from where the r notebook is kept. If you are storing in a spperate parent directory, then
  #'enter that in place of the default. If the file is in the same directory, enter \code{""}. If the file is in a
  #'parent directory, enter \code{".."}
  #'
  #'@note This function relies on a CSV file 'pathogen_list_flagged.csv' which contaings names for the general terms,
  #'provide, and mal-ed terms. Theres also a column indicating whether the variable is a pathogen or not. See
  #'\code{\link{get_target_pathogens}}
  #'
  #' @return The dataframe with only target pathogens.

  # Get the pathogen names
  targets <- get_target_pathogens(parent_dir = parent_dir)

  # Rename to lowercase
  names(d_f) <- tolower(names(d_f))

  if(vers == "60m"){
    id_identifiers <- c("participant_id", "observation_id", "observation_date_eupath_0004991_", "country_obi_0001627_", "age_days_eupath_0000579_",
                        "stool_type_eupath_0010869_")
  }else if(vers == "24m"){
    id_identifiers <- c(
      "participant_id", "observation_id", "age_at_last_contact_days_eupath_0011147_", "last_contact_date_eupath_0011146_",
      "stool_type_gut_function_eupath_0011745_", "age_days_eupath_0000579_", "sample_collection_date_obi_0001619_")
  }
  # Get all the names that are not ct names


  id_cols <- d_f %>% select(all_of(id_identifiers))

  # Find all the names containing ct data
  ct_cols <- d_f %>% select(find_names_matching(names(d_f), "_ct_value"))


  # Get the pathogen columns
  targets <- targets %>% filter(ispathogen > 0)

  if(vers == "60m"){
    ct_cols <- ct_cols %>% select(one_of(targets$maled_60m))
  }else if(vers == "24m"){
    warning("Be advised: Using older version of MAL-ED")
    ct_cols <- ct_cols %>% select(one_of(targets$mal_ed))
  }else{
    stop("Impropoer version given. Options are 24m and 60m")
  }

  final <- cbind(id_cols, ct_cols)

  return(final)
}

