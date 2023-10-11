# For combining the datasets from Mal-Ed and PROVIDE.
# Only interested in a subset of the Provide Pathogens, and
# the intersection of that subset with the Mal-Ed data.

#'
#' get_both_datasets <- function(method = "simple", parent_dir = "../data/", bangladesh_only = F, specific_stool_type = NA){
#'   #' @title Import MAL-ED and PROVIDE qPCR Data
#'   #'
#'   #' @description This is the main function for importing the main dataset. Contains three cidfferent
#'   #'
#'   #' @param method: a string indicating the type of generalisation for the dataset. The available
#'   #'      methods are:
#'   #'      \itemize{
#'   #'      \item simple : whether that subject had that pathogen at anytime for
#'   #'      either asymptomatic or symptomatic diarhhea.
#'   #'      \item simple_noncondensed : the original dataset, with a binary representation
#'   #'      for all the qPCR data, and NOT integrated over the entire dataset. This means
#'   #'      each row is a stoll sample, and reach column is a pathogen
#'   #'      \item simple_sum: number of times that subject had that pathogen  with
#'   #'      either symptomatic or asymptomatic diarhhea.
#'   #'      \item binary_list : for each subject and pathogen, a boolean vector
#'   #'      representing one's for the subject having that pathogen at that
#'   #'      timepoint.
#'   #'      \item qpcr_list : a list of all the qpcr scores.
#'   #'      }
#'   #'
#'   #' @param parent_dir: a string as the parent directory of the file that's calling the method.
#'   #'      When running from a notebook, the parent directory is "../" and running from
#'   #'      a script, the parent directory is "" (empty)
#'   #'
#'   #' @source combining_datasets.R
#'   #'
#'   #' NOTE: This method raiases a
#'   #'
#'   #' @return: both dataframes with only the target pathogens for both.
#'   #'      Also contains that participant ID numbers, stool type, age, and collection date.
#'
#'   # General terms for naming the columns
#'   general_terms <- get_general_terms(parent_dir = parent_dir)
#'
#'   # Full Provide dataset
#'   full_set <- import_provide_dataset(parent_dir = parent_dir)
#'
#'   # All names as a dataframe
#'   name_set <- get_target_pathogens(parent_dir = parent_dir)
#'
#'   # Pathogen names
#'   pathogen_name_set <- name_set %>% filter(ispathogen > 0)
#'
#'
#'
#'   # Read in MAL-ED dataset
#'   mal_ed <- read.csv(paste(parent_dir, "maled_full_reduced.csv", sep = ""), header = T, stringsAsFactors = F)
#'
#'   # If the method is designated as a simple representation or binary non-condensed
#'   if(method == "simple" | method == "simple_noncondensed"){
#'
#'     mal_ed_final <-
#'       # Hamronize E.coli definitions
#'       harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
#'
#'       select_targets_mal_ed(parent_dir = parent_dir) %>%
#'       filter(age_days_eupath_0000579_ <= 371) %>%
#'       convert_maled_binary(ignore_cols = get_nonpathogen_cols(names(.), pathogen_name_set$mal_ed))
#'
#'     # print(get_nonpathogen_cols(names(mal_ed), pathogen_name_set$mal_ed), type = "labels")
#'
#'     provide_final <- provide_binary(full_set, pathogen_name_set$actual) %>%
#'       ungroup() %>%
#'       # set_nas_zero() %>%
#'       mutate_at(vars(pathogen_name_set$actual), binarize_summations)
#'
#'     provide_final <- generalize_pathogen_names(provide_final, name_set$asym_provide, name_set$actual)
#'     mal_ed_final <- generalize_pathogen_names(mal_ed_final, name_set$mal_ed, name_set$actual)
#'
#'   }else if(method == "simple_sum"){
#'     mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
#'       select_targets_mal_ed(parent_dir = parent_dir) %>%
#'       filter(age_days_eupath_0000579_ <= 371) %>%
#'       convert_maled_binary(ignore_cols = get_nonpathogen_cols(names(.), pathogen_name_set$mal_ed)) %>%
#'       group_by(participant_id, stool_type_gut_function_.eupath_0011745.) %>%
#'       summarise_at(pathogen_name_set$mal_ed[!is.na(pathogen_name_set$mal_ed)], sum, na.rm = T) %>%
#'       generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>% remove_na_columns() %>%
#'       ungroup() %>%
#'       arrange(desc(stool_type))
#'
#'     provide_final <- provide_binary(full_set, pathogen_name_set$actual) %>%
#'       ungroup() %>%
#'       set_nas_zero()
#'   }else if(method == "binary_list"){
#'     mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
#'       select_targets_mal_ed(parent_dir = parent_dir) %>%
#'       filter(age_days_eupath_0000579_ <= 371) %>%
#'       generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>%
#'       remove_na_columns() %>%
#'       ungroup() %>%
#'       condense_to_lists(pathogen_name_set$actual[!is.na(pathogen_name_set$mal_ed)], binarize_qpcr)
#'
#'
#'     provide_final <- condense_to_lists(full_set,
#'                                        pathogen_name_set$actual[!is.na(pathogen_name_set$actual)], binarize_qpcr) %>%
#'       mutate(participant = as.character(participant))
#'
#'   }else if(method == "qpcr_list"){
#'     mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
#'       select_targets_mal_ed(parent_dir = parent_dir) %>%
#'       filter(age_days_eupath_0000579_ <= 371) %>%
#'       generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>%
#'       remove_na_columns() %>%
#'       ungroup() %>%
#'       condense_to_lists(pathogen_name_set$actual[!is.na(pathogen_name_set$mal_ed)], return_same)
#'
#'
#'     provide_final <- condense_to_lists(full_set,
#'                                        pathogen_name_set$actual[!is.na(pathogen_name_set$actual)], return_same) %>%
#'       mutate(participant = as.character(participant))
#'   }
#'
#'   # NOTE: I think we might want the sample collection dates as well.
#'   # Instead of the sample collection date, we'll use the age of the individual, since that's more relative and all samples are
#'    # < 53 weeks of age .
#'
#'
#'
#'   provide_final$study <- "provide"
#'   mal_ed_final$study <- "maled"
#'
#'   provide_final$participant <- as.character(provide_final$participant)
#'
#'   if(bangladesh_only){
#'     mal_ed_final <- append_country_vars(mal_ed_final, parent_dir) %>% filter(country == 'Bangladesh')
#'   }
#'
#'
#'   # provide_final$collection_date <- as.Date(tolower(provide_final$collection_date), tryFormats = c("%m/%d/%Y", "%m/%d/%y"))
#'   mal_ed_final$collection_date <- as.Date(tolower(mal_ed_final$collection), format = "%d-%b-%y")
#'
#'   combined_df <- bind_rows(provide_final, mal_ed_final)
#'
#'   if(method == "simple" & is.na(specific_stool_type)){
#'     combined_df <- combined_df %>%
#'       select(-age) %>%
#'       group_by(participant, study) %>%
#'       mutate(stool_type = ifelse(stool_type == "Diarrhea", 1, 0)) %>%
#'       summarise_all(condense) %>%
#'       tidyr::pivot_longer(pathogen_name_set$actual,
#'                           names_to = "pathogen",
#'                           values_to = "lists") %>%
#'       mutate(sum_lists = lapply(lists, sum),
#'              isdiar = lapply(stool_type, sum),
#'              stool_type = ifelse(isdiar > 0, "Diarrhea", "Asymptomatic")) %>%
#'       mutate(is_present = (sum_lists > 0) * 1) %>%
#'       tidyr::pivot_wider(id_cols = c("participant", "study", "stool_type"),
#'                          names_from = "pathogen",
#'                          values_from = "is_present") %>% ungroup()
#'   }else if(method == "simple" & !is.na(specific_stool_type)){
#'     combined_df <- combined_df %>%
#'       select(-age) %>%
#'       filter(stool_type == specific_stool_type) %>%
#'       group_by(participant, study) %>%
#'       mutate(stool_type = ifelse(stool_type == "Diarrhea", 1, 0)) %>%
#'       summarise_all(condense) %>%
#'       tidyr::pivot_longer(pathogen_name_set$actual,
#'                           names_to = "pathogen",
#'                           values_to = "lists") %>%
#'       mutate(sum_lists = lapply(lists, sum),
#'              isdiar = lapply(stool_type, sum),
#'              stool_type = ifelse(isdiar > 0, "Diarrhea", "Asymptomatic")) %>%
#'       mutate(is_present = (sum_lists > 0) * 1) %>%
#'       tidyr::pivot_wider(id_cols = c("participant", "study", "stool_type"),
#'                          names_from = "pathogen",
#'                          values_from = "is_present") %>% ungroup()
#'   }
#'
#'   return(combined_df)
#' }


combine_like_pathogens <- function(d_f){
  #' @title Combine ETEC and EPEC variants
  #'
  #' @description After talking with Mami, it was recommended that we combine the etecs
  #' and rename typical epec as epec.
  #'
  #' @return dataframe with those pathogens combined, fewer rows.
  all_etec <- d_f %>% select(one_of(c("st_etec", "etec", "lt_etec"))) %>% rowSums(.)

  new_etec <- (all_etec > 0) * 1

  final <- cbind(
    d_f %>% select(-one_of(c("st_etec", "lt_etec", "etec"))),
    etec = new_etec
  ) %>%
    rename(
      epec = `typical epec`
    ) %>%
    select(
      - `atypical epec`
    )

  return(final)

}


get_both_datasets <- function(method = "simple", parent_dir = "../data/", bangladesh_only = F, specific_stool_type = NA, maled_version = "60m", age_limit = 372){
  #' @title Import MAL-ED and PROVIDE qPCR Data
  #'
  #' @description This is the main function for importing the main dataset. Contains three cidfferent
  #'
  #' @param method: a string indicating the type of generalisation for the dataset. The available
  #'      methods are:
  #'      \itemize{
  #'      \item simple : whether that subject had that pathogen at anytime for
  #'      either asymptomatic or symptomatic diarhhea.
  #'      \item simple_noncondensed : the original dataset, with a binary representation
  #'      for all the qPCR data, and NOT integrated over the entire dataset. This means
  #'      each row is a stoll sample, and reach column is a pathogen
  #'      \item simple_sum: number of times that subject had that pathogen  with
  #'      either symptomatic or asymptomatic diarhhea.
  #'      \item binary_list : for each subject and pathogen, a boolean vector
  #'      representing one's for the subject having that pathogen at that
  #'      timepoint.
  #'      \item qpcr_list : a list of all the qpcr scores.
  #'      }
  #'
  #' @param parent_dir: a string as the parent directory of the file that's calling the method.
  #'      When running from a notebook, the parent directory is "../" and running from
  #'      a script, the parent directory is "" (empty)
  #'
  #'
  #' @source combining_datasets.R
  #'
  #' NOTE: This method raiases a
  #'
  #' @return: both dataframes with only the target pathogens for both.
  #'      Also contains that participant ID numbers, stool type, age, and collection date.

  AGE_LIMIT <- age_limit

  # General terms for naming the columns in the PROVIDE dataset
  general_terms <- get_general_terms(parent_dir = parent_dir)

  # Full Provide dataset
  full_set <- import_provide_dataset(parent_dir = parent_dir)

  # All names as a dataframe
  name_set <- get_target_pathogens(parent_dir = parent_dir)

  # Pathogen names
  pathogen_name_set <- name_set %>% filter(ispathogen > 0)

  # Read in MAL-ED dataset
  if(maled_version == "60m"){

    reduced_filename <- "maled_samples_60m_reduced.csv"

    if(!file.exists(file.path(parent_dir,reduced_filename))){
      warning("Reduced file not found, taking original MAL-ED datafiles and creating reduced file....")
      reduced_full_maled_dataset(parent_dir)

    }
    mal_ed <- read.csv(file.path(parent_dir, "maled_samples_60m_reduced.csv"), header = T, stringsAsFactors = F)
    # maled_names <- ogen_name_set$maled_60m

    # id_cols <- c("Participant_Id", "Observation_Id", "Observation date [EUPATH_0004991]", "Age (days) [EUPATH_0000579]",
    #              "Country [OBI_0001627]", "Diarrheal stool [EUPATH_0011339]")
    id_columns <- c("Participant_Id", "Observation_Id", "Observation_date_EUPATH_0004991_", "Country_OBI_0001627_",
                    "Age_days_EUPATH_0000579_", "Stool_type_EUPATH_0010869_")

    names(mal_ed) <- fix_names(names(mal_ed))

    # print(names(mal_ed))


  }else if(maled_version == "24m"){
    message("Warning: even though you chose 24m, using the original release of the dataset. Synergy between these two,
            releases has not been confirmed.")
    mal_ed <- read.csv(file.path(parent_dir, "maled_full_reduced.csv"), header = T, stringsAsFactors = F)

    id_cols <- get_nonpathogen_cols(names(mal_ed), pathogen_name_set$mal_ed)
  }




  # If the method is designated as a simple representation or binary non-condensed
  if(method == "simple" | method == "simple_noncondensed"){

    mal_ed_final <- mal_ed %>%
      # Hamronize E.coli definitions
      harmonize_pathogen_data(maled_name_dictionary(maled_version)) %>%
      select_targets_mal_ed(parent_dir = parent_dir, vers = maled_version) %>%
      filter(age_days_eupath_0000579_ <= AGE_LIMIT) %>%
      convert_maled_binary(ignore_cols = id_cols)

    # print(get_nonpathogen_cols(names(mal_ed), pathogen_name_set$mal_ed), type = "labels")

    provide_final <- provide_binary(full_set, pathogen_name_set$actual) %>%
      ungroup() %>%
      # set_nas_zero() %>%
      mutate_at(vars(pathogen_name_set$actual), binarize_summations)

    provide_final <- generalize_pathogen_names(provide_final, name_set$asym_provide, name_set$actual)

    if(maled_version == "60m"){
      mal_ed_final <- generalize_pathogen_names(mal_ed_final, name_set$maled_60m, name_set$actual)
    }else if(maled_version == "24m"){
      warning("Did not select most up-to-date version: selecting first pahse version.")
      mal_ed_final <- generalize_pathogen_names(mal_ed_final, name_set$mal_ed, name_set$actual)
    }

  }else if(method == "simple_sum"){
    warning("Warning: This selection has not been updated for most recent release of MAL-ED")
    mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
      select_targets_mal_ed(parent_dir = parent_dir) %>%
      filter(age_days_eupath_0000579_ <= AGE_LIMIT) %>%
      convert_maled_binary(ignore_cols = get_nonpathogen_cols(names(.), pathogen_name_set$mal_ed)) %>%
      group_by(participant_id, stool_type_gut_function_.eupath_0011745.) %>%
      summarise_at(pathogen_name_set$mal_ed[!is.na(pathogen_name_set$mal_ed)], sum, na.rm = T) %>%
      generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>% remove_na_columns() %>%
      ungroup() %>%
      arrange(desc(stool_type))

    provide_final <- provide_binary(full_set, pathogen_name_set$actual) %>%
      ungroup() %>%
      set_nas_zero()
  }else if(method == "binary_list"){
    mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
      select_targets_mal_ed(parent_dir = parent_dir) %>%
      filter(age_days_eupath_0000579_ <= AGE_LIMIT) %>%
      generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>%
      remove_na_columns() %>%
      ungroup() %>%
      condense_to_lists(pathogen_name_set$actual[!is.na(pathogen_name_set$mal_ed)], binarize_qpcr)


    provide_final <- condense_to_lists(full_set,
                                       pathogen_name_set$actual[!is.na(pathogen_name_set$actual)], binarize_qpcr) %>%
      mutate(participant = as.character(participant))

  }else if(method == "qpcr_list"){
    mal_ed_final <- harmonize_pathogen_data(mal_ed, maled_name_dictionary()) %>%
      select_targets_mal_ed(parent_dir = parent_dir) %>%
      filter(age_days_eupath_0000579_ <= AGE_LIMIT) %>%
      generalize_pathogen_names(name_set$mal_ed, name_set$actual) %>%
      remove_na_columns() %>%
      ungroup() %>%
      condense_to_lists(pathogen_name_set$actual[!is.na(pathogen_name_set$mal_ed)], return_same)


    provide_final <- condense_to_lists(full_set,
                                       pathogen_name_set$actual[!is.na(pathogen_name_set$actual)], return_same) %>%
      mutate(participant = as.character(participant))
  }

  # NOTE: I think we might want the sample collection dates as well.
  # Instead of the sample collection date, we'll use the age of the individual, since that's more relative and all samples are
  # < 53 weeks of age .



  provide_final$study <- "provide"
  mal_ed_final$study <- "maled"

  provide_final$participant <- as.character(provide_final$participant)

  if(bangladesh_only){
    if(maled_version == "60m"){
      mal_ed_final <- mal_ed_final %>% filter(country_obi_0001627_ == 'Bangladesh')


    }else if(maled_version == "24m"){
      mal_ed_final <- append_country_vars(mal_ed_final, parent_dir) %>% filter(country == 'Bangladesh')
      # provide_final$collection_date <- as.Date(tolower(provide_final$collection_date), tryFormats = c("%m/%d/%Y", "%m/%d/%y"))
      mal_ed_final$collection_date <- as.Date(tolower(mal_ed_final$collection), format = "%d-%b-%y")
    }
  }

  if(maled_version == "60m"){

    mal_ed_final <- mal_ed_final %>%
      rename(
        # observation_date = observation_date_eupath_0004991_,
        country = country_obi_0001627_) %>% mutate(collection_date = as.Date(collection_date),
                                                   stool_type = ifelse(stool_type == "Monthly", "Asymptomatic", stool_type))

  }



  combined_df <- bind_rows(provide_final, mal_ed_final)

  if(method == "simple" & is.na(specific_stool_type)){
    combined_df <- combined_df %>%
      select(-age) %>%
      group_by(participant, study) %>%
      mutate(stool_type = ifelse(stool_type == "Diarrhea", 1, 0)) %>%
      summarise_all(condense) %>%
      tidyr::pivot_longer(pathogen_name_set$actual,
                          names_to = "pathogen",
                          values_to = "lists") %>%
      mutate(sum_lists = lapply(lists, sum),
             isdiar = lapply(stool_type, sum),
             stool_type = ifelse(isdiar > 0, "Diarrhea", "Asymptomatic")) %>%
      mutate(is_present = (sum_lists > 0) * 1) %>%
      tidyr::pivot_wider(id_cols = c("participant", "study", "stool_type"),
                         names_from = "pathogen",
                         values_from = "is_present") %>% ungroup()
  }else if(method == "simple" & !is.na(specific_stool_type)){
    combined_df <- combined_df %>%
      select(-age) %>%
      filter(stool_type == specific_stool_type) %>%
      group_by(participant, study) %>%
      mutate(stool_type = ifelse(stool_type == "Diarrhea", 1, 0)) %>%
      summarise_all(condense) %>%
      tidyr::pivot_longer(pathogen_name_set$actual,
                          names_to = "pathogen",
                          values_to = "lists") %>%
      mutate(sum_lists = lapply(lists, sum),
             isdiar = lapply(stool_type, sum),
             stool_type = ifelse(isdiar > 0, "Diarrhea", "Asymptomatic")) %>%
      mutate(is_present = (sum_lists > 0) * 1) %>%
      tidyr::pivot_wider(id_cols = c("participant", "study", "stool_type"),
                         names_from = "pathogen",
                         values_from = "is_present") %>% ungroup()
  }

  return(combined_df)
}


import_complete_datasets <- function(method = "simple_noncondensed", bangladesh_only = T, parent_dir = "../data/",
                                     use_pathogen_subset = TRUE, age_limit = 372){
  #' @title Import Complete Cases of Original Data.
  #'
  #' @description Import only the complete cases of the original dataset for both MAL-ED and PROVIDE. Similar to
  #' \code{\link{get_both_datasets}} except only the complete cases are returned. Parameters are the same as
  #' \code{\link{get_complete_cases}}. Also include the age bins that certain observations fall into, the season, and the
  #' months. See \code{\link{get_pathogen_selection_selection_targets}} for details on the pathogens used.
  #'
  #'
  #'   @param method: a string indicating the type of generalisation for the dataset. The available
  #'      methods are:
  #'      \itemize{
  #'      \item simple : whether that subject had that pathogen at anytime for
  #'      either asymptomatic or symptomatic diarhhea.
  #'      \item simple_noncondensed : the original dataset, with a binary representation
  #'      for all the qPCR data, and NOT integrated over the entire dataset. This means
  #'      each row is a stoll sample, and reach column is a pathogen
  #'      \item simple_sum: number of times that subject had that pathogen  with
  #'      either symptomatic or asymptomatic diarhhea.
  #'      \item binary_list : for each subject and pathogen, a boolean vector
  #'      representing one's for the subject having that pathogen at that
  #'      timepoint.
  #'      \item qpcr_list : a list of all the qpcr scores.
  #'      }
  #'
  #' @param parent_dir: a string as the parent directory of the file that's calling the method.
  #'      When running from a notebook, the parent directory is "../" and running from
  #'      a script, the parent directory is "" (empty)
  #'
  #' @param bangladesh_only: A boolean: whether only to include observations from MAL-ED collected in Bangladesh.
  #' Does not affect PROVIDE samples.
  #'
  #' use
  #'
  #'
  #' @source combining_datasets.R
  #'
  #' NOTE: This method raiases a
  #'
  #' @return: both dataframes with only the target pathogens for both.
  #'      Also contains that participant ID numbers, stool type, age, and collection date.

  # Get the pathogen targets to use for selection.
  path_select <- get_pathogen_selection_targets(parent_dir)
  print(age_limit)

  # Get the final the dataset, only get the kids from PROVIDE with corresponding asympotmatic samples.
  final <- get_both_datasets(method = method, bangladesh_only = bangladesh_only, parent_dir = parent_dir, age_limit = age_limit) %>%
    only_kids_with_asym_samples(.,.) %>%
    mutate(
      bin_12 = floor(age/(372/12)),
      bin_26 = floor(age/(372/26)),
      bin_6 = floor(age/(372/6))) %>%
    assign_season() %>%
    mutate(month = as.numeric(format(month_day, "%m")) - 1,
           observation_id = ifelse(study == "provide", paste(participant, collection_date, sep = "_"), observation_id)) %>%
    select("participant", "study", "stool_type", "age", "collection_date", "observation_id", "country", "bin_12", path_select) %>%
    mutate(country = case_when(study == "provide" ~ "Bangladesh",
                             study != "provide" ~ country)) %>%
    filter(complete.cases(.))

  return(final)

}


get_pathogen_selection_targets <- function(parent_dir = "../data/", use_subset = TRUE){
  #' @title Get the Target Pathogens
  #'
  #' @description Get a character vector of the target pathogens to use within a \code{dplyr::select} function.
  #' As the option to exclude pathogens where the data is MNAR or inconsistent between studies.
  #'
  #' @param parent_dir: character of the filepath to use for the folder containing the file as the dataset for the different variable neamse
  #'
  #' @param use_subset: boolean, whether to use the subset of pathogens to elimiante MNAR data and OPV subsets. Eliminated pathogens include
  #' all OPVs, Girardia spp., C.Jejuni/coli specific marker, and the cryptolib parvum/hominus specific markers.

  name_set <- get_target_pathogens(parent_dir = parent_dir)

  # Remove OPV samples
  name_set_no_opv <- name_set[!grepl("opv", name_set$actual),]

  pathogen_selection <- name_set_no_opv$actual[name_set_no_opv$ispathogen > 0] %>% gsub("\\s+$", "", .)

  path_select <- setdiff(name_set_no_opv$actual[name_set_no_opv$ispathogen > 0] %>%
                           gsub("\\s+$", "", .), c("crypto_lib parvum", "crypto_lib hominis", "giardia", "c.jejuni/coli"))

  return(path_select)

}

