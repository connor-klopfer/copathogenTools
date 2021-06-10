#' These are all the tools need for data cleaning for the PROVIDE dataset.
#'

# fix_nines <- function(x){
#
# }


only_kids_with_both_stools <- function(d_f, original){
  #' Only get the kids from the PROVIDE dataset with both symptomatic and asymptomatic stools.
  #'
  #' @description We only want a subset of the provide dataset with participants
  kids_w_both <- original %>% filter(study == "provide") %>% group_by(participant, stool_type) %>% summarise(count = n())%>%
    ungroup() %>% group_by(participant) %>% summarise(count = n())%>%
    filter(count > 1) %>% select('participant') %>% unlist()

  final <- d_f %>% filter((study == 'maled') | (study == 'provide' & participant %in% kids_w_both))
  return(final)
}

only_kids_with_asym_samples <- function(d_f, original){
  kids_w_asym <- original %>% filter(study == "provide", stool_type == "Asymptomatic") %>%
    distinct(participant) %>% unlist()

  final <- d_f %>% filter((study == 'maled') | (study == 'provide' & participant %in% kids_w_asym))
  return(final)
}


provide_binary <- function(df, path_names, limit = 35, summarise_by = "sum"){
  #' Convert the the dataframe to a binary representation, where a one means
  #' the value was below the threshold while a zero means it was above the threshold.
  #'
  #' df: the dataframe of qPCR data
  #' path_names: the pathogen names in the dataframe. Used to exclude ID information
  #'     from binarization.
  path_names <- path_names[!is.na(path_names)]
  pathogen_columns <- which(names(df) %in% path_names)
  df <- df %>% mutate_at(pathogen_columns, na_if, -9)
  temp <- df[pathogen_columns] < 35 & df[pathogen_columns] > 0
  # print(df[which(setdiff(names(df), path_names) %in% names(df))])
  final <- cbind(df %>% select(one_of(c(setdiff(names(df), path_names)))),
                 temp)
  final <- final %>% select(one_of(c('participant', 'stool_type','age', "collection_date", path_names))) #%>%
  # group_by(stool_type, participant) %>%
  # summarise_each(sum)
  # print("PROVIDE NAMES")
  # # print(names(final))
  # print(c(setdiff(names(df), path_names)))
  # final <-
  return(final)
}


import_provide_dataset <- function(parent_dir = "../data/"){
  #' @title Import the qPCR data from PROVIDE
  #'
  #'
  #' @description Import the datasets for both symptomatic and asymptomatic datasets, then append
  #' them together with a variable sesignating them as symptomatic or asymptomatic.
  #'
  #'@param parent_dir String for the parent directory: defaults to \code{"../data"} which allows the function to be
  #'used in notebooks and scripts. This means the file is in a directory called \code{data}, which is a seperate
  #'directory from where the r notebook is kept. If you are storing in a spperate parent directory, then
  #'enter that in place of the default. If the file is in the same directory, enter \code{""}. If the file is in a
  #'parent directory, enter \code{".."}
  #'
  #'@return the qPCR data from PROVIDE, with both Symptomatic (Diarrhea) and Asymptomatic datasets combined into one dataset
  #'with a variable, "study" to distibguish them

  # Import the original data contained within the parent directory passed as the argument
  file_location <- paste(parent_dir, "TAC asymptomatic subset.csv", sep = "")

  # Asymptomatic Dataset.
  asymptomatic_raw <- read.csv(paste(parent_dir, "TAC asymptomatic subset.csv", sep = ""),
                               header = T, sep = ',', stringsAsFactors = F)

  # Symptomatic Dataset
  symptomatic_raw <- read.csv(paste(parent_dir, "TAC diarrheal YR1.csv", sep = "") ,
                              header = T, sep = ',', stringsAsFactors = F,
  )

  # Convert the dataset to a Date object
  asymptomatic_raw$specdt <- as.Date(asymptomatic_raw$specdt, format = "%m/%d/%Y")
  symptomatic_raw$specdt <- as.Date(symptomatic_raw$specdt, format = "%m/%d/%y")

  # Get the dataset as the variable name "key", with each row contains different
  # versions of the same variable axcroos different studies (columns)
  all_names <- get_target_pathogens(parent_dir)

  # Rename and set the names to more general terms
  asymptomatic_raw <- generalize_pathogen_names(asymptomatic_raw, all_names$asym_provide, all_names$actual)
  symptomatic_raw <- generalize_pathogen_names(symptomatic_raw, all_names$diar_provide, all_names$actual)

  # Make a new variable
  asymptomatic_raw$stool_type <- 'Asymptomatic'
  symptomatic_raw$stool_type <- "Diarrhea"

  full_set <- dplyr::bind_rows(asymptomatic_raw, symptomatic_raw)

  # Remove columns without a name: this means that pathogen is not needed.
  full_set[which(is.na(names(full_set)))] <- NULL


  return(full_set)
}





set_nas_zero <- function(df){
  df[is.na(df)] <- 0
  return(df)
}



