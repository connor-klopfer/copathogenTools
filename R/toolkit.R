# Tool kit for functions



binarize <- function(x){
  #' @title Convert Yes/No to 1's and 0's
  #'
  #' @description Convert variables with a "Yes" or "No" entry into 1's and 0's.
  #' Makes things easy for co-occurance matrices and math down the road.
  #'
  #' @param x a column in a dataframe
  #'
  #' @return the column as a binary representation.
  #'
  #' @examples sample_column <- c("Yes", "No, "Yes")
  #'
  #' > binarize(sample_column)
  #' > c(1, 0, 1)
  #'
  x[which(is.na(x))] <- "No"
  x <- factor(x, levels = c("No", "Yes"), labels = c(F,T))
  x <- as.integer(x) - 1
  return(x)
}

binarize_qpcr <- function(x){
  #' @title Convert qPCR data to binary representation
  #'
  #' @description Function used by \code{transmute_all} to convert the qPCR data to binary,
  #' representing ones or zeros. Each of the columns have different limits,
  #' so it looks to see what values are below the max limit, since all
  #' negative values would have the max number, either 35 or 50. This
  #' function will not work if everybody is positive for a function. In
  #' which case it will be off by one.
  #'
  #' @note TODO: While this should work in 99.9% of cases, it has the potential
  #' to break. I would like to use dplyr because it's fast and succinct, so
  #' I need to find a way to pass a vector of variables to tansmute_all()
  #' to allow for dynamic limits depending on the columns used.
  #'
  #' @note The PROVIDE protocol give a value of -9 to missing values, if this changes, this code
  #' will need to be updated.
  #'
  #' @param x the current column passed by \code{dplyr::transmute_all}. I beleive this
  #' is simply a vector of numbers.
  #'
  #' @return the transformed vector to replace the passed \code{x} in the dataframe.
  max_value <- max(x, na.rm = T)
  final <- x
  # UPDATE: If the final is < 0, that designates a missing value, so we'll keep
  # it as NA.

  # If the value is -9, this is a missing value in the PROVIDE protcol
  final <- ifelse(final < 0, NA, ifelse(final < max_value, 1, 0))

  # If all the subjects are positive, this method doesn't work, and something else must
  # used. This is a rare occurance, so I'll adapt accordingly if it happens.
  if(sum(final, na.rm = T) == sum(!is.na(x))){
    stop("ERROR: All Subjects are positive, binarize_qpcr is broken.")
  }
  return(final)
}


binarize_summations <- function(x){
  #' @title Binarize using dplyr functions.
  #'
  #' @description Base function called with \code{summarise_at}, to convert the summations to 1 or zero. This is used for
  #' intergrating pathogen incidence for a subject over the whole study. Without this binarisation, simply summing the incidence
  #' can lead to sampling error.
  #'
  #' @param x the column passed by summarise_at as a vector from the column in the df
  #' @return  the boolean (as 1's and 0's) vector representation is x > 0
  #'
  #' @seealso dplyr::summarise_at
  return(as.numeric(x > 0))
}


# To be removed, commenting out to make sure there's no dependency: it shouldn't
# there's no return value.
# build_co_matrix <- function(df){
#    Build a cube for the co-occurance matrix.
# }


collapse_array <- function(my_array, method = "same_timepoint_binary"){
  #'@title Collapse a 3 dimensional array to a 2 dimensional matrix.
  #'
  #'@description To collapse a 3 dimensional array into a single matrix depending on
  #'the preferred method. Uses to build the confidence intervals and convert the
  #'df with lists in them to rewiriable values.
  #'
  #' @param my_array The 3 dimensional array, with rows for kids kids columns equal
  #' to the number of pathogens and slices (z) equal to the number of timepoints,
  #'
  #' @param method the type of collapse the following are available.
  #'      same_timepoint_binary: Do the two pathoges accor togther in the same timepoint at any time?
  #'      TODO: Should we include same age?
  #'      TODO: Maybe one in front, one in the back?
  #' returns: the collapsed 2D matrix. with the same dinemsions as the DFs of lists.

  # The new array, we'll be doing matrix mutiplication so I'll reshape the new matrix accordingly.
  my_array_dims <- dim(my_array)

  co_array <- array(data = NA, dim = c(my_array_dims[1], my_array_dims[2], my_array_dims[3]))

  # We have to set the NA's as zero, otherwise the matrix math won't work
  my_array[is.na(my_array)] <- 0

  #Co-occurance slices at each timepoint
  for(z in 1:dim(my_array)[3]){
    # co_array[,,z] <- t(my_array[,,z]) %*% my_array[,,z]
    co_array[,,z] <- t(my_array[,,z]) %*% my_array[,,z]
  }

  #Collapse the all the values by summations.
  collapsed <- rowSums(co_array, na.rm = FALSE, dims = 2)

  # Boolean, is there at least one co-occurance?
  final <- collapsed > 0

  return(final)

}


condense <- function(x){
  #' To condense the dataframe to a list.
  #' x: the grouped variable
  #' return: the dataframe with a list of all the occurances
  return(list(x))
}


condense_to_lists <- function(df, qpcr_cols, fun){
  #' Converts the qPCR cycle numbers to an exponential value. Accoring to the
  #' Platts mills paper, every change in cyle corresponds with a doubling in
  #' in the relative amount of amplicon. To do this with the qCPR data, we'll
  #' set it to 2 to the paower of the difference minus 1
  #'
  #' df: the dataframe to convert
  #' returns: The dataframe with all the qpcr columns converted to the continuous column.
  #'      values are set to the gemotric mean.
  final <- df
  final[qpcr_cols] <- df %>% transmute_at(vars(qpcr_cols), fun)
  final <- final %>%
    select(one_of(c('participant', 'stool_type', 'age', qpcr_cols))) %>%
    # group_by(stool_type, sid) %>% summarise_all(psych::geometric.mean) %>% ungroup()
    group_by(participant) %>% summarise_all(condense) %>% ungroup()
  return(final)
}


convert_matrix <- function(df){
  # Converts a dataframe to a matrix and records the names from the origianl
  # dataframe into a vector. Then combines the matrix and the name vector into
  # a list to return from the
  # df: dataframe
  # Returns: a matrix representation of the dataframe
  new_matrix <- matrix(nrow = nrow(df), ncol = length(df))
  name_list <- vector(mode = "character", length(length(names(df))))
  for(i in 1:length(df)){
    name_list[i] <- names(df)[i]

    # NOTE: Is this a temp solution, we might have to do something different with NA
    # values other than set them to zero
    df[,i][is.na(df[,i])] <- 0

    new_matrix[,i] <- unlist(df[,i])
    # new_matrix[,i] <- df[,i]
  }

  matrix_name_list <- list()
  matrix_name_list[[1]] <- new_matrix
  matrix_name_list[[2]] <- name_list
  return(matrix_name_list)
}


convert_matrix_lists <- function(d_f){
  # Converts a dataframe to a matrix and records the names from the origianl
  # dataframe into a vector. Then combines the matrix and the name vector into
  # a list to return from the
  # d_f: dataframe
  # Returns: a matrix representation of the dataframe
  new_matrix <- array(data = NA, dim = c(nrow(d_f), length(d_f), 20))
  name_list <- names(d_f)

  # TODO: lets add in another function here to preporcess the array. One that
  # finds the maximum nuber of slices needed and sets any non lists as lists


  for(r in 1:nrow(d_f)){
    for(c in 1:ncol(d_f)){
      if(!is.null(d_f[[c]][[r]])){
        new_matrix[r, c,c(1:length(d_f[[c]][[r]]))] <- as.numeric(d_f[[c]][[r]])
      }else{
        new_matrix[r, c,c(1:length(d_f[[1]][[r]]))] <- (rep(NA, length(d_f[[1]][[r]])))
      }
    }
  }
  new_matrix <- remove_empty_slices(new_matrix)
  #TODO: Need a dynamic way to control for the number of slices in the array
  #TODO:


  matrix_name_list <- list()
  matrix_name_list[[1]] <- new_matrix
  matrix_name_list[[2]] <- name_list
  return(matrix_name_list)
}


convert_qpcrcols_to_exponential <- function(df, qpcr_cols){
  #' Converts the qPCR cycle numbers to an exponential value. Accoring to the
  #' Platts mills paper, every change in cyle corresponds with a doubling in
  #' in the relative amount of amplicon. To do this with the qCPR data, we'll
  #' set it to 2 to the paower of the difference minus 1
  #'
  #' df: the dataframe to convert
  #' returns: The dataframe with all the qpcr columns converted to the continuous column.
  #'      values are set to the gemotric mean.
  final <- df
  final[qpcr_cols] <- df %>% select(one_of(qpcr_cols)) %>% transmute_all(qpcr_continuous)
  final <- final %>% select(one_of(c('sid', 'stool_type', qpcr_cols))) %>%
    # group_by(stool_type, sid) %>% summarise_all(psych::geometric.mean) %>% ungroup()
    group_by(stool_type, sid) %>% summarise_all(sum) %>% ungroup()
  return(final)
}


count_empties <- function(x){
  # Counts the number of empty entries in the column of the
  # dataframe. Used for determining if a column has all empyt entries.
  return(sum(is.na(x)))
}

# Count the number of levels when you factorize a column
count_levels <- function(x){
  return(length(levels(factor(x))))
}

# Filters out empty columns, uses the count_empties() function.
# If the column is empty then the function returns a value that's
# equal to the number of rows in the dataset. Therefore, a column
# with even one entry will pass
# df: dataframe in question.
# Returns: the
eliminate_empty_columns <- function(df, ignore_columns){
  # I'm removing these since apply functions are slow. We'll use dplyr instead.
  # counts <- c(apply(df, 2, count_empties)) < nrow(df)
  if(!is.na(ignore_columns)){
    ignored_columns <- df[ignore_columns]
    df <- df[-ignore_columns]
  }
  df$Tally <- rowSums(df, na.rm = T)
  if(!is.na(ignore_columns)){
    temp_df <- cbind(ignored_columns[which(df$Tally > 0),], df[df$Tally > 0,])
  }else{
    temp_df <- df[df$Tally > 0]
  }
  temp_df$Tally <- NULL
  return(temp_df)
}

eliminate_identical_columns <- function(df){
  # Some columns have the exact same entry for all observations in the
  # dataset. Eliminate those columns, they have nothing useful at the moment.
  # df: dataset
  # Returns: Datasets with non-identical entries.
  counts <- c(apply(df, 2, count_levels)) > 1
  return(df[,counts])
}

extract_binary_matrix <- function(df){
  # The Mal-ED dataset has a lot of variables that show results, then cumulative sums, and
  # results.
  # df : dataframe
  # returns the binary representation of the result variables
  return(df[grep("_by", names(df))] %>%
           .[-grep("result", names(.))] %>%
           transmute_all(binarize))
}

# Searching for multiple terms
# TODO: Possibly cut out loop
find_intersection_of_terms <- function(df, terms_list, use_perl = F){
  intersections <- find_names_matching(names(df), terms_list[1])
  for(x in 2:length(terms_list)){
    intersections <- dplyr::intersect(intersections, find_names_matching(names(df), terms_list[x], use_perl = use_perl))
  }
  return(intersections)
}

# Used for finding variables in the dataset that match a criteria.
# USed mostly for speed than typing in the entire grep function.
# col_names: the names of the dataset in question
# pattern: pattern to look for
# returns: the names of the variables with mathcing pattern
find_names_matching <- function(col_names, pattern, indices = F, use_perl = F, ignore_case = F){
  if(indices){
    return(grep(pattern, col_names, ignore.case = ignore_case, perl = use_perl))
  }else{
    return(col_names[grep(pattern, col_names, ignore.case = ignore_case, perl = use_perl)])
  }

}

fix_names <- function(col_name){
  # Eliminated spaces out of the column names and replace them with underscores.
  # Keeps the nameing convention consistent. Utilizes Regex without PERL
  # col_name: either a vector of the column names of a single column name

  # returns: a vector of the column names to be assigned outside of the function.

  return(gsub("(\\.|   |\\,|\\[|\\]|\\(|\\))", "_", col_name) %>%
           gsub("( )", "_", .) %>% gsub("(_{2,})", "_", .))
}


compare_degree_dist <- function(matrix1, matrix2, by = "row"){
  #' Compare the degree distributions between two matrices. Uses the get_deg_distr
  #' function to find the distr. Then returns the distributions which are not the
  #' same.
  #'
  #' matrix1: matirx of nXm
  #' matrix2: matrix of the same dimensions
  #'
  #' returns: a vector of the degree that are mismatched.
  distr1 <- get_deg_distr(matrix1, by)
  distr2 <- get_deg_distr(matrix2, by)

  return(distr1[which(distr1[,2] != distr2[,2]),1])
}


get_deg_distr <- function(matrix1, by){
  #' Get the degree distribution for two matrices. USes DPLYR and returns a matrix
  #' of the distribution.
  #'
  #' returns: matrix of the distrbutios, with the first col the degree of the nodes
  #'      and the second col is the counts.
  #'
  if(by == "row"){
    distr1 <- rowSums(matrix1, na.rm = T) %>% plyr::count()
  }else if(by == "col"){
    distr1 <- colSums(matrix1, na.rm = T) %>% plyr::count()
  }else{
    stop("Improper Dimension")
  }
  return(distr1)
}


get_general_terms <- function(parent_dir = "../data/", alt_filename){
  #' @title Import general variable names to use on both datasets
  #'
  #'
  #' @description Imports the file: \code{provide_pathogen_names.csv}, which contains the
  #' symptomatic and asymptomatic pathogen names. To be used on the respective
  #' datasets.
  #'
  #'@param parent_dir String for the parent directory: defaults to \code{"../data"} which allows the function to be
  #'      used in notebooks and scripts. This means the file is in a directory called \code{data}, which is a seperate
  #'      directory from where the r notebook is kept. If you are storing in a spperate parent directory, then
  #'      enter that in place of the default. If the file is in the same directory, enter \code{""}. If the file is in a
  #'      parent directory, enter \code{".."}
  #'@return a dataframe where each row contains the variable's general terms, the version in the maled, asymptomatic,
  #'and symptomatic PROVIDE datasets.

  # Get the original dataset
  all_names <- read.csv(file.path(parent_dir, "provide_pathogen_names.csv"), header = T, stringsAsFactors = F)

  # The asymptomatic name/general name
  asym_names <- all_names[,1:2]
  names(asym_names) <- c("asym_provide", "actual")
  sym_names <- all_names[, 3:4]
  names(sym_names) <- c("diar_provide", "actual")

  # Full join by actual: might not be needed but ensures the correct alignment
  all_names <- full_join(asym_names, sym_names, by = "actual")
  all_names <- all_names %>% mutate_at(c("asym_provide", "diar_provide"), tolower)
  return(all_names)
}


# For unique encoding of the pathogen clusters. Generates a vector of numbers to
# that get multiplied by the actual binary vector to generate an integer encoding.
generate_binary_vec <- function(vec_length){
  binary_scaffold <- vector(mode = "numeric", length = 36)
  for(i in 0:35){
    binary_scaffold[i + 1] <- 2^(i)
  }
  return(binary_scaffold)
}

generate_binary_num <- function(df){
  binary <- generate_binary_vec(36)
  binary_values <- vector(mode = "numeric", length = nrow(df))
  df_m <- as.matrix(df[1:36])

  counts <- df_m %*% as.matrix(binary)
  return(counts)
}

generate_co_occurance <- function(df){
  m <- as.matrix(df)
  return(t(m) %*% m)
}


generalize_pathogen_names <- function(df, name_set_specific, name_set_general){
  #' @title Normalise variable names for dataset.
  #'
  #' @description For normalizing the names between the two datasets. Matches the
  #' original name in the dataset with the general term provided in the CSV "key"
  #' file.
  #'
  #' @note
  #' IMPORTANT: The names need to be matched up in the dataframe for this to work.
  #' If the name in the original file or the desired name of of the general variable
  #' name changes, then that name needs to be changed in the respective source file containing the
  #' variable name keys. See \code{\link{get_target_pathogens}} for source file name and how keys are obtained.
  #'
  #'@seealso \code{\link{get_target_pathogens}}
  #'
  #' @param df The dataframe to be renamed
  #' @param name_set_specific A list of the names in the original dataset, passed as a column
  #' from a dataframe.
  #'
  #' @param name_set_general a list of the general terms to be used, passed as a column from a
  #' dataframe
  #'
  #' @return The dataframe with the names set to the more general terms.

  # Convert all to lower case
  names(df) <- tolower(names(df))

  # Names from the original dataset
  name_set_specific <- tolower(name_set_specific)
  # General names
  name_set_general <- tolower(name_set_general)

  # The names of the dataset
  df_col_names <- names(df)

  # Loop through each name in the specific name set
  for(name in name_set_specific){
    # Get the index of the names that mathcies the corrent name in the loop
    matching_idx <- which(df_col_names == name)
    # If there's no matches then skip and go to the next name
    if(length(matching_idx) == 0 | is.na(name)){
      next()
      # IF there's multiple names, then we have a hash collision.
    }else if(length(matching_idx) > 1){
      stop("HASH COLLISION ERROR: Multiple mathcing names")
    }

    # Rename the variable as the general term for the specific term in the loop.
    names(df)[matching_idx] <- name_set_general[which(df_col_names[matching_idx] == name_set_specific)]
  }

  return(df)
}


get_nonpathogen_cols <- function(all_cols, pathogen_cols, type = 'index'){
  #' Retreive the non-pathogen columns to be used as the ID columns. Can either
  #' return the index of the ID cols or the actual labels.
  #'
  #' all-cols: all the names of the dataframe
  #' pathogen_cols: the names of the cols containing qPCR data
  #' type: "index" or "label" if "index" returns the index of the columns. If
  #'      "label", returns the column names
  #' returns: a vector of the indices or labels depending on the type specified.
  if(type == "index"){
    return(which(all_cols %in% setdiff(all_cols, pathogen_cols)))
  } else if(type == "labels"){
    return(all_cols[which(all_cols %in% setdiff(all_cols, pathogen_cols))])
  }else{
    stop("ERROR: Invalid type specified")
  }
}


get_max_stacks <- function(my_array){
  # I
  my_array_dims <- dim(my_array)
  for(x in 1:my_array_dims[3]){
    if(sum(is.na(my_array[,,x])) == my_array_dims[1] * my_array_dims[2]){
      return(x)
    }

  }

}


get_pathogenic_results <- function(){
  all_data <- read.csv("../data/pathogen_positive.csv")
  names(all_data) <- fix_names(names(all_data))
  # It's eliminating all values that are punctuation, need only to remove "."
  pathogens <- c("Salmonella", "Shigella", "Vibrio", "Yersinia", "Aeromonas",
                 "Campylobacter", "Plesiomonas", "Escherichia_coli",
                 "rotavirus", "norovirus", "adenovirus", "astrovirus",
                 "Cryptosporidium", "Giardia", "Entamoeba", "histolytica",
                 "Ascaris", "Trichuris", "Strongyloides", "Cyclospora", "Isospora",
                 "hookworm","EPEC", "ETEC", "EAEC", "tETEC", "tEPEC","EIEC", "EPEC")
  pathogenic_results <- isolate_pathogen_columns(
    all_data,pathogens,c(1:5,grep("Stool_type_gut_function|Sample_collection_date",
                                  names(all_data), ignore.case = T)))

  unique_binary <- generate_binary_num(pathogenic_results)

  pathogenic_results_encoded <- cbind(pathogenic_results, unique_binary)

  # names(pathogenic_results)[43] <- "Binary Combos"

  return(pathogenic_results_encoded)
}


get_pathogen_names <- function(m){
  #' To get the pagen names for the dataframe, m is dataset containing only the pathogen
  #' names.


}
#'
#' get_target_pathogens <- function(parent_dir = "../data/"){
#'   #' @title Get pathogen names.
#'   #'
#'   #' @description Get the target pathogen names as indicated by the Pathogen_List.xlsx file
#'   #' located in the copathogen analysis folder. Gets a dataframe with the pathogen names needed for both anaylsis.
#'   #'
#'   #' @details The dataframe contains five columns: One column contains variable names for the PROVIDE Diarrhea dataset, the second
#'   #' provides variable names for the asymptomatic dataset, and subsequent columns provide variable names for the MAL-ED dataset,
#'   #' general terms to be used across all studies, and finally column of 1's and 0's what indicate whether
#'   #' this variable is a pathogen or not. Each row is for a specific variable. This dataframe is used to match specific names in the
#'   #' respective dataset with the genreral term to be used across studies.
#'   #'
#'   #' @source toolkit.R
#'   #'
#'   #' @return dataframe containinng the names for asymptomatic provide data in one column, symptomatic in another,
#'   #'      and the mal-ed
#'
#'   # Import the original datafile.
#'   pathogen_list <- read.csv(paste(parent_dir, 'pathogen_list_flagged.csv', sep = ""),
#'                             # Import Dataset
#'                             header= T, stringsAsFactors = F) %>%
#'     # Include flagged Pathogens. See file.
#'     filter(Include > 0) %>%
#'     # These columns aren't needed
#'     select(-Comment, -Notes, -Include)
#'
#'   #Rename the columns
#'   names(pathogen_list) <- c("diar_provide", "asym_provide", "description", "mal_ed", "pathogen_class")
#'
#'   # Replace spaces with underscores, underscores work better in Rstudio
#'   flagged <- pathogen_list %>%
#'     mutate(mal_ed_fix = fix_names(mal_ed)) %>%
#'     # mal-ed is missing some pathogens, set those to NA. get rid of original column
#'     mutate(mal_ed_fix = ifelse(mal_ed_fix == "not_present", NA, mal_ed_fix)) %>%
#'     select(-mal_ed)
#'
#'   # Rename again.
#'   names(flagged) <- c( "asym_provide", "diar_provide", "actual", "pathogen_class", "mal_ed")
#'
#'   # For the asymptomatic dataset, give it the variable name for the symptomatic (diarrhea) dataset
#'   # and if the diarrhea vairbale name is blank, use the asymptomatic dataset variable name. This is
#'   # to syncronise the variable names across the PROVIDE dataset.
#'   flagged$asym_provide <- ifelse(flagged$diar_provide != "", flagged$diar_provide, flagged$asym_provide)
#'
#'   # convert to lower case
#'   final <- flagged %>% mutate_at(c("asym_provide", "diar_provide"), tolower)
#'
#'   # Change to participant to for all studies.
#'   final$actual[final$asym_provide == "sid"] <- "participant"
#'   # final$actual[final$mal_ed == find_names_matching(final$mal_ed, "Sample")] <- "sample"
#'
#'   # Variable to flag pathogen col names
#'   final$ispathogen <- 1
#'
#'   # Reset to zero for ID cols.
#'   final$ispathogen[final$asym_provide == "sid"] <- 0
#'   final$ispathogen[final$mal_ed == find_names_matching(final$mal_ed, "Sample")] <- 0
#'   final$ispathogen[final$actual == "age"] <- 0
#'
#'   # Add stool type col names. This has to be added seperately, since this variable is generated in the
#'   # scripts. The provide data is in seperate files, so I make the new variable when I bind the dataframes
#'   final <- rbind(final, c("stool_type", "stool_type", "stool_type", NA, "Stool_type_gut_function_EUPATH_0011745_", 0))
#'
#'   # missing in the original file. adding in here.
#'   final$diar_provide[final$asym_provide == "sid"] <- "sid"
#'   final$diar_provide[final$asym_provide == "agedays"] <- "ageatepi"
#'
#'   final <- final %>% transmute_all(tolower)
#'
#'   final$actual <- gsub("\\s+$", "", final$actual)
#'
#'
#'   return(final)
#' }

get_target_pathogens <- function(parent_dir = "../data/"){
  #' @title Get pathogen names.
  #'
  #' @description Get the target pathogen names as indicated by the Pathogen_List.xlsx file
  #' located in the copathogen analysis folder. Gets a dataframe with the pathogen names needed for both anaylsis.
  #'
  #' @details The dataframe contains five columns: One column contains variable names for the PROVIDE Diarrhea dataset, the second
  #' provides variable names for the asymptomatic dataset, and subsequent columns provide variable names for the MAL-ED dataset,
  #' general terms to be used across all studies, and finally column of 1's and 0's what indicate whether
  #' this variable is a pathogen or not. Each row is for a specific variable. This dataframe is used to match specific names in the
  #' respective dataset with the genreral term to be used across studies.
  #'
  #' @source toolkit.R
  #'
  #' @return dataframe containinng the names for asymptomatic provide data in one column, symptomatic in another,
  #'      and the mal-ed

  # Import the original datafile.
  pathogen_list <- read.csv(file.path(parent_dir, 'maled_provide_pathogen_name_key.csv'),
                            # Import Dataset
                            header= T, stringsAsFactors = F) %>%
    # Include flagged Pathogens. See file.
    filter(Include > 0) %>%
    # These columns aren't needed
    select(-Comment, -Notes, -Include)

  #Rename the columns
  names(pathogen_list) <- c("diar_provide", "asym_provide", "description", "mal_ed", "pathogen_class", "maled_60m")

  # Replace spaces with underscores, underscores work better in Rstudio
  flagged <- pathogen_list %>%
    mutate(mal_ed_fix = fix_names(mal_ed),
           mal_ed_60m_fix = fix_names(maled_60m)) %>%
    # mal-ed is missing some pathogens, set those to NA. get rid of original column
    mutate(mal_ed_fix = ifelse(mal_ed_fix == "not_present", NA, mal_ed_fix),
           mal_ed_60m_fix = ifelse(mal_ed_60m_fix == "not_present", NA, mal_ed_60m_fix)) %>%
    select(-mal_ed, -maled_60m)

  # Rename again.
  names(flagged) <- c( "asym_provide", "diar_provide", "actual", "pathogen_class", "mal_ed", "maled_60m")

  # For the asymptomatic dataset, give it the variable name for the symptomatic (diarrhea) dataset
  # and if the diarrhea vairbale name is blank, use the asymptomatic dataset variable name. This is
  # to syncronise the variable names across the PROVIDE dataset.
  flagged$asym_provide <- ifelse(flagged$diar_provide != "", flagged$diar_provide, flagged$asym_provide)

  # convert to lower case
  final <- flagged %>% mutate(across(c("asym_provide", "diar_provide"), ~ tolower(.x)))


  # Change to participant to for all studies.
  final$actual[final$asym_provide == "sid"] <- "participant"
  # final$actual[final$mal_ed == find_names_matching(final$mal_ed, "Sample")] <- "sample"

  # Variable to flag pathogen col names
  final$ispathogen <- 1

  # Reset to zero for ID cols.
  final$ispathogen[final$asym_provide == "sid"] <- 0
  final$ispathogen[final$mal_ed == find_names_matching(final$mal_ed, "Sample")] <- 0
  final$ispathogen[final$actual == "age"] <- 0

  # Add stool type col names. This has to be added seperately, since this variable is generated in the
  # scripts. The provide data is in seperate files, so I make the new variable when I bind the dataframes
  final <- rbind(final, c("stool_type", "stool_type", "stool_type", NA, "Stool_type_gut_function_EUPATH_0011745_", "Stool_type_EUPATH_0010869_", 0))

  # missing in the original file. adding in here.
  final$diar_provide[final$asym_provide == "sid"] <- "sid"
  final$diar_provide[final$asym_provide == "agedays"] <- "ageatepi"

  final <- final %>% transmute_all(tolower)

  final$actual <- gsub("\\s+$", "", final$actual)

  return(final)
}



isolate_lab_results <- function(df, pathogens, additional_columns){
  df_names <- names(df)
  return(df[df_names[c(additional_columns, grep(paste(pathogens, collapse = "|"), df_names, ignore.case = T))]])
}

isolate_pathogens <- function(df, pathogens, additional_columns){
  df_names <- names(df)
  return(df[df_names[c(additional_columns, grep(paste(pathogens, collapse = "|"), df_names, ignore.case = T))]])
}


isolate_pathogen_columns <- function(df, pathogens, id_columns){
  return(df %>% isolate_pathogens(.,pathogens, id_columns) %>% remove_cumulative_sums() %>%
           eliminate_empty_columns()%>% extract_binary_matrix() %>% cbind(.,df[id_columns]))
}

remove_cumulative_sums <- function(df){
  return(df[-grep("Cumulative_sum", names(df), ignore.case = T)])
}

remove_whitespace <- function(term){
  return(gsub(" ", "", term))
}

combine_observation <- function(x){
  if(sum(!is.na(x)) > 1){
    return("Collision")
  }else{
    return(!x[is.na(x)])
  }
}

make_qpcrsimple_tidy <- function(df_long, pathogen_names){
  #' Convertes the Dataframe of the simpleversion of the qPCR Data to a 'tidy' format. Allows for easier graphing.
  df_tidy <- df_long %>% tidyr::pivot_longer(c(pathogen_names), names_to = 'pathogen', values_to = 'is_present')
  return(df_tidy)
}


number_of_repeat_instances <- function(x){
  return(sum(!is.na(x)))
}

only_qpcr_data <- function(df, pathogen_cols){
  return(df[apply(!is.na(df[pathogen_cols]),1, sum) > 1,])
}


qpcr_continuous <- function(x){
  #' SImpler function to be wrapped in the 'transmute_all' function.
  #' Works by setting the value to 2 to the power of the cycle, minus
  #' one, since a power of zero is equal to 1.
  #'
  #' x: the column used by dplyr::transmte_all
  #' returns: the converted value for the column.
  return(ifelse(x > 0, 2^(35 - x) - 1, 0))
}


remove_empty_slices <- function(my_array){
  #' In order to make the indices more dynamic, the matrices are
  #' initially made to be really large. This function then removes
  #' the uneeded columns. First by finding the first slice with all
  #' the missing values, and removing all those indices.
  #'
  #' my_array: the three dimensional array with the extra slices.
  #' returns: the array with any slices that have all missing
  #'      values removed.
  max_val <- get_max_stacks(my_array)
  return(my_array[,,-c(max_val:dim(my_array)[3])])
}


remove_na_columns <- function(df){
  #' When generalizing the column names, there's columns returned with NA
  #' as the column name. These are not part of the 'desired' pathogens.
  #' Remove these columns from the dataset.
  #'
  #' df: the dataframe to reduce
  #' returns: the dataframe with all the NA columns removed.
  df[which(is.na(names(df)))] <- NULL
  return(df)
}


return_same <- function(x){
  #' This function is to save me typing later.
  return(x)
}


get_first_element <- function(x, idx){
  print(x[[idx]])
}


#' Write a function for getting the probability of diarrhea for a copathogen interaction in a network.
#' Get a table of the significant interactions and all the pathogen interactions.

p_of_diarrhea <- function(d_f){
  #' The probability of diarrhea for all the copathogen interaction pairs in the dataframe
  d_f
}

