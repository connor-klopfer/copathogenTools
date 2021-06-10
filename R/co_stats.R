
get_joint_probability <- function(d_f, names_key){
  #' @title Calculate the joint probability for two pathogens in the dataframe.
  #'
  #' @description Take the tidy representation of the dataset, meaning one
  #' observation per pathogen per study. find the probability of two
  #' pathogens showing up together out of all kids for their respective group.
  #'
  #' @param d_f: The dataframe as returned in the original get _both_datasets()
  #' function
  #' @param names_key: dataframe as the "names_key" containing all the names
  #' of the pathogens in question.
  #'
  #' @return A dataframe with pathogen 1, pathogen 2, and the observed joint
  #' probability of occuring together
  all_dfs <- list()
  iter = 1
  # For each level of the stool type, this is specific to segmenting by stool type
  for(stool in unique(d_f$stool_type)){
    # For each level of the study type.
    for(st in unique(d_f$study)){
      # Filter by study and stool type in the loop. Select the pathogen and convert
      # to a matrix.
      temp_l <- d_f %>% filter(study == st, stool_type == stool) %>% ungroup() %>%
        select(one_of(names_key$actual)) %>% convert_matrix_lists()
      # Extract the matrix
      temp_m <- temp_l[[1]]
      # Set NA's to zero, since we're taking the dot product.
      # NOTE: this might be an improper implementation of handling missing data,
      # I might have to change this.
      temp_m[is.na(temp_m)] <- 0
      # Get the co-occurance by multiplying the matrix
      temp_co <- t(temp_m) %*% temp_m
      # The diagonal is the number of times the pathogens appear in all kids

      # BUG: This is not returning the number of kids in the sutdy, this is returning
      # The total number of pathogens datected !!!! BUG FOUND: 23Mar20
      all_kids <- sum(psych::tr(temp_m %*% t(temp_m)))
      # FIX as of 23Mar20: All kids is set as the nmber of rows
      all_kids <- nrow(temp_m)

      # The diagonal should contain all the times that a pathogen shows up total in all the
      # kids.
      pathogen_count <- diag(temp_co)

      # Divide the number of times two pathogens appear by the porportion of kids that haee patogen
      # B by the total number of kids in the study.
      prob_b = matrix(rep(pathogen_count / all_kids, nrow(temp_co)), nrow = nrow(temp_co), byrow = T)

      print(paste("N rows (Should be:", 40, nrow(prob_b)))


      # For the indirect pobability, for any that are present in the

      # A - (A cap B)
      diagongal_repeated <- matrix(rep(pathogen_count, ncol(temp_co)), ncol = ncol(temp_co), byrow = F)

      print(paste("Dimension of dimensional repeated:", dim(diagongal_repeated)))

      conditional_negation <- (diagongal_repeated - temp_co)/matrix(rep(all_kids - pathogen_count, nrow(temp_co)), nrow = nrow(temp_co), byrow = T)






      # Take the counts of the number of times two pathogens show up and divide by the number of
      # kids total within the subgroup.
      temp_co <- temp_co / all_kids

      conditional_probability <- temp_co / prob_b
      # # Rename the column names and row names as the pathogen names from the original
      # # column names in the dataframe. Remember the CO matrix is symmetrical.
      # colnames(temp_co) <- temp_l[[2]]
      # row.names(temp_co) <- temp_l[[2]]

      # all_dfs[[iter]] <- temp_co %>%                        # Convert to dataframe
      #   as.data.frame() %>%                                 # Extract rownames as a variable
      #   mutate(first_path = row.names(.)) %>%
      #   tidyr::pivot_longer(-first_path,                    # Select everything but the first pathogen
      #                       names_to = 'second_path',       # variable and pivot to a longer format, with
      #                       values_to = 'joint_prob',       # each pathogen pair as a row.
      #                       values_drop_na = TRUE) %>%      # Drop NA values.
      #   mutate(study = st, stool_type = stool)              # Set the stool type and study

      all_dfs[[iter]] <- left_join(
        convert_matrix_to_dataframe(temp_co, "joint_prob", temp_l[[2]]),# Convert to dataframe
        convert_matrix_to_dataframe(conditional_probability, "cond_prob", temp_l[[2]]),
                                                by =c("first_path", "second_path")) %>%
        left_join(
          convert_matrix_to_dataframe(conditional_negation, 'negation', temp_l[[2]]),
          by = c("first_path", "second_path")) %>%
        mutate(study = st, stool_type = stool)



      iter = iter + 1
    }
  }

  # Append all the dataframes together
  all_jps <- do.call(rbind, all_dfs)

  # Set wor names to null, makes things look funny otherwise.
  row.names(all_jps) <- NULL
  return(all_jps)
}


convert_matrix_to_dataframe <- function(m, var_name, pathogen_names){
  colnames(m) <- pathogen_names
  row.names(m) <- pathogen_names

  final <- m %>%                        # Convert to dataframe
    as.data.frame() %>%                                 # Extract rownames as a variable
    mutate(first_path = row.names(.)) %>%
    tidyr::pivot_longer(-first_path,                    # Select everything but the first pathogen
                        names_to = 'second_path',       # variable and pivot to a longer format, with
                        values_to = var_name,       # each pathogen pair as a row.
                        values_drop_na = TRUE)

  return(final)
}


find_co_strength <- function(jps, marginals){
  #' Defference between observed joint probabililty and independent
  #'
  #' @description Finds the difference between the observed probability and the independent probability
  #' of two events occuring at the same time according to their marginal probabilites.
  #'
  #' TODO: CONTINUE HERE 23Mar20
  jp_mar <- left_join(left_join(
    jps,
    marginals %>% mutate(first_path = pathogen) %>%
      select(-pathogen), by = c("study", "stool_type", "first_path")) %>%
      mutate(path1_ind = marginal_p) %>%
      select(one_of(c("study", "stool_type", "first_path", "second_path", "path1_ind", "joint_prob"))),
    marginals %>% mutate(second_path = pathogen) %>%
      select(one_of(c("study", "stool_type", "second_path", "marginal_p"))),
    by = c("study", "stool_type", "second_path")) %>%
    mutate(path2_ind = marginal_p) %>% mutate(independence = path1_ind * path2_ind) %>%
    mutate(p_diff = log10((joint_prob)/independence))

  return(jp_mar)

}


negatvie_interaction <- function(d_f, path_names){
  # Of pathogens that are present in group A, which ones are not present in group B

  d_f %>% make_tidy(d_f, path_names)



}






