#' Post-processing to work with results from the configuration model.



find_shared_sig_pairs <- function(results){
  #' @title Find the shared significant pairs between the results in both datasets
  #'
  #' @description Looks in both subsets of the data and find the pathogen pairs that
  #' are significant in both studies for both datasets. NOTE: C_score impolementation
  #' is in the code, but deactivated because I haven't been using it.
  #'
  #' @param results the results of the configuration model for both studies and both
  #' stool types, condesned into one dataframe.


  new_var_names <- c("path1", "path2", "provide", "maled")

  all_intersections <- rbind(

    inner_join(
      results %>% filter(study == "provide", stool_type == "Diarrhea", Flag > 0),
      results %>% filter(study == "maled", stool_type == "Diarrhea", Flag > 0),
      by = c("path1", "path2"),
      suffix = c("_provide", "_maled")
    ) %>%
      select(path1, path2, contains("actual"), -contains("c_actual")) %>%
      rename_with(.fn = ~ new_var_names, .cols = everything()) %>%
      mutate(stool_type = "Diarrhea",
             # int_type = "count"
             ),

    inner_join(
      results %>% filter(study == "provide", stool_type == "Asymptomatic", Flag > 0),
      results %>% filter(study == "maled", stool_type == "Asymptomatic", Flag > 0),
      by = c("path1", "path2"),
      suffix = c("_provide", "_maled")) %>%
      select(path1, path2, contains("actual"), -contains("c_actual")) %>%
      rename_with(.fn = ~ new_var_names, .cols = everything()) %>%
      mutate(stool_type = "Asymptomatic",
             # int_type = "count"
             )

    # inner_join(
    #   results %>% filter(study == "provide", stool_type == "Diarrhea", c_Flag > 0),
    #   results %>% filter(study == "maled", stool_type == "Diarrhea", c_Flag > 0),
    #   by = c("path1", "path2"),
    #   suffix = c("_provide", "_maled")
    # ) %>%
    #   select(path1, path2, contains("c_actual")) %>%
    #   rename_with(.fn = ~ new_var_names, .cols = everything()) %>%
    #   mutate(stool_type = "Diarrhea", int_type = "c_score"),


    # inner_join(
    #   results %>% filter(study == "provide", stool_type == "Asymptomatic", c_Flag > 0),
    #   results %>% filter(study == "maled", stool_type == "Asymptomatic", c_Flag > 0),
    #   by = c("path1", "path2"),
    #   suffix = c("_provide", "_maled")
    # ) %>%
    #   select(path1, path2, contains("c_actual")) %>%
    #   rename_with(.fn = ~ new_var_names, .cols = everything()) %>%
    #   mutate(stool_type = "Asymptomatic", int_type = "c_score")
  ) %>%
    tidyr::pivot_longer(c("provide", "maled"), names_to = "study")


  # NOTE: ADDED by CK on 15Dec2020

  final <-  inner_join(
    all_intersections %>%
      # filter(int_type == "c_score") %>%
      select(-value),
    results,
    by = c("path1", "path2", "stool_type", "study"))

  # return(all_intersections)
  return(final)
}


get_c_actual <- function(d_f){
  d_f$c_actual <- NA
  for(r in 1:nrow(d_f)){
    p1 <- d_f$actual[which(d_f$path1 == d_f$path1[r] & d_f$path2 == d_f$path1[r])]
    p2 <- d_f$actual[which(d_f$path1 == d_f$path2[r] & d_f$path2 == d_f$path2[r])]

    d_f$c_actual[r] <- ((p1 - d_f$actual[r]) * (p2 - d_f$actual[r])) / (p1 * p2)

  }

  return(d_f)
}


get_all_maled_tables <- function(parent_dir = ".."){
  #' @title All maled subset tables
  #'
  #' @description Get all of the tables from the MAL-ED study. COmpile them into a list to return. The tables
  #' will include Observations, saamples, participants, houlsholds, and ontology. This allows for easier data import.
  #'
  #' @param parent_dir: The parent directory containing the dataset text files from the MAL-ED website.
  #'
  #'  @return A list of all the datasets with the keys "ontology", "observations", "samples", "participants", and "housholds"

  all_dfs <- list()
  all_dfs[["ontology"]] <- data.table::fread(
    paste(parent_dir, "data/phase3/ISASimple_Gates_MAL-ED_phase3_RSRC_ontologyMetadata.txt", sep = "/"))

  all_dfs[["samples"]] <- data.table::fread(
    paste(parent_dir, "data/phase3/ISASimple_Gates_MAL-ED_phase3_RSRC_samples.txt", sep = "/"))

  all_dfs[["observations"]] <- data.table::fread(
    paste(parent_dir, "data/phase3/ISASimple_Gates_MAL-ED_phase3_RSRC_observations.txt", sep = "/"))

  all_dfs[["participants"]] <- data.table::fread(
    paste(parent_dir, "data/phase3/ISASimple_Gates_MAL-ED_phase3_RSRC_participant.txt", sep = "/"))

  all_dfs[["households"]] <- data.table::fread(
    paste(parent_dir, "data/phase3/ISASimple_Gates_MAL-ED_phase3_RSRC_households.txt", sep = "/"))

  return(all_dfs)
}


filter_out_symmetric <- function(d_f){
  #' @title Flip symmetric pathogens
  #'
  #' @description Faster, simpler, and safer version of removing pathogens that
  #' appear on eaiter side of the diagonal.
  #'
  #' @param d_f : Dataframe containg pathogen pairs that have pairs on both sides of the
  #' diagonal of a symmetric matrix.
  #'
  #' @return dataframe with redundant pathogen pairs removed.
  #'
  final <- d_f %>%
    mutate(path1_switch = ifelse(path1 > path2, path1, path2),
           path2_switch = ifelse(path1 > path2, path2, path1)) %>%
    select(-path1, -path2) %>%
    rename(path1 = path1_switch, path2 = path2_switch) %>%
    select(path1, path2, everything()) %>%
    filter(path1 != path2) %>%
    distinct()
  return(final)
}

combine_like_pathogens <- function(d_f){
  #' @title Combine Ecoli and EPEC pathogens
  #'
  #' @description Some pathogens are similar, and after talking with domain
  #' experts, we found Typical EPEC should be named as epec, atypical is dropped,
  #' and all the ETEC variants are combined.
  #'
  #' @return Dataframe with the combined like pathogens.
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
