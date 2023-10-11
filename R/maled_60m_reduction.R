#' This is for reducing the original dataset to only positive values
#' Use this for reducing the original dataset, any other simplified method seems
#' to crash R
#'
#' Update: 04Sep20: This was updated to simpleify the dataset from the MAL-ED 60m update


reduced_full_maled_dataset <- function(parent_dir){
  samples <- data.table::fread(file.path(parent_dir, "ISASimple_Gates_MAL-ED_phase3_RSRC_samples.txt"))

  observations <- data.table::fread(file.path(parent_dir, "ISASimple_Gates_MAL-ED_phase3_RSRC_observations.txt"))

  households <- data.table::fread(file.path(parent_dir, "ISASimple_Gates_MAL-ED_phase3_RSRC_households.txt"))

  # Append the observation ID, the date of the observation, and the age of the participant at the time
  # of the observation.

  # Get the names of the columns containing qPCR data.
  qpcr_columns <- find_intersection_of_terms(samples, c(" Ct ", "value"))

  # We want observations where there is qPCR data, make a column to act as a flag.
  samples$qpcr_tested <- samples %>%
    select(all_of(qpcr_columns)) %>%
    rowSums(!is.na(.))

  samples_reduced <- samples %>% filter(qpcr_tested > 0)

  samples_w_dates <- left_join(
    samples_reduced,
    observations %>% select(
      Participant_Id,
      Observation_Id,
      `Observation date [EUPATH_0004991]`,
      `Age (days) [EUPATH_0000579]`),
    by = c("Participant_Id", "Observation_Id")) %>%
    left_join(households %>% # Append the country
                select(Household_Observation_Id, `Country [OBI_0001627]`) %>%
                rename(Household_Id = Household_Observation_Id), by = "Household_Id")

  ready_to_write <- samples_w_dates %>%
    select(
      Participant_Id,
      Observation_Id,
      `Stool type [EUPATH_0010869]`,
      `Observation date [EUPATH_0004991]`,
      `Age (days) [EUPATH_0000579]`,
      `Country [OBI_0001627]`,
      all_of(qpcr_columns)
    )


  write.csv(ready_to_write, file.path(parent_dir, "maled_samples_60m_reduced.csv"))

  rm(samples)
  rm(observations)
  rm(households)

}


