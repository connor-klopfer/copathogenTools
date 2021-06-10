#### For cleaning the dataset according to the algorithm outlined by Ross and Dorothy November
#### Started: 11Nov19
#### By: Connor Klopfer


#' maled_name_dictionary <- function(){
#'   #' @title Convert the e.coli variables to their general terms.
#'   #'
#'   #' @description Generate a dictionary with the general terms as the keys and
#'   #' the sepcific variable names in MAL-ED as the values.
#'   #'
#'   #'
#'   #' @return : named list (dictionary) with the abbreviations used in the code
#'   #'      and the actual names as the values.
#'
#'   # These are the specific variables in MAL-ED
#'   pathogen_labels <- c("EAEC_aaiC_Ct_value_by_TAC_result_EUPATH_0015346_",
#'                        "EAEC_aatA_Ct_value_by_TAC_result_EUPATH_0015348_",
#'                        "EPEC_bfpA_Ct_value_by_TAC_result_EUPATH_0015353_",
#'                        "EPEC_eae_Ct_value_by_TAC_result_EUPATH_0011857_",
#'                        "STEC_stx1_Ct_value_by_TAC_result_EUPATH_0015384_",
#'                        "STEC_stx2_Ct_value_by_TAC_result_EUPATH_0015385_",
#'                        "ETEC_LT_pos_Ct_value_by_TAC_result_EUPATH_0015355_",
#'                        "ETEC_STh_Ct_value_by_TAC_result_EUPATH_0015356_",
#'                        "ETEC_STp_Ct_value_by_TAC_result_EUPATH_0015357_"
#'   )
#'
#'   # The general names
#'   names(pathogen_labels) <- c("aaiC", "aatA", "bfpa", "eae", "stx1", "stx2",
#'                               "eteclt", "sth", "stp")
#'   return(pathogen_labels)
#' }


maled_name_dictionary <- function(vers = "60m"){
  #' @title Convert the e.coli variables to their general terms.
  #'
  #' @description Generate a dictionary with the general terms as the keys and
  #' the sepcific variable names in MAL-ED as the values.
  #'
  #'
  #' @return : named list (dictionary) with the abbreviations used in the code
  #'      and the actual names as the values.

  if(vers == "60m"){
    pathogen_labels <- c("EAEC_aaiC_Ct_value_by_TAC_result_EUPATH_0015346_",
                         "EAEC_aatA_Ct_value_by_TAC_result_EUPATH_0015348_",
                         "EPEC_bfpA_Ct_value_by_TAC_result_EUPATH_0015353_",
                         "EPEC_eae_Ct_value_by_TAC_result_EUPATH_0011857_",
                         "STEC_stx1_Ct_value_by_TAC_result_EUPATH_0015384_",
                         "STEC_stx2_Ct_value_by_TAC_result_EUPATH_0015385_",
                         "ETEC_LT_pos_Ct_value_by_TAC_result_EUPATH_0015355_",
                         "ETEC_STh_Ct_value_by_TAC_result_EUPATH_0015356_",
                         "ETEC_STp_Ct_value_by_TAC_result_EUPATH_0015357_"
    )
  }else if(vers == "24m"){
    pathogen_labels <- c("EAEC_aaiC_Ct_value_by_TAC_result_EUPATH_0015346_",
                         "EAEC_aatA_Ct_value_by_TAC_result_EUPATH_0015348_",
                         "EPEC_bfpA_Ct_value_by_TAC_result_EUPATH_0015353_",
                         "EPEC_eae_Ct_value_by_TAC_result_EUPATH_0011857_",
                         "STEC_stx1_Ct_value_by_TAC_result_EUPATH_0015384_",
                         "STEC_stx2_Ct_value_by_TAC_result_EUPATH_0015385_",
                         "ETEC_LT_pos_Ct_value_by_TAC_result_EUPATH_0015355_",
                         "ETEC_STh_Ct_value_by_TAC_result_EUPATH_0015356_",
                         "ETEC_STp_Ct_value_by_TAC_result_EUPATH_0015357_"
    )

  }

  # These are the specific variables in MAL-ED

  # The general names
  names(pathogen_labels) <- c("aaiC", "aatA", "bfpa", "eae", "stx1", "stx2",
                              "eteclt", "sth", "stp")
  return(pathogen_labels)
}


#'
#' harmonize_pathogen_data <- function(d_f, df_labels){
#'   #' @title Hamonize the E.coli definitions for MAL-ED and PROVIDE datasets
#'   #'
#'   #' @param d_f : The MAL-ED dataset pulled from the \code{maled_full_reduced.csv} file.
#'   #'
#'   #' @param df_labels : a dictionary of variable names for the e.coli variables as returned
#'   #' by \code{\link{maled_name_dictionary}}
#'   #'
#'   #' @description Function for harmonizing the pathogen data with the PROVIDE study. Takes in the MAL-ED data
#'   #' frame after cleaning and adds the adjustements to match.
#'   #'
#'   #' @return The dataframe with the same rules applied to e.coli qPCR variables between MAL-ED and PROVIDE
#'   #'
#'   #' @details Applies the following functions to the MAL-ED dataset to hamronize with PROVIDE
#'   #' \itemize{
#'   #' \item \code{\link{eaec_definition}}
#'   #' \item \code{\link{tepec_definition}}
#'   #' \item \code{\link{atepec_definition}}
#'   #' \item \code{\link{stec_definition}}
#'   #' \item \code{\link{lt_etec_definition}}
#'   #' \item \code{\link{st_etec_definition}}
#'   #' \item \code{\link{etec_definition}}
#'   #' }
#'   #'
#'
#'   # Setting the threshold as 35. vector() initialises as 0's
#'   threshold_vector <- vector(mode = 'numeric', length = nrow(d_f)) + 35
#'
#'   # The reason for using the dictionary is for this section. Since the general terms are the
#'   # keys and the original variable name is the value. When selecting the term from df_labels
#'   # thats general, its returning the specific name, and indexes the dataframe by that specific
#'   # name. This method allows for easier readibility instead of writing the specific variable
#'   # names in this section of the code.
#'
#'   new_df <- d_f %>% dplyr::mutate(
#'     EAEC_ct_value_custom = eaec_definition(.[[df_labels[["aaiC"]]]],
#'                                      .[[df_labels[["aatA"]]]],threshold_vector),
#'     eEPEC_ct_value_custom = tepec_definition(.[[df_labels[["bfpa"]]]],
#'                                        .[[df_labels[["eae"]]]], threshold_vector),
#'     atepec_ct_value_custom = atepec_definition(.[[df_labels[["bfpa"]]]],
#'                                          .[[df_labels[["eae"]]]],
#'                                          .[[df_labels[["stx1"]]]],
#'                                          .[[df_labels[["stx2"]]]],threshold_vector),# Variables needed
#'     tepec_ct_value_custom = tepec_definition(.[[df_labels[["bfpa"]]]],
#'                                        .[[df_labels[["eae"]]]], threshold_vector),
#'     stec_ct_value_custom = stec_definition(.[[df_labels[["bfpa"]]]],
#'                                      .[[df_labels[["eae"]]]],
#'                                      .[[df_labels[["stx1"]]]],
#'                                      .[[df_labels[["stx2"]]]],threshold_vector),
#'     lt_etec_ct_value_custom = lt_etec_definition(.[[df_labels[["eteclt"]]]],
#'                                            .[[df_labels[["sth"]]]],
#'                                            .[[df_labels[['stp']]]], threshold_vector),
#'     st_etec_ct_value_custom = st_etec_definition(.[[df_labels[["eteclt"]]]],
#'                                            .[[df_labels[["sth"]]]], threshold_vector),
#'     etec_ct_value_custom = etec_definition(.[[df_labels[["eteclt"]]]],
#'                                      .[[df_labels[["sth"]]]],
#'                                      .[[df_labels[['stp']]]], threshold_vector)
#'     )
#'
#'   return(new_df)
#' }


harmonize_pathogen_data <- function(d_f, df_labels){
  #' @title Hamonize the E.coli definitions for MAL-ED and PROVIDE datasets
  #'
  #' @param d_f : The MAL-ED dataset pulled from the \code{maled_full_reduced.csv} file.
  #'
  #' @param df_labels : a dictionary of variable names for the e.coli variables as returned
  #' by \code{\link{maled_name_dictionary}}
  #'
  #' @description Function for harmonizing the pathogen data with the PROVIDE study. Takes in the MAL-ED data
  #' frame after cleaning and adds the adjustements to match.
  #'
  #' @return The dataframe with the same rules applied to e.coli qPCR variables between MAL-ED and PROVIDE
  #'
  #' @details Applies the following functions to the MAL-ED dataset to hamronize with PROVIDE
  #' \itemize{
  #' \item \code{\link{eaec_definition}}
  #' \item \code{\link{tepec_definition}}
  #' \item \code{\link{atepec_definition}}
  #' \item \code{\link{stec_definition}}
  #' \item \code{\link{lt_etec_definition}}
  #' \item \code{\link{st_etec_definition}}
  #' \item \code{\link{etec_definition}}
  #' }
  #'

  # Setting the threshold as 35. vector() initialises as 0's
  threshold_vector <- vector(mode = 'numeric', length = nrow(d_f)) + 35

  # The reason for using the dictionary is for this section. Since the general terms are the
  # keys and the original variable name is the value. When selecting the term from df_labels
  # thats general, its returning the specific name, and indexes the dataframe by that specific
  # name. This method allows for easier readibility instead of writing the specific variable
  # names in this section of the code.

  new_df <- d_f %>% dplyr::mutate(
    EAEC_ct_value_custom = eaec_definition(.[[df_labels[["aaiC"]]]],
                                           .[[df_labels[["aatA"]]]],threshold_vector),
    eEPEC_ct_value_custom = tepec_definition(.[[df_labels[["bfpa"]]]],
                                             .[[df_labels[["eae"]]]], threshold_vector),
    atepec_ct_value_custom = atepec_definition(.[[df_labels[["bfpa"]]]],
                                               .[[df_labels[["eae"]]]],
                                               .[[df_labels[["stx1"]]]],
                                               .[[df_labels[["stx2"]]]],threshold_vector),# Variables needed
    tepec_ct_value_custom = tepec_definition(.[[df_labels[["bfpa"]]]],
                                             .[[df_labels[["eae"]]]], threshold_vector),
    stec_ct_value_custom = stec_definition(
                                           .[[df_labels[["eae"]]]],
                                           .[[df_labels[["bfpa"]]]],
                                           .[[df_labels[["stx1"]]]],
                                           .[[df_labels[["stx2"]]]],threshold_vector),
    lt_etec_ct_value_custom = lt_etec_definition(.[[df_labels[["eteclt"]]]],
                                                 .[[df_labels[["sth"]]]],
                                                 .[[df_labels[['stp']]]], threshold_vector),
    st_etec_ct_value_custom = st_etec_definition(.[[df_labels[["eteclt"]]]],
                                                 .[[df_labels[["sth"]]]], threshold_vector),
    etec_ct_value_custom = etec_definition(.[[df_labels[["eteclt"]]]],
                                           .[[df_labels[["sth"]]]],
                                           .[[df_labels[['stp']]]], threshold_vector)
  )

  return(new_df)
}



eaec_definition <- function(aaic, aata, threshold){
  #' @title The specific rule set to apply to EAEC
  #'
  #' @description Based on the results of AAIC and AATA on the TAQ Card, implements
  #' the rules described in 'Ecoli defintions.docx' file.
  #'
  #' @details The rules applied:
  #' \enumerate{
  #' \item If AAIC and AATA is NA: return NA
  #' \item Set EAEC as the lowest or non-missing qPCR value.
  #' \item Negative (Ct = 35) if both AAIC and AATA are both missing.}
  #'
  #' @return a vector of numeric values in the range of [0,35] with
  #' missing values to serve as the column for eaec

  #Set all to negative.
  eaec <- vector(mode = 'numeric', length = length(aaic)) + 35


  # If AAIC and AATA is missing, EAEC is missing.
  eaec[which(is.na(aaic) & is.na(aata))] <- NA


  # REMOVED: Not needed: 23Jun20
    # # If AAIC is missing, set EAEC as the maximum of the threshold and the
  # # value of AATA at that row. This ensure that wherever AAIC is missing,
  # # EAEC is preset to the negative value.
  # eaec[which(is.na(aaic))] <- pmax(aata[which(is.na(aaic))],
  #                                  threshold[which(is.na(aaic))])
  #
  # # Same for AATA
  # eaec[which(is.na(aata))] <- pmax(aaic[which(is.na(aata))],
  #                                  threshold[which(is.na(aata))])

  # If either are below the treshold, set as that miniumum value. Otherwise,
  # leave it, since it's already set at 35
  eaec[which(aaic < threshold | aata < threshold)] <- pmin(aaic[which(aaic < threshold | aata < threshold)],
                                                           aata[which(aaic < threshold | aata < threshold)], na.rm = T)

  return(eaec)
}


# TODO: CONFIRM THAT EAE IS THE SAME AS EPECEA
tepec_definition <- function(bfpa, eae, threshold){
  #' @title tEPEC defintion
  #'
  #' @description Setting the values of TEPEC based ont he qPCR results
  #' of bfpA and etecea
  #'
  #' @details uses the following rules: \itemize{
  #' \item Use bfpA is etecea is negative or missing.
  #' \item Use minuimum of bfpA and etecea if both are present.
  #' \item Missing if bfpA is missing
  #' \item neagtive if bfpA is negative.}
  #'
  #' @return A numeric vector of values in range [0,35] with missing values
  #' for tEPEC

  # Set as the negative value (35) for the vector.
  tepec <- vector(mode = "numeric", length = length(bfpa)) + 35

  # If BFPA is negative set as NA
  tepec[is.na(bfpa)] <- NA

  # Use minimum of bfpa and eae if bfpa is positive, otherwise will stay as the
  # preset 35 or NA for missing values.
  tepec[which(bfpa < threshold)] <- pmin(bfpa[which(bfpa < threshold)] ,
                                         eae[which(bfpa < threshold)] , na.rm = T)

  return(tepec)
}

# TODO: What is eae
# TODO: The problem here is that theres an eae, bfpa, and an atypical marker
# EPEC, but this is not factored into Dorothy's instructions. I need to figure
# out what to do with this when I get my hands on the codebook.
# Revision: I think there's some redundancy here, lets test to see if this definition
# lines up with the MaL-ED vairables for the respective pathogens
atepec_definition <- function(bfpa, eae, stx1, stx2, threshold){
  #' @title Set the rules for atEPEC
  #'
  #' @description Applying the rules provided in 'Ecoli definitions.docx'
  #' for atEPEC and setting as a new column.
  #'
  #' @details Rules applied are as follows: \itemize{
  #' \item if epecea positive and bfpa, stx1, and stx2 are negative then positive
  #' \item missing if bfpa, stx1, or stx2 is missing
  #' \item negative otherwise}
  #'
  #' @return A vector of values in range [0,35] with missing values to represent the
  #' column foe atEPEC

  # Set all values as negative.
  atepec <- vector(mode = "numeric", length = length(bfpa)) + 35

  # This is a long one, so I'm just going to save the boolean index and call it
  # If epecea is positive and bfpa is negative and stx1 and stx2 is negative
  eae_index <- which(eae < threshold & bfpa >= threshold &
                  stx1 >= threshold & stx2 >= threshold)

  # for all values where the index is TRUE, set as the eae value.
  atepec[eae_index] <- eae[eae_index]

  # If any of these values are missing, set as missing.
  atepec[which(is.na(eae) | is.na(bfpa) | is.na(stx1) | is.na(stx2))] <- NA


  return(atepec)
  #
  # if(eae < threshold & bpfa >= threshold & stx1 >= threshold & stx2 >= threshold){
  #   return(eae)
  # }else if (is.na(eae) | is.na(bpfa) | is.na(stx1) | is.na(stx2)){
  #   return(NA)
  # }else{
  #   return(threshold)
  # }
}



stec_definition <- function(eae, bfpa, stx1, stx2, threshold){
  #' @title Ruleset for STEC
  #'
  #' @description Apply the ruleset in 'Ecoli definitions.docx' in
  #' the copathogen folder. Applies the ruleset to MAL-ED to so we can
  #' merge with PROVIDE
  #'
  #' @param eae : column containing qPCR results measuring EPECEA
  #' @param bpfa : column containing qPCR results measuring bfpA
  #' @param stx1 : column containing qPCR results measuring STX1
  #' @param stx2 : column containing qPCR results measuring STX2
  #' @param threshold : a vector with length equal to the number
  #' of rows in another column containing the value for a negative
  #' result (35)
  #'
  #' @details applies the following rules: \itemize{
  #' \item positive if epecea is positive and bfpA is negative and stx1 and stx2 is
  #' negative \itemize{
  #' \item If condition met, use the lowest Ct of stx1 and stx2, and the highest Cq
  #' value between that number and epecea}
  #' \item Negative if epecea negative or bfpA positive and stx1 or stx2 non-missing
  #' \item missing if epecea or bfpA missing or stx1 and stx2 missing
  #' }
  #'
  #' @return a numeric vector in range [0,35] to be used for stec_definition

  # Set the values as negative
  stec <- vector(mode = "numeric", length = length(bfpa)) + 35

  # Create the first boolean vector, if epecea positive. bfpa
  # negative, and at least stx1 or stx2 poisitive
  condition1 <- which(eae < threshold & bfpa >= threshold &
                        (stx1 < threshold | stx2 < threshold))

  # If condition 1 is met, use the lowest value for stx and the hgihest between lowest
  # stx and eae
  stec[condition1] <- pmax(eae[condition1], pmin(stx1[condition1], stx2[condition1]))

  # Unsure if this is needed but just to be sure in case the two are not exhaustive
  # This is just ensure one of the rules are followed.
  condition2 <- which((eae >= threshold | bfpa < threshold) & (!is.na(stx1) | !is.na(stx2)))

  # If condition 2 is met, then set to negative value
  stec[condition2] <- threshold[condition2]

  # If any are missing, set as missing.
  stec[which(is.na(eae) | is.na(bfpa) | (is.na(stx1) & is.na(stx2)))] <- NA

  return(stec)
}

lt_etec_definition <- function(eteclt, STh, STp, threshold){
  #' @title Rule set for LT-ETEC
  #'
  #' @description Implement rule set used in PROVIDE and appply to the
  #' MAL-ED dataset according to the rules provided in 'Ecoli definitions.docx'
  #'
  #'@param eteclt : column containing the qPCR results for eteclt
  #'@param STh : column containing qpCR results for STh
  #'@param STp : column containing qPCR results for STp
  #'
  #'@details the rules applied are: \itemize{
  #'\item positive if eteclt positive and STh, STp both negative (use Cq of eteclt)
  #'\item Negative if Eteclt negative
  #'\item Negative if Eteclt positive and and with STh or STp positive
  #'\item Negative if Eteclt missing and STh or STp positive
  #'\item missing if Eteclt, STh, and STp missing}
  #'
  #'@return A nummeric vector contianing the results for LT ETEC based on these
  #'applied rules.

  # Set the negative values
  lt_etec <- vector(mode = "numeric", length = length(eteclt)) + 35

  # IF eteclt positive and STh and STp negative.
  condition1 <- which(eteclt < threshold & STh >= threshold & STp >= threshold)

  # IF the condition is met, set as eteclt, since the value is < 35 anyways as per
  # the above condition.
  lt_etec[condition1] <- eteclt[condition1]

  # If all are missing, set as missing
  lt_etec[which(is.na(eteclt) & is.na(STh) & is.na(STp))] <- NA

  # IF eteclt negative or STh or STp positive, set as negative
  condition2 <- which(eteclt >= threshold | STh < threshold | STp < threshold)

  lt_etec[condition2] <- threshold[condition2]

  return(lt_etec)
}

st_etec_definition <- function(eteclt, STh, threshold){
  #' @title Apply ST-ETEC ruleset
  #'
  #' @description Apply the ruleset for ST-ETEC  to the MAL-ED
  #' dataset as defined in in 'Ecoli definitions.docx' file.
  #'
  #' @param eteclt : column containing qPCR results for eteclt
  #' @param STh : column containing qPCR results for STh
  #'
  #' @details Applies the following rules: \itemize{
  #' \item Positive if STh positive, use lower Ct
  #' \item Positive if STh positive, eteclt negative of missing, use STh Cq
  #' \item Negatvie if STh negative (etec positive, negative or missing)
  #' \item Missing if STh missing}
  #'
  #' @return Numeric column in range [0,35] with missing values to be used
  #' as the ST-ETEC column
  #'

  # Set the negative values
  st_etec <- vector(mode = "numeric", length = length(eteclt)) + 35

  # if STEC negative, set as negative
  st_etec[which(STh >= threshold)] <- 35

  # if STh missing, then set as missing
  st_etec[which(is.na(STh))] <- NA

  # Where the STh and etecltis negative
  condition1 <- which(STh < threshold & eteclt < threshold)

  # Set as the minium between the two
  st_etec[condition1] <- pmin(STh[condition1], eteclt[condition1])

  # IF STh is positive and eteclt missing or negative, use STh. I'm pretty
  # sure that the above condition will cover this, but I included this to be explicit
  condition2 <- which(STh < threshold & (eteclt >= threshold | is.na(eteclt)))

  st_etec[condition2] <- STh[condition2]

  return(st_etec)
}

# TODO: ETEC
etec_definition <- function(eteclt, STh, STp, threshold){
  #' @title Apply ETEC ruleset to MAL-ED
  #'
  #' @description Apply the ruleset the MAL-ED dataset using STh
  #' STp and eteclt markers to synchroise the MAL-ED with PROVIDE
  #'
  #' @param eteclt : column containing qPCR results for eteclt
  #' @param STh : column containing qPCR results for STh
  #' @param STp : column containing qPCR results for STp
  #'
  #' @details the rulesets applied are: \itemize{
  #' \item positive if any are positive, use lowest Ct
  #' \item negative if all are negative or missing,
  #' \item missing if all are missing
  #' }
  #'
  #' @return numeric vector in range [0,35] with missing values to
  #' serve as the column for ETEC

  # Initialize missing values
  etec <- vector(mode = "numeric", length = length(eteclt)) + 35

  # Get the minimum of the three, excluding missing values
  etec <- pmin(STh, STp, eteclt, na.rm = T)

  # If they're all missing, then set as a missing value
  etec[which(is.na(STh) & is.na(STp) & is.na(eteclt))] <- NA

  return(etec)
}


