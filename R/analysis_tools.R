#' Non-normality testing. using the kolmogornov smirnov test to analyze the normality of a distribution
#' of copathogen interactions. 


get_ks_df <- function(configuration_file, parent_dir = "../"){
  num_pairs <- ifelse("time_rep" %in% names(configuration_file), 
                                            dim(configuration_file)[2] - 1, 
                                            dim(configuration_file)[2])
  
  if("time_rep" %in% names(configuration_file)){
    pair_names <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
    ks_score <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
    ks_p <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
    sums <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
    aves <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
    time_point <- vector(mode = "numeric", length = num_pairs * length(unique(configuration_file$time_rep)))
  }else{
    pair_names <- vector(mode = "numeric", length = num_pairs)
    ks_score <- vector(mode = "numeric", length = num_pairs)
    ks_p <- vector(mode = "numeric", length = num_pairs)
    sums <- vector(mode = "numeric", length = num_pairs)
    aves <- vector(mode = "numeric", length = num_pairs)
    time_point <- vector(mode = "numeric", length = num_pairs)
  }
  
if('time_rep' %in% names(configuration_file)){
  for(num in 0:max(as.numeric(configuration_file$time_rep))){
    subset_df <- configuration_file %>% filter(time_rep == num)
      for(pair in 1:num_pairs){
          iterator <- (num * num_pairs) + pair
          pair_names[iterator] <- names(subset_df)[pair]
          ks_temp <- ks.test(x = subset_df[,pair], y = (rnorm(1000, mean = mean(subset_df[,pair]), sd = sd(subset_df[,pair]))))
          ks_score[iterator] <- ks_temp$statistic
          ks_p[iterator] <- ks_temp$p.value
          aves[iterator] <- mean(subset_df[,pair])
          sums[iterator] <- sum(subset_df[,pair])
          time_point[iterator] <- num
          
          ks_df <- data.frame(pair_names, ks_score, ks_p, sums, aves, time_point) %>% 
            tidyr::separate(pair_names, c('path1', 'path2'), "\\+") %>% 
            mutate(path1 = gsub("\\s{2}", "", path1), path2 = gsub("\\s{2}", "", path2)) %>% 
            filter(path1 != path2) %>%  # Remove items along the diagonal 
            mutate(combined = paste(path1, path2, sep = "+"), sig = ks_p < 0.05) 
      }
  }
  }else{
    for(pair in 1:num_pairs){
      pair_names[pair] <- names(config_reps)[pair]
      ks_temp <- ks.test(x = config_reps[,pair], y = (rnorm(1000, mean = mean(config_reps[,pair]), sd = sd(config_reps[,pair]))))
      ks_score[pair] <- ks_temp$statistic
      ks_p[pair] <- ks_temp$p.value
      aves[pair] <- mean(config_reps[,pair])
      sums[pair] <- sum(config_reps[,pair])
    
      ks_df <- data.frame(pair_names, ks_score, ks_p, sums, aves) %>% 
        tidyr::separate(pair_names, c('path1', 'path2'), "\\+") %>% 
        mutate(path1 = gsub("\\s{2}", "", path1), path2 = gsub("\\s{2}", "", path2)) %>% 
        filter(path1 != path2) %>%  # Remove items along the diagonal 
        mutate(combined = paste(path1, path2, sep = "+"), sig = ks_p < 0.05) 
    }
    }
  
  return(ks_df)
}


iterage_through_values <- function(ks_df, original_df){
  names(original_df) <- gsub("\\s{2}" , "", names(original_df))
  names(original_df) <- trimws(names(original_df), which = "right")
  ks_df$diar_p <- NA
  for(r in 1:nrow(ks_df)){
    ks_df$diar_p[r] <- apply_probability(ks_df, original_df, ks_df$path1[r], ks_df$path2[r], ks_df$time_point[r])
  }
  return(ks_df)
}


apply_probability <- function(ks_df, original_df, var1, var2, timepoint){
  var1 <- trimws(var1, which = "right")
  var2 <- trimws(var2, which = "right")
  if("age_bin" %in% names(original_df)){
    # print(2 %in% (original_df[[var1]] + original_df[[var2]]))
    if(2 %in% (original_df[[var1]] + original_df[[var2]])){
      # print(names(ks_df))
      # print(paste("TimePoint:", timepoint))
      # View(ks_df)
      
      # # print(sum(original_df[[var1]] == 1))
      # print(sum(original_df[[var2]] == 1))
      # print(sum(original_df[['age_bin']] == as.numeric(timepoint)))
    final <- original_df %>% select(var1, var2, "stool_type", "age_bin") %>%
      filter(.[[var1]] == 1, .[[var2]] == 1, .[['age_bin']] == as.numeric(timepoint)) %>% 
      summarise(counts = sum(stool_type == 1)/n()) %>% as.numeric()
    }else{
      final <- 0
    }
    
  }else{
    final <- original_df %>% select(var1, var2, "stool_type") %>%
      filter(.[[var1]] == 1, .[[var2]] == 1) %>%
      summarise(counts = sum(stool_type == "Diarrhea")/n()) %>% as.numeric()
  }
  # print(final)
  return(final)
}




bin_timeseries_dplyr <- function(x, n_bins){
  MAX_VALUE = 365 + 7
  bin_size = MAX_VALUE / n_bins
  
  final = floor(x/bin_size)
  
  return(final)
}