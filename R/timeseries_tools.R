#' The time series tools, for breaking the datset down into different time slices. 
#' 
add_messy_matrix <- function(m1, m2){
  #' For adding two matrices that might have missing values. Uses a series of 
  #' ifelse statements to add two matrices together that might have missing values. 
  #' 
  #' m1: a matrix that might have missing  values. 
  #' m2: a metrix that might have missing values. 
  #' 
  #' returns: the sum of the two matrices, where any missing values are dropped if 
  #'      the other matrix has a value. otherwise it's NA
  
  return(ifelse(is.na(m1), ifelse(is.na(m2), NA, m2), ifelse(is.na(m2), m1, m1+m2)))
}



bin_timeseries <- function(m, ages, name_array, n_bins, time_length = 372){
  #' Bin the 

  # Create a multidimensional vector, where each partition is a number of nothes, 
  days_partition <- bin_days(time_length, n_bins)
  
  new_m <- array(NA, dim = c(dim(m)[c(1:2)], n_bins))
  
  # for each time slice in the array 
  for(t_s in 1:dim(m)[3]){
    # for each bin in the timeseries. 
    for(bin in 1:n_bins){
      # if the age falls within the bin, then add it to the array. Otherwise, remove it. 
      # -> Problem: theres some NAs within the list. If two arrays are in the same bin, 
      # and one of them has an NA, take the non_NA value, if they're both NA, then the restult is NA, 
      temp <- m[,,t_s]
      temp[which(!(ages[,t_s] %in% days_partition[[bin]])), ] <- NA
      
      new_m[,,bin] <- add_messy_matrix(new_m[,,bin], temp)
      new_m[,,bin] <- ifelse(is.na(new_m[,,bin]), NA, ifelse(new_m[,,bin] > 0, 1, 0))
    }
  }
  return(new_m)
}

bin_days <- function(t_length, n_bins){
  days <- seq(t_length)
  
  return(split(days, ceiling(days/(t_length/n_bins))))
}


