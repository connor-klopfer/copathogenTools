#' For building a plotly network from the CI tables.  Needs to be written exactly as in the the edge_rewiring.Rmd programs
#' Apadted from the plotly for R tutorials.


build_network <- function(edge_df, my_title = ""){


  # my_edge_list <- practice_ci %>% select(one_of(c("Path1", "Path2")))
  G <- graph_from_edgelist(as.matrix(edge_df))
  # L <- layout.circle(G)
  L <- layout.fruchterman.reingold(G)



  verts <- V(G)
  edges_from_layout <- as.data.frame(as_edgelist(G, names = F))


  Xn <- L[,1]
  Yn <- L[,2]


  my_network <- plot_ly(x = Xn,
                     y = Yn,
                     mode = "markers",
                     text =  verts$name,
                     hoverinfo = "text",
                     type = "scatter")

  edge_lines <- build_edges(edges_from_layout, Xn, Yn)

  my_axis <- list(title = my_title,
                  showgrid = F,
                  showticklabels = F,
                  zeroline = F)

  final <- layout(my_network,
                  xaxis = my_axis,
                  yaxis = my_axis,
                  shapes = edge_lines)

  return(final)
}


combine_study_results <- function(df1, df2, study1, study2){
  df1$study <- study1
  df2$study <- study2
  return(rbind(df1, df2))
}


combine_study_results_4 <- function(study1_sym, study2_sym, study1_asym, study2_asym, study_var1, study_var2){
  study1_sym$study <- study_var1
  study2_sym$study <- study_var2

  sym <- rbind(study1_sym, study2_sym)
  sym$stool <- "Symptomatic"

  study1_asym$study <- study_var1
  study2_asym$study <- study_var2
  asym <- rbind(study1_asym, study2_asym)
  asym$stool <- "Asymptomatic"

  return(rbind(sym, asym))
}


# co_strength_bar <- function(ci_plot){
#   # Plot the co strength plots wil the bar graphs as a plot instread of the fill in the matrix
#
#   if(!("pr_dist" %in% names(ci_plot))){
#     print("Boolean Broken")
#     # d_f <- ci_plot %>% mutate(
#     #   Dist = pmin(abs(d_f$actual - d_f$low), abs(d_f$actual - d_f$high)))
#     #
#     # d_f$Flag <- ifelse(
#     #   (d_f$actual > d_f$high | d_f$actual < d_f$low), 1, 0)
#     #
#     # d_f2 <- d_f %>% mutate(pr_dist = (1 * (Dist/abs(actual - average))) * sign(actual - average)) %>%
#     #   mutate(pr_dist = ifelse((Flag > 0 & actual > 0), pr_dist, NA))
#     ci_plot <- process_ci_df(ci_plot)
#   }
#
#   final_vis <- ci_plot %>%
#     ggplot(aes(x = paste(path1, path2), y = actual)) +
#     geom_point()
#   return(final_vis)
# }


get_diverging_color_palette <- function(){
  #'@title Get the Divering Color pallete for figures. Should be the up to date version
  #'
  #'@description  Get the current divergin color palette to use on all appropiate
  #'figures. The Basic scheme is orange <-> gray <-> pruple. This is assumed to be the best
  #'asthetic/color-blind-friandly version.
  #'
  #'@return Names list of  3 members contain

  return(list(
    high = "#998ec3", # Purple
    low = "#f1a340", # Orange
    mid = "#d6d6d6") # Gray
)
}

build_edges <- function(edge_list, Xn, Yn){
  edge_lines <- list()
  for(e in 1:nrow(edge_list)){

    origin <- edge_list[e,]$V1
    dest <- edge_list[e,]$V2

    edge_line <- list(
      type = "line",
      line = list(color = "#030303", width = 0.3),
      x0 = Xn[origin],
      y0 = Yn[origin],
      x1 = Xn[dest],
      y1 = Yn[dest]
    )

    edge_lines[[e]] <- edge_line

  }

  return(edge_lines)
}


symmetric_empty_df <- function(d_f, path1, path2){
  all_names <- union(d_f[[path1]], d_f[[path2]])
  print(all_names)
  name_m <- matrix(NA, length(all_names), length(all_names))
  tracker <- lower.tri(name_m, diag = F)


  for(x in 1:length(all_names)){
    for(y in 1:length(all_names)){
      # if(tracker[x, y]){
      name_m[x, y] <- paste(all_names[x], all_names[y], sep = "_")
      # }
    }
  }

  final <- as.data.frame(name_m)  %>% tidyr::pivot_longer(everything(), names_to = "new_col", values_drop_na = T) %>%
    tidyr::separate(value, c("Path1", "Path2"), sep = "_") %>% select(one_of(c("Path1", "Path2")))
  return(final)
}


symmetric_df <- function(df, path1, path2, target){
  #' Take the tidy CI table and expand to a matrix, filled with the values given in target
  #'
  #' @description Add description here.
  original <- df %>% select(one_of(c(path1, path2, target)))

  inverse <- df %>% select(one_of(c(path2, path1, target)))
  names(inverse) <- c(path1, path2, target)

  symmetric <- rbind(original, inverse)

  full <- symmetric %>% select(one_of(path1, path2, target)) %>% tidyr::spread(path2, target)

  row.names(full) <- full[[path1]]

  full_m <- full %>% select(-path1) %>% as.matrix()


  # full_m[lower.tri(full_m, diag = T)] <- NA

  final<- reshape2::melt(full_m, na.rm = T)

  names(final) <- c(path1, path2, target)

  # m_size <- length(union(df[[path1]], df[[path2]]))
  #
  # tracker <- upper.tri(matrix(NA, m_size, m_size))
  #
  # for(x in 1:m_size){
  #   for(y in 1:m_size){
  #     if(tracker[x, y]){
  #       symmetric[[target]][which(row_val == x & col_val == y)] <- NA
  #     }
  #   }
  # }

  return(symmetric)
}


symmetric_df_multiple <- function(d_f, path1, path2, target, bins){
  bin_levels <- unique(d_f[[bins]])
  # final <- data.frame(c(), c(), c(), c(), c())

  final <- symmetric_empty_df(d_f, path1, path2)

  all_names_p1 <- unique(d_f[[path1]])
  all_names_p2 <- unique(d_f[[path2]])
  # names(final) <- c(path1, path2, target, bins)
  for(b in 1:length(bin_levels)){
    df_subset <- d_f[which(d_f[[bins]] == bin_levels[b]),] %>% symmetric_df(path1, path2, target)
    # df_subset <- df_subset %>% mutate(bins = bin_levels[b])
    names(df_subset) <- c(path1, path2, paste("Bin: ", b))
    # df_subset$bothVals <- c(df_subset$bins, df_subset[[target]])
    # if(b == 1){
    #   final <- df_subset
    # }else{
    final <-left_join(final, df_subset, by = c(path1, path2))
    # }
    # final <- full_join(final)
    # final <- rbind(final, df_subset)
  }

  final <- final %>%
    tidyr::pivot_longer(cols = starts_with("Bin:"), names_to = "bins", values_to = target, values_drop_na = F) %>%
    mutate(bins = as.numeric(gsub("Bin: ", "", .$bins)))

  return(final)

}


melt_df <- function(d_f, path1, path2, target){

  final <- d_f %>% select(one_of(path1, path2, target)) %>% tidyr::spread(path2, target)

  return(final)

}


#' process_ci_df <- function(ci_plot, confidence_interval = 0.95){
#'   #' @title Calculate the Confidence interval and Significance on Configuration Results
#'   #'
#'   #' @description Calculate the confidence interval from standard deviation and the mean of the
#'   #' ensemble generated by the configuration model. Calculates the significance strength, measured
#'   #' on a scale from [0,1]. Includes a binary flag of whter the actual value falls ourstide the
#'   #' confidence interval. Can specify the interval to cover 99% and 95% of the distribution.
#'   #'
#'   #' @param ci_plot : A dataframe of containing the results from the configuration model produced
#'   #' in the python code. Must contain a column containing the "average", standard deviation (stdev).
#'   #'
#'   #' @param confidence_interval :  (optional) A floating point number (either 0.95 or 0.99) indicating
#'   #' the width of the interval.
#'   #'
#'   #' @details The confidence intervals are also calucalted within the python code for the configuration model,
#'   #' but overwriting this column here allows for a little more flexibility. While the distribution is technically a poisson,
#'   #' I use a normal approximation for distributions with an average > 10, and a possion Chi Sqraured tool to
#'   #' find the High Density Region for a poisson distribution. Can change the level of significance to widen or narrrow
#'   #' that high density interval.
#'   #'
#'   #' @return The original dataframe with 5 additional columns: One indicating the higher and lower intervals, a binary
#'   #' indicator of significance, the distance outside that interval, and a normalised distance in the range of [0,1].
#'
#'   if("c_mean" %in% names(ci_plot)){
#'     c_score <- TRUE
#'     names(ci_plot)[1:9] <- c("path1", "path2", "actual", "average", "stdev", "low", "high", "c_mean", "c_stdev")
#'   }else{
#'     names(ci_plot)[1:7] <- c("path1", "path2", "actual", "average", "stdev", "low", "high")
#'     c_score  <- FALSE
#'   }
#'
#'
#'
#'   # This is the final version, if the average is below 10, then use the normal approximation, otherwise use the poisson
#'   # distribution.
#'   # ci_plot <- ci_plot %>% mutate(low = ifelse(average > 10, average - (1.96 * stdev), qchisq(0.025, 2* average)/2),
#'   #                               high = ifelse(average > 10, average + (1.96 * stdev), qchisq(0.975, 2* (average + 1)/2)))
#'   #
#'
#'   # Set the confidence interval
#'   if(confidence_interval == 0.95){
#'     normal_ci <- 1.96
#'   }else if(confidence_interval == 0.99){
#'     normal_ci <- 2.576
#'   }else{
#'     stop("Incorrect Confidence interval provided, options are 0.95 and 0.99")
#'   }
#'
#'   poisson_range <- confidence_interval / 2
#'
#'
#'
#'   ci_plot <- ci_plot %>% mutate(low = ifelse(average > 10, average - (normal_ci * stdev), stats::qpois(0.5 - poisson_range, average)),
#'                                 high = ifelse(average > 10, average + (normal_ci * stdev), stats::qpois(0.5 + poisson_range, average)))
#'
#'
#'   if(c_score){
#'     ci_plot <- ci_plot %>% mutate(c_low = c_mean - (normal_ci * c_stdev),
#'                                   c_high = c_mean + (normal_ci * c_stdev))
#'   }
#'
#'
#'   # ci_plot <- ci_plot %>% mutate(low = stats::qpois(0.025, average), high = stats::qpois(0.975, average))
#'   #
#'   # simple_results_maled %>% mutate(low_p = qchisq(0.025, 2* average)/2, high_p = qchisq(0.975, 2* (average + 1)/2))
#'
#'   # ci_plot <- ci_plot %>% mutate(low = exactPoiCI(average, level = "lower"), high = exactPoiCI(average, level = "higher"))
#'
#'   # Distance outside the interval
#'   d_f <- ci_plot %>% mutate(Dist = pmin(abs(ci_plot$actual - ci_plot$low), abs(ci_plot$actual - ci_plot$high)))
#'
#'   # Binary indicator for significance
#'   d_f$Flag <- ifelse(
#'     (d_f$actual > d_f$high | d_f$actual < d_f$low), 1, 0)
#'
#'   #  A nromalised distance indicator
#'   d_f2 <- d_f %>% mutate(
#'     pr_dist = (1 * (Dist/abs(actual - average))) * sign(actual - average)) %>%
#'     mutate(pr_dist = ifelse((Flag > 0 & actual > 0), pr_dist, NA))
#'
#'   d_f2 <- d_f2 %>% mutate(ci_width = high - low)
#'
#'   return(d_f2)
#' }


process_ci_df <- function(ci_plot, confidence_interval = 0.95, central_measurment = "mean"){
  #' @title Calculate the Confidence interval and Significance on Configuration Results
  #'
  #' @description Calculate the confidence interval from standard deviation and the mean of the
  #' ensemble generated by the configuration model. Calculates the significance strength, measured
  #' on a scale from [0,1]. Includes a binary flag of whter the actual value falls ourstide the
  #' confidence interval. Can specify the interval to cover 99% and 95% of the distribution.
  #'
  #' @param ci_plot : A dataframe of containing the results from the configuration model produced
  #' in the python code. Must contain a column containing the "average", standard deviation (stdev).
  #'
  #' @param confidence_interval :  (optional) A floating point number (either 0.95 or 0.99) indicating
  #' the width of the interval.
  #'
  #' @details The confidence intervals are also calucalted within the python code for the configuration model,
  #' but overwriting this column here allows for a little more flexibility. While the distribution is technically a poisson,
  #' I use a normal approximation for distributions with an average > 10, and a possion Chi Sqraured tool to
  #' find the High Density Region for a poisson distribution. Can change the level of significance to widen or narrrow
  #' that high density interval.
  #'
  #' @return The original dataframe with 5 additional columns: One indicating the higher and lower intervals, a binary
  #' indicator of significance, the distance outside that interval, and a normalised distance in the range of [0,1].

  names(ci_plot)[1:7] <- c("path1", "path2", "actual", "average", "stdev", "low", "high")

  if("c_mean" %in% names(ci_plot)){
    c_score <- TRUE
    # names(ci_plot)[ncol(ci_plot) - 1:ncol(ci_plot)] <- c("c_mean", "c_stdev")
  }else{
    c_score  <- FALSE
  }

  # This is the final version, if the average is below 10, then use the normal approximation, otherwise use the poisson
  # distribution.
  # ci_plot <- ci_plot %>% mutate(low = ifelse(average > 10, average - (1.96 * stdev), qchisq(0.025, 2* average)/2),
  #                               high = ifelse(average > 10, average + (1.96 * stdev), qchisq(0.975, 2* (average + 1)/2)))
  #

  # Set the confidence interval
  if(confidence_interval == 0.95){
    normal_ci <- 1.96
  }else if(confidence_interval == 0.99){
    normal_ci <- 2.576
  }else{
    stop("Incorrect Confidence interval provided, options are 0.95 and 0.99")
  }

  poisson_range <- confidence_interval / 2



  ci_plot <- ci_plot %>% mutate(low = ifelse(average > 10, average - (normal_ci * stdev), stats::qpois(0.5 - poisson_range, average)),
                                high = ifelse(average > 10, average + (normal_ci * stdev), stats::qpois(0.5 + poisson_range, average)))


  if(c_score){
    ci_plot <- ci_plot %>% mutate(c_low = c_mean - (normal_ci * c_stdev),
                                  c_high = c_mean + (normal_ci * c_stdev))
  }


  # ci_plot <- ci_plot %>% mutate(low = stats::qpois(0.025, average), high = stats::qpois(0.975, average))
  #
  # simple_results_maled %>% mutate(low_p = qchisq(0.025, 2* average)/2, high_p = qchisq(0.975, 2* (average + 1)/2))

  # ci_plot <- ci_plot %>% mutate(low = exactPoiCI(average, level = "lower"), high = exactPoiCI(average, level = "higher"))

  # Distance outside the interval
  d_f <- ci_plot %>% mutate(Dist = pmin(abs(ci_plot$actual - ci_plot$low), abs(ci_plot$actual - ci_plot$high)))

  # Binary indicator for significance
  d_f$Flag <- ifelse(
    (d_f$actual > d_f$high | d_f$actual < d_f$low), 1, 0)

  #  A nromalised distance indicator
  d_f2 <- d_f %>% mutate(
    pr_dist = (1 * (Dist/abs(actual - average))) * sign(actual - average)) #%>%
    # mutate(pr_dist = ifelse((Flag > 0 & actual > 0), pr_dist, NA))

  d_f2 <- d_f2 %>% mutate(ci_width = high - low)

  if(c_score){
    d_f2 <- get_c_actual(d_f2)

    d_f2$c_Flag <- ifelse(
      (d_f2$c_actual > d_f2$c_high | d_f2$c_actual < d_f2$c_low), 1, 0)

    d_f2 <- d_f2 %>% mutate(c_dist = pmin(abs(c_actual - c_low), abs(c_actual - c_high)) * sign(c_actual - c_mean))
  }


  return(d_f2)
}



exactPoiCI <- function (X, conf.level=0.95, level = "lower") {
  alpha = 1 - conf.level
  if(level == "lower"){
    return(0.5 * qchisq(alpha/2, (2*X +2)))
  }else if(level == "higher"){
    return(0.5 * qchisq((1-(alpha/2)), (2*X)))
  }else{
    stop("ERROR Improper")
  }
  # upper <- 0.5 * qchisq((1-(alpha/2)), (2*X))
  # return(c(lower, upper))
}



build_ci_plot <- function(d_f, path1, path2, target, bins, multiple = F){
  #' Build a CI plot with pploty for Configuration model ensemble.
  #'
  #' @description This function takes the dataframe generated by the enfoguration model
  #' ensemble generated in R to make a plot similar to a colored correlation matrix.
  #'
  #' @param d_f: the dataframe generated by the CI experiment function in rewire.R
  #' @param path1: A string that's the column name containing the names fo the first pathogen
  #' @param path2: see path1, for the second pathogen
  #' @param target: a string as the column header for the variable to color fill for.
  #' @param multiple: where this is a binned timeseries, while will draw a slider bar, building a plot
  #' for each time frame.
  #'
  #' @return a ggplot object, a geom tile plot with frame values to be passed to a ggplot function

  # Create the proportional distance
  d_f2 <- d_f %>% mutate(pr_dist2 = (1 * (Dist/abs(actual_col - mean_col))) * sign(actual_col - mean_col))

  # Take the dataframe, in tidy format, and convert to a motrix, filled with the target values given
  # in the arguments.
  if(multiple){
    melted_df <- symmetric_df_multiple(d_f2, path1, path2, target, bins)
  }else{
    melted_df <- symmetric_df(d_f2, path1, path2, target)
  }

  # Get all the names in the series.
  all_names <- union(d_f[[path1]], d_f[[path2]])

  all_names_y <- sort(all_names, decreasing = T)
  all_names_x <- sort(all_names, decreasing = F)
  # final <- ggplot(melted_df, aes(x =melted_df[[path1]], y = melted_df[[path2]], frame = bins)) +
  #   geom_tile(aes(fill = melted_df[[target]])) +
  #
  final <- melted_df %>% ggplot(aes(x =Path2, y = Path1, frame = bins)) +
    geom_tile(aes(fill = pr_dist2)) +
    coord_fixed()+
    # facet_wrap(.~bins)+
    # scale_y_reverse()+
    scale_x_discrete(limits = all_names_x, drop = FALSE)+
    scale_y_discrete(limits = all_names_y, drop = FALSE)+
    scale_fill_gradient2()+
    theme_dark()+
    # labs(x = path1, y = path2, fill = "ProDist")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1))

  return(final)

}


# Build the ci_matrix for the Vacc-generated results.
vacc_ci_plot <- function(d_f, path1, path2, target, bins, multiple = F, plot_title = NA, joined = F, n_panels = 2){
  #' Build a CI plot with pploty for Configuration model ensemble.
  #'
  #' @description This function takes the dataframe generated by the enfoguration model
  #' ensemble generated in R to make a plot similar to a colored correlation matrix.
  #'
  #' @param d_f: the dataframe generated by the CI experiment function in rewire.R
  #' @param path1: A string that's the column name containing the names fo the first pathogen
  #' @param path2: see path1, for the second pathogen
  #' @param target: a string as the column header for the variable to color fill for.
  #' @param multiple: where this is a binned timeseries, while will draw a slider bar, building a plot
  #' for each time frame.
  #'
  #' @return a ggplot object, a geom tile plot with frame values to be passed to a ggplot function

  # Create the proportional distance

  if(!("pr_dist" %in% names(d_f))){
    if(multiple){
      names(d_f)[1:8] <- c("path1", "path2", "actual", "average", "stdev", "low", "high", "bins")

    }else{
      names(d_f)[1:7] <- c("path1", "path2", "actual", "average", "stdev", "low", "high")
    }

    # d_f <- d_f %>% mutate(
    #   Dist = pmin(abs(d_f$actual - d_f$low), abs(d_f$actual - d_f$high)))
    #
    # d_f$Flag <- ifelse(
    #   (d_f$actual > d_f$high | d_f$actual < d_f$low), 1, 0)
    #
    # d_f2 <- d_f %>% mutate(pr_dist = (1 * (Dist/abs(actual - average))) * sign(actual - average)) %>%
    #   mutate(pr_dist = ifelse((Flag > 0 & actual > 0), pr_dist, NA))

    d_f2 <- process_ci_df(d_f)
  }else{
    d_f2 <- d_f
    if(!multiple){
      d_f2$bins <- 1
    }
  }

  # Get all the names in the series.
  all_names <- union(d_f[[path1]], d_f[[path2]])

  all_names_y <- sort(all_names, decreasing = T)
  all_names_x <- sort(all_names, decreasing = F)

  if(joined){
    if(n_panels == 2){
      final <- d_f2 %>% ggplot(aes(x =path1, y = path2, frame = bins)) +
        geom_tile(aes(fill = pr_dist)) +
        coord_fixed()+
        facet_wrap(.~study)+
        # scale_y_reverse()+
        scale_x_discrete(limits = all_names_x, drop = FALSE)+
        scale_y_discrete(limits = all_names_y, drop = FALSE)+
        scale_fill_gradient2()+
        theme_dark()+
        labs(x = path1, y = path2, fill = "ProDist", title = plot_title)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1))
      return(final)

    }else if(n_panels == 4){
      final <- d_f2 %>% ggplot(aes(x =path1, y = path2, frame = bins)) +
        geom_tile(aes(fill = pr_dist)) +
        coord_fixed()+
        facet_grid(study~stool)+
        # scale_y_reverse()+
        scale_x_discrete(limits = all_names_x, drop = FALSE, expand = c(0, 0))+
        scale_y_discrete(limits = all_names_y, drop = FALSE, expand = c(0, 0))+
        scale_fill_gradient2()+
        theme_dark()+
        labs(x = "", y = "", fill = "ProDist", title = plot_title)+
        theme(axis.text.x = element_text(angle = 45,
                                         size = rel(0.8),
                                         hjust = 1.2, vjust = 1.1),
              axis.text.y = element_text(
                size = rel(0.8)
                ),
              legend.margin=margin(t = 10),
              legend.position = c(0.95, 0.95))
    }

  }else{
    final <- d_f2 %>% ggplot(aes(x =path1, y = path2, frame = bins)) +
      geom_tile(aes(fill = pr_dist2)) +
      # geom_tile(aes(fill = weighted_interaction)) +
      coord_fixed()+
      # facet_wrap(.~bins)+
      # scale_y_reverse()+
      scale_x_discrete(limits = all_names_x, drop = FALSE)+
      scale_y_discrete(limits = all_names_y, drop = FALSE)+
      scale_fill_gradient2()+
      theme_dark()+
      labs(x = path1, y = path2, fill = "ProDist", title = plot_title)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1))

    return(final)

  }

}


filter_out_nonsymmetric <- function(d_f){
  #' Filter out redundant paris from symmetric matrix.
  #'
  #' @description Iterate through the upper triangle in a symmetric matrix and remove those
  #' entries from the results tables. Removes redudant pairs. Performs a row-wise deletion.
  #'
  #' @param d_f: The dataframe of the results from the configuration model, containing the actual,
  #' the mean, and the condifence interval.
  #'
  #' @return A dataframe with the redundant pairs removed.

  # Get all the unique pairs of pathogens.
  path_names <- unique(d_f$path1)
  # Use a matrix with the upper triangle (diagonal included) set to FALSE
  tracker <- upper.tri(matrix(nrow = length(path_names), ncol = length(path_names)), diag = T)
  # If this is an empty dataframe, then return the original empty dataframe
  if(nrow(d_f) == 0){
    return(d_f)
  }
  # Make a copy
  d_f_copy <- d_f
  # Iterate through the names, if the tracker matrix at that position is set to true,, filter
  # out all the pathogens that match both names.
  for(x in 1:length(path_names)){
    for(y in 1:length(path_names)){
      if(tracker[x, y]){
        d_f_copy <- d_f_copy %>% filter(path1 != path_names[x] | path2 != path_names[y])
      }
    }
  }

  return(d_f_copy)

}


filter_out_redundancies <- function(d_f){
  #' Filter out some of thepathogen pair redundancies
  #'
  #'  @description Filter out the pairs that might be considered redudant, such
  #'  as campy pan and campy c.jejuni, especially for visualisations, since many redudant
  #'  pairs are considered significant.
  #'
  #'  @param d_f: Results dataframe with the results from the configuration model such as
  #'  actual value, ensemble mean, and confidence interval.
  #'
  #'  @return A dataframe with the rows containing the redundant pairs removed.
  bad_rows <- c()
  # If this is an empty dataframe, return the original empty dataframe
  if(nrow(d_f) == 0){
    return(d_f)
  }
  # For every row in the dataframe, if any of these conditions are met, record the index.
  for(r in 1:nrow(d_f)){
    if(grepl("opv", d_f$path1[r]) & grepl("opv", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }else if(grepl("etec", d_f$path1[r]) & grepl("etec", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }else if(grepl("campy pan", d_f$path1[r]) & grepl("c.jejuni", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }else if(grepl("c.jejuni", d_f$path1[r]) & grepl("campy pan", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }else if(grepl("crypto", d_f$path1[r]) & grepl("crypto", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }else if(grepl("epec", d_f$path1[r]) & grepl("epec", d_f$path2[r])){
      bad_rows <- c(bad_rows, r)
    }
  }
  # Return the dataframe with the bad rows removed.
  return(d_f[-bad_rows,])
}


co_strength_bar <- function(ci_plot){
  # Plot the co strength plots wil the bar graphs as a plot instread of the fill in the matrix

  if(!("pr_dist" %in% names(ci_plot))){
    print("Booelan fail")

    names(ci_plot)[1:7] <- c("path1", "path2", "actual", "average", "stdev", "low", "high")
    d_f2 <- process_ci_df(ci_plot)
    d_f2 <- d_f2 %>% mutate(ci_width = high - low)
  }else{
    d_f2 <- ci_plot
  }


  # d_f <- ci_plot %>% mutate(Dist = pmin(abs(ci_plot$actual - ci_plot$low), abs(ci_plot$actual - ci_plot$high)))
  #
  # d_f$Flag <- ifelse(
  #   (d_f$actual > d_f$high | d_f$actual < d_f$low), 1, 0)
  #
  # d_f2 <- d_f %>% mutate(pr_dist = (1 * (Dist/abs(actual - average))) * sign(actual - average)) %>%
  #   mutate(pr_dist = ifelse((Flag > 0 & actual > 0), pr_dist, NA))
  #
  # # View(d_f2)

  print(names(d_f2))

  final_vis <- d_f2 %>%
    filter_out_nonsymmetric() %>%
    arrange(abs(high - average)) %>%
    top_n(100, actual) %>%
    mutate(color_flag = tanh(actual - average)) %>%
    ggplot(aes(x = reorder(paste(path1, path2, sep = " - "), abs(actual - average)))) +
    geom_point(aes(y = actual))+
    geom_errorbar(aes(ymin = low, ymax = high))+
    # geom_bar(aes(y = abs(actual - average), fill = color_flag), stat = "identity", color = "black")+
    geom_bar(aes(y = weighted_interaction, fill = color_flag), stat = "identity", color = "black")+
    theme(axis.text = element_text(size = rel(0.7)))+
    scale_fill_gradient2(midpoint = 0, low = 'red', high = 'blue')+
    coord_flip()#+
  # theme()
  return(final_vis)
}


sub_function <- function(x){
  x <- unlist(x)
  x <- ifelse(x == "Diarrhea", 1, 0)
  return(x)
}



get_diar_fraction <- function(d_f, source_df, study, tp = NA){

  study_specific <- ifelse((study == "provide") | (study == "maled"), T, F)

  for(col_name in 1:length(names(source_df))){
    names(source_df)[col_name] <- gsub("\\s+$", "", names(source_df)[col_name])
  }

  d_f$diar_co <- NA
  d_f$co_all <- NA
  d_f$expected_diar <- NA
  if(study_specific){
    d_f$diar_fraction <- source_df[which(source_df$study == study & source_df$stool_type == "Diarrhea"),] %>% nrow() /
      source_df[which(source_df$study == study),] %>% nrow()
  }else{
    d_f$diar_fraction <- source_df[which(source_df$stool_type == "Diarrhea"),] %>% nrow() /
      nrow(source_df)
  }
  if(("bins" %in% names(d_f)) |("timepoint" %in% names(d_f) & sum(is.na(d_f[['timepoint']])) == 0)){
    d_f$time_point <- NA
    d_f$n_obs <- NA
  }else{
    if(study_specific){
      d_f$n_obs <- source_df[source_df$study == study,] %>% nrow()
    }else{
      d_f$n_obs <- nrow(source_df)
    }
  }
  # d_f_subset <- d_f %>% filter(Falg > 0)
  for(r in 1:nrow(d_f)){
    # if(("bins" %in% names(d_f)) | ("timepoint" %in% names(d_f) & sum(is.na(d_f[['timepoint']])) == 0)){
    if(!is.na(tp)){
      d_f$diar_co[r] <- source_df[which(source_df$study == study &
                                          source_df$country == d_f$country[r] &
                                          source_df[[d_f$path1[r]]] > 0 &
                                          source_df[[d_f$path2[r]]] > 0 &
                                          source_df$stool_type == "Diarrhea" &
                                          source_df[[tp]] == d_f$timepoint[r]),] %>%
        # distinct(participant) %>%
        nrow()

      d_f$co_all[r] <- source_df[which(source_df$study == study &
                                         source_df$country == d_f$country[r] &
                                         source_df[[d_f$path1[r]]] > 0 &
                                         source_df[[d_f$path2[r]]] > 0 &
                                         source_df[[tp]] == d_f$timepoint[r]),] %>%
        # distinct(participant) %>%
        nrow()

      # Why are there distinct observations for building these graphs?
      d_f$n_obs[r] <- source_df[which(source_df$study == study &
                                        source_df$country == d_f$country[r] &
                                        source_df[[tp]] == d_f$timepoint[r]),] %>%
        # distinct(participant) %>%
        nrow()

      d_f$time_point[r] <- d_f$timepoint[r]

      if(!study_specific){
        d_f$diar_co[r] <- source_df[which(
          source_df[[d_f$path1[r]]] > 0 &
            source_df[[d_f$path2[r]]] > 0 &
            source_df$country == d_f$country[r] &
            source_df$stool_type == "Diarrhea" &
            source_df[[tp]] == d_f$timepoint[r]),] %>%
          # distinct(participant) %>%
          nrow()

        d_f$co_all[r] <- source_df[which(
          source_df[[d_f$path1[r]]] > 0 &
            source_df[[d_f$path2[r]]] > 0 &
            source_df$country == d_f$country[r] &
            source_df[[tp]] == d_f$timepoint[r]),] %>%
          # distinct(participant) %>%
          nrow()

        # Why are there distinct observations for building these graphs?
        d_f$n_obs[r] <- source_df[which(source_df[[tp]] == d_f$timepoint[r] & source_df$country == d_f$country[r]),] %>%
          # distinct(participant) %>%
          nrow()
      }

    }else{
      # if(d_f$Flag[r] > 0){
      d_f$diar_co[r] <- source_df[which(source_df$study == study &
                                          source_df[[d_f$path1[r]]] > 0 &
                                          source_df[[d_f$path2[r]]] > 0 &
                                          source_df$stool_type == "Diarrhea"),] %>% nrow()

      d_f$co_all[r] <- source_df[which(source_df$study == study &
                                         source_df[[d_f$path1[r]]] > 0 &
                                         source_df[[d_f$path2[r]]] > 0),] %>% nrow()

      d_f$expected_diar[r] <- d_f$co_all[r] * d_f$diar_fraction[r]
      # }
    }


    # if(!study_specific){
    #   d_f$diar_co[r] <- source_df[which(
    #     source_df[[d_f$path1[r]]] > 0 &
    #       source_df[[d_f$path2[r]]] > 0 &
    #       source_df$stool_type == "Diarrhea"),] %>%
    #     # distinct(participant) %>%
    #     nrow()
    #
    #   d_f$co_all[r] <- source_df[which(
    #     source_df[[d_f$path1[r]]] > 0 &
    #       source_df[[d_f$path2[r]]] > 0),] %>%
    #     # distinct(participant) %>%
    #     nrow()
    #
    #   # Why are there distinct observations for building these graphs?
    #   # d_f$n_obs[r] <- source_df[which(source_df[[tp]] == d_f$timepoint[r]),] %>%
    #   #   # distinct(participant) %>%
    #   #   nrow()
    # }

  }
  return(d_f)
}


get_bangladesh_seasons <- function(){
  seasons <- list(
    'basanta' = list(
      "start" = as.Date("2020-02-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-04-15", format = "%Y-%m-%d"),
      'id' = 0,
      'name' = "basanta"
    ),
    'grisma' = list(
      "start" = as.Date("2020-04-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-06-15", format = "%Y-%m-%d"),
      'id' = 1,
      "name" = "grisma"
    ),
    'barsa' = list(
      "start" = as.Date("2020-06-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-08-15", format = "%Y-%m-%d"),
      'id' = 2,
      "name" = "barsa"
    ),
    'sharat' = list(
      "start" = as.Date("2020-08-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-10-15", format = "%Y-%m-%d"),
      'id' = 3,
      "name" = "sharat"
    ),
    'hemanta' = list(
      "start" = as.Date("2020-10-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-12-15", format = "%Y-%m-%d"),
      'id' = 4,
      "name" = "hemanta"
    ),
    'shit' = list(
      "start" = as.Date("2020-12-15", format = "%Y-%m-%d"),
      "end" = as.Date("2020-02-15", format = "%Y-%m-%d"),
      'id' = 5,
      "name" = "shit"
    )
  )
  return(seasons)
}

assign_season <- function(d_f){
  d_f <- d_f %>% mutate(month_day = as.Date(paste(format(d_f$collection_date, "%m"), format(d_f$collection_date, "%d"), "2020", sep = "-"), format = "%m-%d-%Y"))
  seasons <- get_bangladesh_seasons()
  d_f$season <- NA
  for(s in seasons){
    if(s[['name']] == "shit"){
      d_f$season[which((s[['start']] <= d_f$month_day) | (d_f$month_day < s[['end']]))] <- s[['id']]
    }else{
      d_f$season[which((s[['start']] <= d_f$month_day) & (d_f$month_day < s[['end']]))] <- s[['id']]
    }
  }
  return(d_f)
}



############### i got another code version working, so I'm mightr delete this 10Jan20 #############


#' add_vertices_tograph <- function(g, verts){
#'   for(v in verts){
#'     add_vertices(g, 1, label = v, color = 'red')
#'     # print(verts[[v]])
#'   }
#'   return(g)
#'
#' }
#'
#' build_graph <- function(df){
#'   #' The main function for building the network. This gets called from
#'   #' the notebook.
#'   #'
#'   vertices <- get_vertices(df)
#'   g <- make_empty_graph(n = length(vertices), directed = F)
#'   g <- add_vertices_tograph(g, vertices)
#'   g <- add_pathogen_edges(g, df, vertices)
#'   return(g)
#' }
#'
#'
#' get_vertices <- function(df){
#'   #' To get all the vertices from the graph, looks at the unique instances of either
#'   #' Path1 or Path2
#'   unique_labels <- unique(c(df$Path1, df$Path2))
#'   label_list <- list()
#'   for(p in 1:length(unique_labels)){
#'     label_list[[unique_labels[p]]] <- p
#'   }
#'   return(label_list)
#'
#'
#' }
#'
#' add_pathogen_edges <- function(g, df, v_list){
#'   #'Adding edges between vertices, loops through rows in the dataframe,
#'   #'adds edges for each row. Row edges are added as red if they are
#'   #'below the confidence interval, and blue edges are above the
#'   #'confidence interval.
#'   #'
#'   for(i in 1:nrow(df)){
#'     g <- add_edges(g, c(v_list[[df$Path1[i]]], v_list[[df$Path2[i]]]),
#'                    attr = list(color = ifelse(df$actual_col[i] < df$low_CI[i], "red", "blue"), weight = df$actual_col[i]))
#'   }
#'   return(g)
#' }


readable_path_names <- function(d_f, formatted = F){
  #'@title Change pathogen variable names
  #'
  #'@description  Change the pathogen names from the variable names, to more the figure friendly
  #'names using a neamed list (R's version of a dictionary)
  #'
  #'@param d_f the dataframe, as the result from the configuration model, usually called within
  #'plotting functions in the package.
  #'
  #'@param formatted : Boolean, should the names be formatted appropiately for scientific names? The
  #'default is FALSE, where there is no formatting. Otherwise, must be used in concert with \code{element_markdown()}
  #'from the \code{ggtext} package.
  #'
  #'@return dataframe, d_f with only the pathogen names changed.

  if(!formatted){
    path_names <- c(
      "ascaris lumbricoides"          = "A. lumbroides",
      "aeromonas"                     = "Aeromonas",
      "ancyclostoma"                  = "Ancyclostoma",
      "trichuris trichiura"           = "Trichuris",
      "e.bieneusi"                    = "E. bieneusi",
      "e.intestinalis"                = "E. intestinalis",
      "cryptosporidium"               = "Cryptosporidium spp.",
      "salmonella"                    = "Salmonella",
      "h.pylori"                      = "H. pylori",
      "c.jejuni/coli"                 = "C. jejuni/coli",
      "campy pan"                     = "Campylobacter spp.",
      "b.fragilis"                    = "B. fragilis",
      "c.difficile"                   = "C. difficile",
      "adenovirus f"                  = "Adenovirus 40/41",
      "norovirus gi"                  = "Norovirus GI",
      "norovirus gii"                 = "Norovirus GII",
      "astrovirus"                    = "Astrovirus",
      "necator"                       = "Necator",
      "strongyloides"                 = "Strongyloides",
      "cyclospora"                    = "Cyclospora",
      "isospora"                      = "Isospora",
      "e.histolytica"                 = "E. histolytica",
      "m.tb"                          = "M. tuberculosis",
      "v.cholerae"                    = "V. cholerae",
      "shigella & eiec"               = "Shigella spp.",
      "sapovirus"                     = "Sapovirus",
      "rotavirus"                     = "Rotavirus",
      "eaec"                          = "EAEC",
      "atypical epec"                 = "aEPEC",
      "typical epec"                  = "tEPEC",
      "stec"                          = "STEC",
      "lt_etec"                       = "ETEC lt",
      "st_etec"                       = "ETEC st",
      "etec"                          = "ETEC",
      "epec"                          = "EPEC"
    )
  }else{
    path_names <- c(
      "ascaris lumbricoides"          = "*A. lumbroides*",
      "aeromonas"                     = "*Aeromonas*",
      "ancyclostoma"                  = "*Ancyclostoma*",
      "trichuris trichiura"           = "*Trichuris*",
      "e.bieneusi"                    = "*E. bieneusi*",
      "e.intestinalis"                = "*E. intestinalis*",
      "cryptosporidium"               = "*Cryptosporidium* spp.",
      "salmonella"                    = "*Salmonella*",
      "h.pylori"                      = "*H. pylori*",
      "c.jejuni/coli"                 = "*C. jejuni/coli*",
      "campy pan"                     = "*Campylobacter* spp.",
      "b.fragilis"                    = "*B. fragilis*",
      "c.difficile"                   = "*C. difficile*",
      "adenovirus f"                  = "Adenovirus 40/41",
      "norovirus gi"                  = "Norovirus GI",
      "norovirus gii"                 = "Norovirus GII",
      "astrovirus"                    = "Astrovirus",
      "necator"                       = "*Necator*",
      "strongyloides"                 = "*Strongyloides*",
      "cyclospora"                    = "*Cyclospora*",
      "isospora"                      = "*Isospora*",
      "e.histolytica"                 = "*E. histolytica*",
      "m.tb"                          = "*M. tuberculosis*",
      "v.cholerae"                    = "*V. cholerae*",
      "shigella & eiec"               = "*Shigella* spp.",
      "sapovirus"                     = "Sapovirus",
      "rotavirus"                     = "Rotavirus",
      "eaec"                          = "EAEC",
      "atypical epec"                 = "aEPEC",
      "typical epec"                  = "tEPEC",
      "stec"                          = "STEC",
      "lt_etec"                       = "ETEC lt",
      "st_etec"                       = "ETEC st",
      "etec"                          = "ETEC",
      "epec"                          = "EPEC"
    )

  }

  final <- d_f %>% mutate(
    path1 = plyr::revalue(path1, path_names, warn_missing = F),
    path2 = plyr::revalue(path2, path_names, warn_missing = F)
  )

  return(final)
}


plot_simple_results <- function(results, mymaintitle, mytitle1, mytitle2,
                                fig_stool_type,
                                bar_alpha = 0.4,
                                error_bar_size = 0.9,
                                actual_limit = 5,
                                point_size = 1.5,
                                axis_text_size = 1.5,
                                axis_title = 1.6,
                                plot_title = 1.7,
                                sup_title_size = 30){
  #' Plot Configuration Results for a by-observation model
  #'
  #' @description Plot the results from running the configuration model on a per-observation basis. There is no
  #' timserires componenet here. On the left is the bar plots showing the proportion of stools with that stool type and on the
  #' left is the error+point plot showing the significance strength. Shows all of the significant pairs, with redundent and symmetrical
  #' pairs removed.
  #'
  #' @param results The results from the confiruamtion model, contains results from all 4 data subsets (MAL-ED, PROVIDE, Diarrhea, Asymptomatic).
  #' @param mymaintitle The Super title for the entire facet grid.
  #' @param mytitle1 The title on the left, with the bars.
  #' @param mytitle2 The title on the right, with the point and the errorbars showing significance strength.
  #' @param fig_stool_type The stool type, used to determine what fills the bars in the bar plots.
  #' @param error_bar_size The thickness of the errorbars, defaults is 0.9

  bar_fill_alpha <- bar_alpha
  plot_base <- results %>%
    filter_out_nonsymmetric() %>% filter(Flag > 0, actual > actual_limit) %>% filter_out_redundancies() %>%
    mutate(study = ifelse(study == "maled", "MAL-ED", "PROVIDE")) %>%
    mutate(color_flag = pr_dist
    ) %>%
    readable_path_names() %>%
    ggplot(aes(x = reorder(paste(path1, path2, sep = " + "), pr_dist)))

  my_diverging_c <- get_diverging_color_palette()

  if(fig_stool_type == "Diarrhea"){
    first_panel <- plot_base +
      geom_bar(aes(y = (diar_co/n_obs)*100, fill = color_flag),
               alpha = bar_fill_alpha, stat = "identity", color = "White")+
      geom_errorbar(aes(ymin = (diar_fraction) * (co_all/n_obs)*100, ymax = (diar_fraction) * (co_all/n_obs)*100), size = error_bar_size)#+
  }else if(fig_stool_type == "Asymptomatic"){
    first_panel <- plot_base +
      geom_bar(aes(y = (co_all - diar_co)/n_obs * 100,
                   fill = color_flag),
               alpha = bar_fill_alpha,
               stat = "identity",
               color = "White")+
      geom_errorbar(aes(ymin = (1-diar_fraction) * (co_all/n_obs)*100,
                        ymax = (1-diar_fraction) * (co_all/n_obs)*100), size = error_bar_size)

  }else{
    stop("Error: You need to give a stool type.")
  }

  grid.arrange(
    first_panel +
      geom_bar(aes(y = (co_all / n_obs) * 100, color = color_flag),
               alpha = .1, stat = "identity", fill = "white", size = 1)+
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))+
      scale_fill_gradient2(midpoint = 0,
                           low = my_diverging_c[['low']],
                           high = my_diverging_c[['high']],
                           mid = my_diverging_c[['mid']])+
      scale_color_gradient2(midpoint = 0,
                            low = my_diverging_c[['low']],
                            high = my_diverging_c[['high']],
                            mid = my_diverging_c[['mid']], guide = "none")+
      coord_flip()+

      theme(axis.text = element_text(size = rel(axis_text_size)),
            axis.title = element_text(size = rel(axis_title)),
            plot.title = element_text(size = rel(plot_title)),
            panel.background = element_blank(),
            panel.grid = element_line(color = "lightgray"),
            legend.position = "right")+
      labs(title = mytitle1,
           x = "", y = "Percent of all stools",
           fill = "Significance\nStrength"),

    results %>%
      filter_out_nonsymmetric() %>% filter(Flag > 0, actual > actual_limit) %>%
      filter_out_redundancies() %>%
      mutate(study = ifelse(study == "maled", "MAL-ED", "PROVIDE")) %>%
      mutate(color_flag = pr_dist) %>%
      readable_path_names() %>%
      ggplot(aes(x = reorder(paste(path1, path2, sep = " + "), pr_dist))) +
      geom_point(aes(y = actual), size = point_size)+
      geom_errorbar(aes(ymin = low, ymax = high), size = error_bar_size)+
      coord_flip()+
      theme(axis.text = element_text(size = rel(axis_text_size)),
            axis.title = element_text(size = rel(axis_title)),
            plot.title = element_text(size = rel(plot_title)),
            panel.background = element_blank(), panel.grid = element_line(color = "lightgray"))+
      labs(title = mytitle2,
           x = "", y = "Number of stool samples with co-occurrence",
           fill = "tanh( Actual  - Average)"),
    ncol = 2,
    top = grid::textGrob(mymaintitle, gp=grid::gpar(fontsize = sup_title_size))
  )

}

get_diar_barplot <- function(d_f, plot_title, y_axis_title, level, chosen_measurement = "count", fig_stool_type = "Diarrhea"){

  temp_df <- d_f %>%
    mutate(interaction_type = ifelse(
      (measurement == "count" & pr_dist >= 0) | (measurement == "c_score" & c_dist >= 0), "Higher", "Lower")) %>%
    filter(interaction_type == level) %>% filter_out_nonsymmetric() %>% filter_out_redundancies() %>%
    readable_path_names() %>%
    mutate(study = ifelse(study == "maled", "MAL-ED", "PROVIDE"))

  if(measurement == "count"){
    plot_obj <- ggplot(data = temp_df, aes(x = reorder(paste(path1, path2, sep = " + "), abs(pr_dist))))

  }else if(measurement == "c_score"){
    plot_obj <- ggplot(data = temp_df, aes(x = reorder(paste(path1, path2, sep = " + "), abs(c_dist))))
  }else{
    stop("ERROR: Improper measurement argument given")
  }

  shared_direction <- temp_df %>% group_by(path1, path2, interaction_type) %>%
    summarise(direction_size = n()) %>%
    filter(direction_size > 1) %>% select(path1, path2)

  temp_df <- inner_join(temp_df, shared_direction, b = c("path1", "path2"))

  if(fig_stool_type == "Diarrhea"){
    bar_fill_opts <- geom_bar(data = temp_df, aes(y = (diar_co / n_obs) * 100, fill = study),
                              alpha = .5, stat = "identity", color = "White",
                              position = "dodge", width = bar_width)
    error_bar_opts <- geom_errorbar(data = temp_df, aes(ymin = diar_fraction * (co_all/n_obs)*100,
                                                        ymax = diar_fraction * (co_all/n_obs)*100,
                                                        group = study),
                                    size = .9, position  =  "dodge", width = bar_width)

  }else if(fig_stool_type == "Asymptomatic"){
    bar_fill_opts <- geom_bar(data = temp_df, aes(y = ((co_all - diar_co) / n_obs) * 100, fill = study),
                              alpha = .5, stat = "identity", color = "White",
                              position = "dodge", width = bar_width)
    error_bar_opts <- geom_errorbar(data = temp_df, aes(ymin = (1 - diar_fraction) * (co_all/n_obs)*100,
                                                        ymax = (1 - diar_fraction) * (co_all/n_obs)*100,
                                                        group = study),
                                    size = .9, position  =  "dodge", width = bar_width)
  }


  final <- plot_obj + bar_fill_opts + error_bar_opts+
    geom_bar(data = temp_df, aes(y = (co_all / n_obs) * 100, color = study),
             alpha = .1, stat = "identity", fill = "white", size = 1,
             position = "dodge", width = bar_width)+
    scale_y_sqrt(
      breaks = seq(from = 0, to = 26, 4), labels = seq(from = 0, to = 26, 4), limits = c(0, 20)
    )+
    scale_fill_manual(values = two_var_coloring)+
    scale_color_manual(values = two_var_coloring, guide = FALSE)+
    coord_flip()+
    theme(axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.4)),
          plot.title = element_text(size = rel(1.5)),
          panel.background = element_blank(),
          panel.grid = element_line(color = "lightgray"),
          legend.position = "right",
    )+
    labs(title = plot_title,
         x = "", y = y_axis_title,
         color = "",
         fill = "Significance\nStrength")
  return(final)
}


get_sig_bars <- function(d_f, plot_title, y_axis_title, level){
  #' @title PLot Significant Pairs Shared in both studies
  #'
  #' @description Plot the significant pairs and the distribution of the 95%
  #' Credible interval and the actual number of observed co-occurences. Color is
  #' done by study and partitions by higher and lower. NOTE: The ability to plot
  #' the checkerboard score (C-Score) was removed and archived on 29Jan21.
  #'
  #' @param plot_title The title for the figure
  #' @param y_axis_title The title for the Y axis
  #' @param level String as the direction, c("Higher", "Lower")
  point_shape <- 19
  point_size <- 3
  bar_size <- 1.2

  d_f <-  find_shared_sig_pairs(d_f)

  temp_df <- d_f %>%
    mutate(interaction_type = ifelse(
      pr_dist >= 0, "Higher", "Lower")) %>%
    filter(interaction_type == level)%>%
    filter_out_nonsymmetric() %>% filter_out_redundancies() %>%
    readable_path_names() %>%
    mutate(study = ifelse(study == "maled", "MAL-ED", "PROVIDE"))

  shared_direction <- temp_df %>%
    group_by(path1, path2, interaction_type) %>%
    group_by(path1, path2) %>%
    summarise(direction_size = n()) %>%
    filter(direction_size > 1) %>% select(path1, path2)

  temp_df <- inner_join(temp_df, shared_direction, b = c("path1", "path2"))

  plot_obj <- ggplot(
    data = temp_df,
    aes(x = reorder(paste(path1, path2, sep = " + "), abs(pr_dist))))+
    geom_point(data = temp_df, aes(y = actual, color = study, group = study),
               shape = point_shape, size = point_size,
               position = position_dodge(width = 0.9))+
    geom_errorbar(
      data = temp_df,
      aes(ymin = low, ymax = high, color = study),
      size = bar_size,
      position = "dodge",
      width = bar_width)+
  scale_y_continuous(limits  =c(0,200))

  final <- plot_obj +
    coord_flip()+
    scale_color_manual(values = two_var_coloring)+
    theme(axis.text = element_text(size = rel(1.3)),axis.title = element_text(size = rel(1.4)),
          plot.title = element_text(size = rel(1.5)),
          panel.background = element_blank(), panel.grid = element_line(color = "lightgray"))+
    labs(title = plot_title,
         color = "",
         x = "", y = y_axis_title,
         fill = "tanh( Actual  - Average)")
  return(final)
}

# shared_pairs_data <- rbind(
#   inner_join(
#     all_intersections %>% filter(int_type == "c_score") %>% select(-value),
#     all_simple,
#     by = c("path1", "path2", "stool_type", "study")),
#
#   inner_join(
#     all_intersections %>% filter(int_type == "count") %>% select(-value),
#     all_simple,
#     by = c("path1", "path2", "stool_type", "study"),
#   ))


# Diar Results, Segmenting out by study

plot_shared_pairs_simple <- function(results, chosen_stool_type, chosen_direction, chosen_int_type = "c_score"){

  if(!(int_type %in% names(results))){
    results <- results %>% mutate(int_type = chosen_int_type)
  }

  sup_title <- paste0(chosen_stool_type, " Samples by Frequency (|P| < .05)")


  grid.arrange(
    get_diar_barplot(
      results %>% filter(stool_type == chosen_stool_type, int_type == chosen_int_type),
      "\nPercent of co-occurrence in stools\nHigher than Expected",
      "Percent of all stools",
      "Higher",
      fig_stool_type = chosen_stool_type,
      measurement = chosen_int_type),

    get_sig_bars(
      results %>% filter(stool_type == chosen_stool_type, int_type == chosen_int_type) ,
      "\nNumber of co-occurrence in stools\nHigher than Expected",
      "Occurence in stools",
      "Higher",
      measurement = chosen_int_type),

    get_diar_barplot(
      results %>% filter(stool_type == chosen_stool_type, int_type == chosen_int_type) ,
      "\nPercent of co-occurrence in stools\nLower than Expected",
      "Percent of all stools",
      "Lower",
      fig_stool_type = chosen_stool_type,
      measurement = chosen_int_type),

    get_sig_bars(
      results %>% filter(stool_type == chosen_stool_type, int_type == chosen_int_type) ,
      "\nNumber of co-occurrence in stools\nLower than Expected",
      "Occurence in stools",
      "Lower",
      measurement = chosen_int_type),

    ncol = 2,
    nrow = 2,
    top = grid::textGrob(sup_title, gp=grid::gpar(fontsize=30))
  )

}




