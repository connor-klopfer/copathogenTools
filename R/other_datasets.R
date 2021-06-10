# Outside datasets

import_weather_data <- function(parent_dir){

  weather_data <- read.csv(paste0(parent_dir, "2198236.csv"), stringsAsFactors = F)
  weather_data$DATE <- as.Date(weather_data$DATE)

  weather_data_complete <- weather_data %>%
    tidyr::complete(DATE = seq.Date(min(DATE), max(DATE), by="day"))

  weather_year_ave <- weather_data_complete %>%
    mutate(month_day = as.Date(paste(format(DATE, "%m"), format(DATE, "%d"), sep = "-"), format = "%m-%d"),
           year = format(DATE, "%Y")) %>%
    group_by(month_day) %>%
    mutate(TAVG = tidyr::replace_na(TAVG, mean(TAVG, na.rm = T)),
           PRCP = tidyr::replace_na(PRCP, mean(PRCP, na.rm = T))) %>%
    mutate(TAVG = ifelse(is.nan(TAVG), NA, TAVG),
           PRCP = ifelse(is.nan(PRCP), NA, PRCP)) %>%
    ungroup() %>%
    tidyr::fill(TAVG, PRCP) %>%
    group_by(month_day) %>%
    summarise(ave_temp = mean(TAVG), ave_precp = mean(PRCP))


  return(weather_year_ave)
}
