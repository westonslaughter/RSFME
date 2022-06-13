library(lubridate)
library(tidyverse)
library(lfstat)
# method description: mean weekly [con]*mean weekly streamflow

# create function of weekly flux estimation
estimate_flux_fernow_weekly <- function(chem_df, q_df, ws_size){
  # initialize output df
  startDate <- min(chem_df$date)
  endDate <- max(chem_df$date)
  
  dateRange <- seq(startDate, endDate, by = "days")
  out_df <- tibble(date = dateRange, con_est = NA) %>%
    full_join(., chem_df, by = 'date')
  
 # make weekly aggregation of chem and q
  weekly_chem_df <- chem_df %>%
    group_by(week = paste0(year(date), '-', week(date))) %>%
    summarize(date = min(date),
              mean_weekly_con = mean(con, na.rm = T))
  
  weekly_q_df <- q_df %>%
    filter(date >= startDate, # trim data to match chem
           date <= endDate) %>%
    group_by(week = paste0(year(date), '-', week(date))) %>%
    summarize(date = min(date),
              mean_weekly_q = mean(q_lps, na.rm = T))
  
  flux_df <- weekly_chem_df %>%
    full_join(weekly_q_df, by = 'week') %>%
    filter(!is.na(date.x)) %>%
    mutate( q_lpw = (mean_weekly_q*604800),
            flux_weekly_kg_ha = q_lpw*mean_weekly_con*1e-6/ws_size,
            method = 'fernow') %>%
    select(date = date.y, flux_weekly_kg_ha, method) %>%
      na.omit() %>%
      mutate(wy = water_year(date, origin = 'usgs'))
  return(flux_df)
  
}