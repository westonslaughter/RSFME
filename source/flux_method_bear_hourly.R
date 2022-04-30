library(lubridate)
library(tidyverse)
library(lfstat)

# method description: linear of chem to match (hourly) discharge
source('source/helper_functions.R')

# create function of weekly flux estimation
estimate_flux_bear_hourly <- function(chem_df, q_df, ws_size){
  #make q hourly
  hourly_q <- q_df %>%
    group_by(hour = paste0(hour(date), '-', as_date(date))) %>%
    summarize(q_lph = sum(q_lps*900),
              date = min(date)) %>%
    select(date, q_lph)
  
  # interpolate con to match and calculate flux
  flux_df <- interpolate_con_to_q(chem_df = chem_df, 
                                  q_df = hourly_q) %>%
   mutate(flux_hourly_kg_ha = (con*q_lph)/(ws_size*1e6),
          method = 'bear') %>%
   select(date, flux_hourly_kg_ha, method)
  return(flux_df)
  
}