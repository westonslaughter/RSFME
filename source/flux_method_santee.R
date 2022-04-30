library(lubridate)
library(tidyverse)
library(lfstat)
# method description: weekly/biweekly sample*cummlative q over that time, no interpolation

# create function of weekly flux estimation
estimate_flux_santee <- function(chem_df, q_df, ws_size){
  
 # Join chem and q
 join_df <- full_join(chem_df, q_df, by = 'date') %>%
   arrange(date)
 
 # make index of chem samples
 ind <- which(!is.na(join_df$con))
 
 # sum all q between samples and reduce
 out_df <- join_df %>%
   mutate(flux_kg_ha = NA)
 for(i in 1:length(ind)){
   con_mg_l <- join_df$con[i]
   sampleDate <- join_df$date[i]
   
   # handle different first case
   if(i == 1){
     q_total <- join_df %>%
       filter(date <= sampleDate) %>%
       summarize(q = sum(q_lps*900)) %>%
       .$q
   }else{
    #find previous sample date
    prevSampleDate <- 
    q_total <-  
   }
   
   #
   out_df$flux_kg_ha[i] <- (con_mg_l*q_total)/(ws_size*1e6)
 }
  
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
    select(date = date.y, flux_weekly_kg_ha, method)
  return(flux_df)
  
}