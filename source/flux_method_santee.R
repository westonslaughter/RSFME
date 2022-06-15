# method description: weekly/biweekly sample*cumulative q over that time, no interpolation

estimate_flux_santee <- function(chem_df, q_df, ws_size){
  
 # Join chem and q
 join_df <- full_join(chem_df, q_df, by = 'date') %>%
   arrange(date)
 
 # make index of chem samples
 ind <- which(!is.na(join_df$con))
 
 # sum all q between samples and reduce
 out_df <- join_df %>%
   mutate(flux_kg_ha = NA,
          periodStartDate = as_datetime('9999-09-09 09:09:09'),
          sampleDate = as_datetime('9999-09-09 09:09:09'))
 
 for(i in 1:length(ind)){
   con_mg_l <- join_df$con[ind[i]]
   sampleDate <- join_df$date[ind[i]]
   
   # handle different first case
   if(i == 1){
     q_total <- join_df %>%
       filter(date <= sampleDate) %>%
       summarize(q = sum(q_lps*900)) %>%
       .$q
     out_df$periodStartDate[ind[i]] <- join_df$date[[1]]
   }else{
    #find previous sample date
    prevSampleDate <- join_df$date[[ind[i-1]]] 
    q_total <-  join_df %>%
      filter(date <= sampleDate,
             date > prevSampleDate) %>%
      summarize(q = sum(q_lps*86400)) %>%
      .$q
    out_df$periodStartDate[ind[i]] <- prevSampleDate
   }
   
   # calculate flux over the period between samples
   out_df$flux_kg_ha[ind[i]] <- (con_mg_l*q_total)/(ws_size*1e6)
   out_df$sampleDate[ind[i]] <- sampleDate
 }

 flux_df <- out_df %>%
    mutate(method = 'santee',
           wy = water_year(date, origin = 'usgs')) %>%
    select(date, flux = flux_kg_ha, periodStartDate, sampleDate, method, wy) %>%
    filter(!is.na(flux))
 
 #make daily to match q assumptions
daily_out <- tibble(date = as.numeric(), flux = as.numeric())
for(i in 1:nrow(flux_df)){
    startDate <- flux_df$periodStartDate[i]
    endDate <- flux_df$sampleDate[i]
    daySeq <- seq(startDate, endDate, by = "days")
    
    agg_flux <- flux_df$flux[i]
    
    chunk_out <- tibble(date = daySeq, flux = agg_flux/length(daySeq))
    
    daily_out <- rbind(daily_out, chunk_out)
    
}
daily_out <- daily_out %>%
    mutate(method = 'santee', wy = water_year(date, origin = 'usgs'))
    
  return(daily_out)
  
}