
#method description: q*con, with average con used between two points

# create function of hbef flux estimation
estimate_flux_hbef_daily <- function(chem_df, q_df, ws_size){
  # initialize output df
  startDate <- min(chem_df$date)
  endDate <- max(chem_df$date)
  
  dateRange <- seq(startDate, endDate, by = "days")
  out_df <- tibble(date = dateRange, con_est = NA) %>%
    full_join(., chem_df, by = 'date')
  
  # find indexes of measured values
  ind <- which(is.na(out_df$con)==FALSE)
  
  # create table of averages between measured values
  avg_df <- tibble(ind, average = (out_df$con[ind]+out_df$con[lead(ind)])/2)
  
  # populate average values between measured values
  for(i in 2:length(avg_df$ind)){
    out_df$con_est[avg_df$ind[i-1]:avg_df$ind[i]] <- avg_df$average[i-1]
    out_df$con_est[avg_df$ind[i-1]] <- out_df$con[avg_df$ind[i-1]]
    out_df$con_est[avg_df$ind[i]] <- out_df$con[avg_df$ind[i]]
  }
  
  # create daily q ts to match chem ts
  daily_q <- q_df %>%
    mutate(day = lubridate::date(date)) %>%
    group_by(day) %>%
    summarize(q_lps = mean(q_lps)) %>%
    mutate(q_lpd = q_lps*86400) %>%
    select(day, q_lpd) %>%
    filter(day >= startDate,
           day <= endDate) %>%
    rename(date = day)
  
  # compute daily flux
  flux_df <- out_df %>%
    full_join(., daily_q, by = 'date') %>%
    mutate(con_est = as.numeric(con_est),
           flux = con_est*q_lpd*1e-6/ws_size,
           method = 'hbef') %>%
    select(date, flux, method)

  return(flux_df)
}


