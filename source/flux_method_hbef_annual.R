# create function of hbef flux estimation for annual
estimate_flux_hbef_annual <- function(chem_df, q_df, ws_size){
  startDate <- min(chem_df$date)
  endDate <- max(chem_df$date)
  source('source/flux_method_hbef_daily.R')
  flux_df <- estimate_flux_hbef_daily(chem_df = chem_df, q_df = q_df, ws_size = ws_size)
  # compute total flux in kg/ha and return it
  out <- flux_df %>%
    group_by(wy = water_year(date, origin = "usgs")) %>%
    summarize(flux = sum(flux, na.rm = TRUE))
    
  
  out <- out %>%
    mutate(method = 'hbef')
  return(out)
}
