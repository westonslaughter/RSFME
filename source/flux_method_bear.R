# method description: linear interp of chem to match discharge
source('source/helper_functions.R')
# create function of weekly flux estimation
estimate_flux_bear <- function(chem_df, q_df, ws_size){
  # interpolate con to match and calculate flux
  flux_df <- interpolate_con_to_q(chem_df = chem_df, 
                                  q_df = q_df) %>%
   mutate(flux = (con*q_lps*86400)/(ws_size*1e6),
          method = 'bear') %>%
   select(date, flux, method)
   
  return(flux_df)
  
}