# create function of hbef flux estimation for annual
estimate_flux_santee_annual <- function(chem_df, q_df, ws_size){
  source('source/flux_method_santee.R')
  flux_df <- estimate_flux_santee(chem_df = chem_df, q_df = q_df, ws_size = ws_size)
  # compute total flux in kg/ha and return it
  out <- flux_df %>%
    group_by(wy = water_year(date, origin = "usgs")) %>%
    #filter(n() > 360) %>%
    summarize(flux_annual_kg_ha = sum(flux_kg_ha))
  
  out <- out %>%
    mutate(method = 'santee')
  return(out)
}