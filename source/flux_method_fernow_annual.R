# create function of fernow flux estimation for annual
estimate_flux_fernow_annual <- function(chem_df, q_df, ws_size){
  startDate <- min(chem_df$date)
  endDate <- max(chem_df$date)
  source('source/flux_method_fernow_weekly.R')
  flux_df <- estimate_flux_fernow_weekly(chem_df = chem_df, q_df = q_df, ws_size = ws_size)
  # compute total flux in kg/ha and return it
  out <- flux_df %>%
    mutate(wy = ifelse(month(date) %in% c(10, 11, 12), year(date) + 1, year(date))) %>%
    group_by(wy) %>%
    na.omit() %>%
    filter(n() >= 45) %>%
    summarize(flux_annual_kg_ha = sum(flux_weekly_kg_ha)) %>%
    mutate(method = 'fernow')
  return(out)
}
