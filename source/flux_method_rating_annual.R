estimate_flux_rating_annual <- function(chem_df, q_df, ws_size){
    
    annual <- estimate_flux_rating_daily(chem_df = conc_data_prep,
                               q_df = q_data_prep,
                               ws_size = area) %>%
        group_by(wy) %>%
        summarize(flux = sum(flux)) %>%
        mutate(method = 'rating')
    
    return(annual)
}