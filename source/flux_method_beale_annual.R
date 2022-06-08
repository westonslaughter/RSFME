estimate_flux_beale_annual <- function(q_df, chem_df, ws_size){
    area = ws_size
    
    db <- prep_usgs_for_riverload(q_df = q_df, chem_df = chem_df)
    
    out <- beale.ratio(db, 1, period = 'year') %>%
        tibble()
    out$year <- dimnames(out$.)[[1]]
    colnames(out) <- c('con', 'year')
    conv_out <- out %>%
        mutate(con = con/(1000*area))
    
    out_formatted <- tibble(wy = conv_out$year,
                            flux = as.numeric(conv_out$con))
    
    return(out_formatted)
}