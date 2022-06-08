estimate_flux_beale_monthly <- function(q_df, chem_df, ws_size){
    area = ws_size

    db <- prep_usgs_for_riverload(q_df = q_df, chem_df = chem_df)
    
    out <- beale.ratio(db, 1, period = 'month') %>%
        tibble()
    out$month <- dimnames(out$.)[[1]]
    colnames(out) <- c('con', 'month')
    conv_out <- out %>%
            mutate(con = con/(1000*area))
    
    out_formatted <- tibble(datetime = ym(conv_out$month),
        flux = as.numeric(conv_out$con)) %>%
        mutate(wy = water_year(datetime))
    
    return(out_formatted)
}
