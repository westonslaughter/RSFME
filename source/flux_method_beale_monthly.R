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
    
    #make daily to match q assumptions
    daily_out <- tibble(date = as.numeric(), flux = as.numeric())
    for(i in 1:nrow(out_formatted)){
        monthDays = as.numeric(days_in_month(out_formatted$datetime[i]))
        
        dailyDates = seq(out_formatted$datetime[i], by = 'days', length.out = monthDays)
        
        agg_flux <- flux_df$flux[i]
        
        chunk_out <- tibble(date = dailyDates, flux = agg_flux/monthDays)
        
        daily_out <- rbind(daily_out, chunk_out)
        
    }
    daily_out <- daily_out %>%
        mutate(method = 'beale', wy = water_year(date, origin = 'usgs'))   
    
    return(daily_out)
}
