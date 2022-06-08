# simple linear interpolation between con samples at frequency of q
# in:
# out: 
interpolate_con_to_q <- function(chem_df, q_df){
 # Join chem and q
 join_df <- full_join(chem_df, q_df, by = 'date') %>%
   arrange(date) %>%
   mutate(interp = 0) %>%
    select(date, con, q_lps, interp)
 join_df$interp[is.na(join_df$con)] <- 1
 
 #interpolate chem to match q
 interp_df <- na_interpolation(join_df, option = "linear", maxgap = Inf) %>%
     mutate(wy = water_year(date, origin = 'usgs'))
 return(interp_df)
}

# take usgs data and make it ready for riverload functions
# in: 
# out:
prep_usgs_for_riverload <- function(chem_df, q_df){
    conv_q <- q_df %>%
        mutate(site_code = !!site_code,
               flow = Q/35.3147) %>% # convert cfs to cubic meters per second)
        select(datetime, flow) %>%
        data.frame()
    
    conv_c <- chem_df %>%
        mutate(datetime = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
        select(datetime, con) %>%
        data.frame()
    
    db <- full_join(conv_q, conv_c, by = "datetime") %>%
        #filter(!is.na(flow)) %>%
        arrange(datetime)
    
    return(db)
}

