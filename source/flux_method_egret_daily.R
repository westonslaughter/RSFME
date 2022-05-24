
adapt_ms_egret <- function(chem_df, q_df, ws_size, kalman = FALSE){
  
  ms_chem <- chem_df %>%
    mutate(date = date(date)) %>%
    group_by(date) %>%
    summarise(nitrate_mgL = mean(nitrate_mgL, na.rm = T)) %>%
    mutate(site_code = 'none',
           var = 'IS_NO3',
           ms_status = 0,
           ms_interp = 0) %>%
    rename(val = nitrate_mgL,
           datetime = date)
    
    ms_q <- q_df %>%
      mutate(date = date(date)) %>%
      group_by(date) %>%
      summarise(q_cfs = mean(q_cfs, na.rm = T)) %>%
      mutate(site_code = 'none',
             var = 'IS_discharge',
             ms_status = 0,
             ms_interp = 0) %>%
      rename(val = q_cfs,
             datetime = date)
    
    site_data <- tibble(site_code = 'none',
                        ws_area_ha = 47655.8,
                        latitude = 35.461389,
                        longitude = -83.353611)
    
    egret_results <- macrosheds::ms_run_egret(stream_chemistry = ms_chem, discharge = ms_q,
                                              prep_data = FALSE, site_data = site_data,
                                              kalman = kalman)
    
    return(egret_results)
}

# run_output <- adapt_ms_egret(chem_df = ocon_chem, q_df = ocon_q, ws_size = ocon_ws_area_ha)
# EGRET::plotConcHist(run_output)
# EGRET::plotConcTime(run_output)
# EGRET::plotFluxQ(run_output)
# EGRET::plotConcPred(run_output)
