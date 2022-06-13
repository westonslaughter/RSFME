#### Set up ####
# Load in packages 
library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(glue)
library(feather)
library(zoo)
library(here)
library(lfstat)
library(macrosheds)

site_data <- ms_download_site_data()

Mode <- function(x, na.rm = TRUE){
    
    if(na.rm){
        x <- na.omit(x)
    }
    
    ux <- unique(x)
    mode_out <- ux[which.max(tabulate(match(x, ux)))]
    return(mode_out)
    
}

discharge <- ms_load_product(macrosheds_root = 'data/ms_data/',
                prodname = 'discharge',
                warn = F)

stream_gauges <- unique(discharge$site_code)

vars <- list(NO3_N = c('NO3_N', 'NO3_NO2_N'), Na = 'Na', Cl = 'Cl', DOC = 'DOC')
for(v in 1:length(vars)){
    
    var <- names(vars[v])
    filter_vars <- vars[v][[1]]
    
    chemistry <- ms_load_product(macrosheds_root = 'data/ms_data/',
                                 prodname = 'stream_chemistry',
                                 filter_vars = filter_vars,
                                 warn = F) %>%
        mutate(var = ms_drop_var_prefix(var))
    
    chemistry <- chemistry %>% 
        filter(site_code %in% !!stream_gauges)
    
    #### Get good site years ####
    sampling_intervals <- tibble()
    for(i in 1:length(stream_gauges)){
        look <- chemistry %>%
            filter(site_code == !!stream_gauges[i],
                   ms_interp == 0) %>%
            filter(!is.na(val)) %>%
            arrange(datetime) 
        
        t_dif <- diff(look$datetime)
        
        look$t_dif <- c(NA, t_dif)
        
        sampling_dif <- Mode(look$t_dif)
        
        fin_tib <- tibble(site_code = stream_gauges[i],
                          interval = sampling_dif)
        
        sampling_intervals <- rbind(sampling_intervals, fin_tib)
    }
    
    good_chem_sites <- sampling_intervals %>%
        filter(interval < 20)
    
    good_chem_site_years <- tibble()
    for(i in 1:nrow(good_chem_sites)){
        
        site_code <- good_chem_sites[[i,1]]
        interval <-  good_chem_sites[[i,2]]
        
        chem_site <- chemistry %>%
            filter(site_code == !!site_code)
        
        if(nrow(chem_site) == 0) next
        
        chem_site <- chem_site %>% 
            filter(ms_interp == 0,
                   !is.na(val))
        
        if(length(unique(chem_site$var)) != 1){
            chem_with_most <- chem_site %>% 
                group_by(var) %>% 
                summarise(n = n())
            
            var_we_want <- chem_with_most %>%
                filter(n == max(chem_with_most$n)) %>%
                pull(var)
            
            chem_site <- chem_site %>%
                filter(var == !!var_we_want)
        } else{
            var_we_want <- unique(chem_site$var)
        }
        
        if(interval <= 8){
            sample_filter <- 45
        } else{
            sample_filter <- 23
        }
        
        chem_site <- chem_site %>%
            mutate(wy = water_year(datetime, origin = 'usgs')) %>%
            group_by(wy) %>%
            summarise(n = n()) %>%
            filter(n >= !!sample_filter) %>%
            mutate(site_code = !!site_code) %>%
            mutate(var = !!var_we_want)
        
        good_chem_site_years <- rbind(good_chem_site_years, chem_site)
    }
    
    good_chem_site_years <- good_chem_site_years %>%
        rename(chem_n = n)
    
    potential_sites <- unique(good_chem_site_years$site_code)
    
    good_site_years <- tibble()
    # Get good discharge site years 
    for(i in 1:length(potential_sites)) {
        
        site_code <- potential_sites[i]
        
        site_q <- discharge %>%
            filter(site_code == !!site_code, 
                   !is.na(val)) %>%
            mutate(wy = water_year(datetime, origin = 'usgs')) %>%
            group_by(wy) %>%
            summarise(n = n()) %>%
            filter(n >= 320)
        
        if(nrow(site_q) == 0) {
            next
        } else{
            site_q <- site_q %>%
                mutate(site_code = !!site_code) %>%
                rename(q_n = n)
            
            good_chem_site_years_site <- good_chem_site_years %>%
                filter(site_code == !!site_code)
            
            good_fin <- inner_join(site_q, good_chem_site_years_site, by = c('wy', 'site_code'))    
            
        }
        good_site_years <- rbind(good_site_years, good_fin)
    }
    
    good_site_years <- good_site_years %>%
        arrange(site_code, wy, .keep_all = T)
    
    good_sites <- unique(good_site_years$site_code)
    
    # Calculate fluxes
    source('source/flux_method_egret_daily.R')
    source('source/flux_method_hbef_annual.R')
    source('source/flux_method_hbef_daily.R')
    source('source/flux_method_fernow_annual.R')
    source('source/flux_method_fernow_weekly.R')
    source('source/flux_method_bear_annual.R')
    source('source/flux_method_bear_hourly.R')
    source('source/flux_method_santee_annual.R')
    source('source/flux_method_santee.R')
    
    for(i in 1:length(good_sites)){
        site_code <- good_sites[i]
        area <- site_data %>%
            filter(site_code == !!site_code,
                   site_type == 'stream_gauge') %>%
            pull(ws_area_ha)
        
        lat <- site_data %>%
            filter(site_code == !!site_code,
                   site_type == 'stream_gauge') %>%
            pull(latitude)
        
        long <- site_data %>%
            filter(site_code == !!site_code,
                   site_type == 'stream_gauge') %>%
            pull(longitude)
        
        
        chem_var <- good_site_years %>%
            filter(site_code == !!site_code) %>%
            pull(var) %>%
            unique()
        
        good_site_years_ <- good_site_years %>%
            filter(site_code == !!site_code)
        
        good_years <- good_site_years_$wy
        
        site_q <- discharge %>%
            filter(site_code == !!site_code) %>%
            mutate(wy = water_year(datetime, origin = 'usgs')) %>%
            filter(wy %in% !!good_site_years_$wy)
        
        site_chem <- chemistry %>%
            filter(site_code == !!site_code) %>%
            mutate(wy = water_year(datetime, origin = 'usgs')) %>%
            filter(wy %in% !!good_site_years_$wy) %>%
            filter(var == !!chem_var)
        
        conc_data_prep <- site_chem %>%
            select(date = datetime, con = val)
        
        q_data_prep <- site_q %>%
            group_by(site_code, datetime, var) %>%
            summarise(val = mean(val, na.rm = T)) %>%
            ungroup() %>%
            select(date = datetime, q_lps = val)
        
        # Methods 
        # HBEF
        daily_directory <- glue('data/ms_fluxes/daily/{var}/hbef/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        hbef_flux_daily <- try(estimate_flux_hbef_daily(chem_df = conc_data_prep, 
                                                        q_df = q_data_prep, 
                                                        ws_size = area) %>%
                                   mutate(wy = water_year(date, origin = 'usgs')) %>%
                                   filter(wy %in% !!good_years))
        
        if(inherits(hbef_flux_daily, 'try-error')){
        } else{
            write_feather(hbef_flux_daily, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/ms_fluxes/annual/{var}/hbef/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        
        hbef_flux_annual <- try(estimate_flux_hbef_annual(chem_df = conc_data_prep, 
                                                          q_df = q_data_prep, 
                                                          ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(hbef_flux_annual, 'try-error')){
        } else{
            write_feather(hbef_flux_annual, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # EGRET (WRTDS) 
        directory_raw <- glue('data/ms_fluxes/daily/{var}/egret/')
        if(!dir.exists(directory_raw)){
            dir.create(directory_raw, recursive = TRUE)
        }
        directory <- glue('data/ms_fluxes/annual/{var}/egret/')
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
        }
        
        egret_flux <- try(adapt_ms_egret(chem_df = conc_data_prep,
                                         q_df = q_data_prep,
                                         ws_size = area,
                                         lat = lat,
                                         long = long))
        
        if(inherits(egret_flux, 'try-error')){
            
        } else{
            write_rds(egret_flux, glue('{directory_raw}/{site_code}.rds'))
            
            # Get annual flux 
            egret_annual <- egret_flux$Daily %>%
                mutate(wy = water_year(Date, origin = 'usgs')) %>%
                filter(wy %in% !!good_years) %>%
                group_by(wy) %>%
                summarise(flux_annual_kg_ha = sum(FluxDay, na.rm = TRUE)) %>%
                mutate(flux_annual_kg_ha = (flux_annual_kg_ha/!!area)) %>%
                mutate(method = 'wrtds')
            
            write_feather(egret_annual, glue('{directory}/{site_code}.feather'))
            
        }
        
        
        # Fernow 
        daily_directory <- glue('data/ms_fluxes/daily/{var}/fernow/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        fernow_flux_daily <- try(estimate_flux_fernow_weekly(chem_df = conc_data_prep, 
                                                             q_df = q_data_prep, 
                                                             ws_size = area) %>%
                                     mutate(wy = water_year(date, origin = 'usgs')) %>%
                                     filter(wy %in% !!good_years))
        
        if(inherits(fernow_flux_daily, 'try-error')){
            
        } else{
            write_feather(fernow_flux_daily, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/ms_fluxes/annual/{var}/fernow/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        
        fernow_flux_annual <- try(estimate_flux_fernow_annual(chem_df = conc_data_prep, 
                                                              q_df = q_data_prep, 
                                                              ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(fernow_flux_annual, 'try-error')){
            
        } else{
            write_feather(fernow_flux_annual, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # Bear 
        daily_directory <- glue('data/ms_fluxes/daily/{var}/bear/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        q_data_hourly <- site_q %>%
            ms_synchronize_timestep(., desired_interval = '15 min', impute_limit = Inf) %>%
            select(site_code, date = datetime, q_lps = val)
        
        bear_flux_daily <- try(estimate_flux_bear_hourly(chem_df = conc_data_prep, 
                                                         q_df = q_data_hourly, 
                                                         ws_size = area) %>%
                                   mutate(wy = water_year(date, origin = 'usgs')) %>%
                                   filter(wy %in% !!good_years))
        
        if(inherits(bear_flux_daily, 'try-error')){
            
        } else{
            write_feather(bear_flux_daily, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/ms_fluxes/annual/{var}/bear/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        
        bear_flux_annual <- try(estimate_flux_bear_annual(chem_df = conc_data_prep, 
                                                          q_df = q_data_hourly, 
                                                          ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(bear_flux_annual, 'try-error')){
        } else{
            write_feather(bear_flux_annual, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # Santee
        daily_directory <- glue('data/ms_fluxes/daily/{var}/santee/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        daily_santee_flux <- try(estimate_flux_santee(chem_df = conc_data_prep, 
                                                      q_df = q_data_hourly, 
                                                      ws_size = area) %>%
                                     mutate(wy = water_year(date, origin = 'usgs')) %>%
                                     filter(wy %in% !!good_years))
        
        if(inherits(daily_santee_flux, 'try-error')){
        } else{
            write_feather(daily_santee_flux, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/ms_fluxes/annual/{var}/santee/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        annual_santee_flux <- try(estimate_flux_santee_annual(chem_df = conc_data_prep, 
                                                              q_df = q_data_hourly, 
                                                              ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(annual_santee_flux, 'try-error')){
        } else{
            write_feather(annual_santee_flux, glue('{annual_directory}/{site_code}.feather'))
        }
        
    }
    
}


# site_data %>%
#     filter(site_code %in% !!good_sites) %>%
#     sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>%
#     mapview::mapview()
