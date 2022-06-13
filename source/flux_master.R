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
library(RiverLoad)
library(lfstat)
library(imputeTS)

# thinning_intervals <- c('monthly', 'bi_weekely', 'weekely')
thinning_intervals <- c('daily', 'weekly', 'biweekly', 'monthly', 'quarterly')
vars <- c('nitrate_nitrite_mgl', 'spcond_uscm')

# Read in site data 
# (Sites were identified with the identify_usgs_gauges.R script)
site_var_data <- read.csv('data/general/site_var_data.csv', colClasses = 'character') %>%
    filter(parm_cd %in% c('00095', '99133')) %>%
    distinct(site_code, parm_cd, .keep_all = T)

n_sites <- site_var_data %>%
    filter(parm_cd == '99133')

# Only want sites with continuous N 
site_var_data <- site_var_data %>%
    filter(site_code %in% !!n_sites$site_code)


site_data <- read.csv('data/general/site_data.csv', colClasses = 'character')
variable_data <- read.csv('data/general/variable_data.csv', colClasses = 'character')


#### download raw data #### 

# get site areas, lat, and long
areas <- dataRetrieval::readNWISsite(site_data$site_code) %>%
    select(site_no, drain_area_va, dec_lat_va, dec_long_va) %>%
    mutate(ws_area_ha = drain_area_va*259) %>%
    select(site_code = site_no, ws_area_ha, lat = dec_lat_va, long = dec_long_va)
site_data <- left_join(site_data, areas)
write_csv(site_data, 'data/general/site_data.csv')

# get site ecogeographic region

# Variables
failed_sites <- c()
for(i in 1:nrow(site_var_data)){
    
    # Load in site info 
    site_info <- site_var_data[i,]
    site_code <- site_info$site_code
    parm_cd <- site_info$parm_cd
    
    var <- variable_data %>%
        filter(usgs_parm_cd == !!parm_cd) %>%
        pull(var)
    
    # Download variable  
    var_data <- readNWISuv(site_code, parameterCd = parm_cd)
    
    if(nrow(var_data) == 0) next
    
    if(!paste0('X_', parm_cd, '_00000') %in% names(var_data)) {
        print('weird var found')
        failed_sites <- c(failed_sites, site_code)
        next
    }
    
    var_data <- var_data %>%
        select(site_code = site_no, datetime = dateTime, val = !!paste0('X_', parm_cd, '_00000')) %>%
        mutate(val = na.approx(val, rule = 2)) %>%
        mutate(var = !!var) %>%
        filter(val >= 0)
    
    directory <- glue('data/raw/{v}',
                      v = var)
    
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    
    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)
    
    write_feather(var_data, file_path)
    
}

# Discharge 
sites <- unique(site_var_data$site_code)
sites <- sites[!sites %in% failed_sites]
for(i in 1:length(sites)){
    # Download Q data 
    site_code <- sites[i]
    q_data <- readNWISuv(site_code, parameterCd = '00060')
    
    # If there is no unit data check for daily values 
    if(nrow(q_data) == 0) {
        q_data <- readNWISdv(site_code, parameterCd = '00060')
        date_name <- 'Date'
    } else{
        date_name <- 'dateTime'
    }
    
    if(nrow(q_data) == 0) {
        print('no Q found')
        next
    }
    
    var_names <- names(q_data)
    var_names <- grep('X_00060', var_names, value = T)
    var_names <- var_names[!grepl('_cd', var_names)]
    
    q_data <- q_data %>%
        select(site_code = site_no, datetime = !!date_name, val = !!var_names) %>%
        mutate(val = na.approx(val, rule = 2)) %>%
        mutate(var = 'q_cfs')
    
    q_directory <- glue('data/raw/q_cfs')
    
    if(!dir.exists(q_directory)){
        dir.create(q_directory, recursive = TRUE)
    }
    
    file_path <- glue('{q_directory}/{s}.feather',
                      s = site_code)
    
    write_feather(q_data, file_path)
    
}

# skip to here if rerunning
# get sites that have both Q and chem 
# q_sites <- list.files('data/raw/q_cfs/')
# chem_sites <- list.files('data/raw/nitrate_nitrite_mgl/')
# spcond_sites <- list.files('data/raw/spcond_uscm/')
# common_sites <- q_sites[q_sites %in% chem_sites]
# common_sites <- common_sites[common_sites %in% spcond_sites]
# common_sites <- str_split_fixed(common_sites, '\\.', n = Inf)[,1]

# Filter out sites that failed to download for all parameters (Q, nitrate, and spcond)
# site_var_data <- site_var_data %>%
#     filter(site_code %in% !!common_sites)

#### Thin data to desired intervals ####
# You have more built out code for this from what I remember but I am just 
# doing a simple thinning here as an example 
site_var_data <- site_var_data %>%
    filter(site_code %in% !!sites)


# get site flashiness index
site_RBI_df <- data.frame()

for(i in 1:nrow(site_var_data)) {

    site <- site_var_data$site_code[i]
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site}.feather')))

    if(inherits(q_data, 'try-error')) next

    parm_cd <- site_var_data$parm_cd[i]
    var <- variable_data %>%
        filter(usgs_parm_cd == !!parm_cd) %>%
      pull(var)

    site_RBI <- ContDataQC::RBIcalc(q_data$val)
    output <- c(site, site_RBI)
    site_RBI_df <- rbind(site_RBI_df, output)
}

write.csv(site_RBI_df, "data/general/site_RBI_data.csv")

for(i in 1:nrow(site_var_data)) {
    
    site <- site_var_data$site_code[i]
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site}.feather')))
    
    if(inherits(q_data, 'try-error')) next
    
    parm_cd <- site_var_data$parm_cd[i]
    var <- variable_data %>%
        filter(usgs_parm_cd == !!parm_cd) %>%
        pull(var)
    
    chem_data <- try(read_feather(glue('data/raw/{var}/{site}.feather')))
    
    if(inherits(chem_data, 'try-error')) next
    
    # Loop through thinning intervals (can add if statements for each of the thinning 
    # intervals)
    ## p <- 5
    for(p in 1:length(thinning_intervals)){
        
        if(thinning_intervals[p] == 'daily'){
            chem_data_thin <- chem_data %>%
                # Idealy we would not want times to be selected that are in the 
                # middle of the night, this could be improved by including local time
                filter(hour(datetime) %in% c(13:18)) %>%
                mutate(date = lubridate::date(datetime)) %>%
                distinct(date, .keep_all = T) %>%
                select(-datetime)
            
            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])
            
            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }
            
            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))
        }
        
        if(thinning_intervals[p] == 'weekly'){
            chem_data_thin <- chem_data %>%
                # mutate(date = date(datetime)) %>%
                # group_by(site_code, var, date) %>%
                # summarise(val = mean(val, na.rm = T)) %>%
                filter(hour(datetime) %in% c(13:18)) %>%
                filter(lubridate::wday(datetime) == 1) %>%
                mutate(date = lubridate::date(datetime)) %>%
                distinct(date, .keep_all = T) %>%
                select(-datetime)
            
            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])
            
            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }
            
            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))
        }
        
        if(thinning_intervals[p] == 'biweekly'){
            chem_data_thin <- chem_data %>%
                # mutate(date = date(datetime)) %>%
                # group_by(site_code, var, date) %>%
                # summarise(val = mean(val, na.rm = T)) %>%
                filter(hour(datetime) %in% c(13:18)) %>%
                filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
                mutate(date = lubridate::date(datetime)) %>%
                distinct(date, .keep_all = T) %>%
                select(-datetime)
            
            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])
            
            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }
            
            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))
        }
        
        if(thinning_intervals[p] == 'monthly'){
            
            chem_data_thin <- chem_data %>%
                filter(hour(datetime) %in% c(13:18)) %>%
                mutate(date = date(datetime)) %>%
                filter(day(date) == 1)
            
            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])
            
            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }
            
            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))
            
        }

        if(thinning_intervals[p] == 'quarterly'){

            chem_data_thin <- chem_data %>%
                filter(hour(datetime) %in% c(13:18)) %>%
                mutate(quarter = quarters(datetime)) %>%
                group_by(quarter) %>%
                top_n(1, datetime)

            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])

            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }

            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))

        }
    }
}

#load('20220608_environment.RData')
#### Calculate fluxes with various methods ####
source('source/helper_functions.R')
source('source/flux_method_egret_daily.R')
source('source/flux_method_hbef_annual.R')
source('source/flux_method_hbef_daily.R')
source('source/flux_method_fernow_annual.R')
source('source/flux_method_fernow_weekly.R')
source('source/flux_method_bear_annual.R')
source('source/flux_method_bear.R')
source('source/flux_method_santee_annual.R')
source('source/flux_method_santee.R')
source('source/flux_method_beale_monthly.R')
source('source/flux_method_beale_annual.R')
source('source/flux_method_rating_daily.R')
source('source/flux_method_rating_annual.R')

for(i in 10:nrow(site_var_data)){
    
    site_code <- site_var_data[i,1]
    parm_cd <- site_var_data[i,3]
    var <- variable_data %>%
        filter(usgs_parm_cd == !!parm_cd) %>%
        pull(var)
    
    site_info <- site_data %>%
        filter(site_code == !!site_code) 
    
    area <- site_info %>%
        pull(ws_area_ha) %>%
        as.numeric()
    
    lat <- site_info %>%
        pull(lat) %>%
        as.numeric()
    
    long <- site_info %>%
        pull(long) %>%
        as.numeric()
    
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site_code}.feather')))
    
    if(inherits(q_data, 'try-error') || nrow(q_data) == 0) next
    
    q_check <- q_data %>%
        mutate(date = date(datetime)) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(datetime, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)
    
    conc_data <- try(read_feather(glue('data/thinned/{var}/daily/{site_code}.feather')))
    if(inherits(conc_data, 'try-error') || nrow(conc_data) == 0) next
    
    conc_check <- conc_data %>%
        mutate(water_year = water_year(date, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)
    
    q_good_years <- q_check$water_year
    conc_good_years <- conc_check$water_year
    
    if(length(q_good_years[q_good_years %in% conc_good_years]) == 0) {
        print('No good years')
        next
    } else{
        good_years <- q_good_years[q_good_years %in% conc_good_years]
    }
    
    conc_data_raw <- try(read_feather(glue('data/raw/{var}/{site_code}.feather')))
    if(inherits(conc_data_raw, 'try-error') || nrow(conc_data_raw) == 0) next
    
    q_data_raw <- try(read_feather(glue('data/raw/q_cfs/{site_code}.feather')))
    if(inherits(q_data_raw, 'try-error') || nrow(q_data_raw) == 0) next
    
    conc_data_raw <- conc_data_raw %>%
        rename(!!var := val) %>%
        select(datetime, !!var)
    
    q_data_raw <- q_data_raw %>%
        rename(Q = val) %>%
        select(datetime, Q)
    
    # These units are kg/ha/day
    true_flux_daily <- full_join(q_data_raw, conc_data_raw) %>%
        mutate(Q = Q*28.316847) %>%
        mutate(flux = ((Q*.data[[var]])/1000000)*60*15) %>%
        mutate(wy = water_year(datetime, origin = "usgs")) %>%
        filter(wy %in% !!good_years) 
    
    # These units are kg/ha/yr
    true_flux_annual <- true_flux_daily %>%
        group_by(wy) %>%
        summarise(flux_annual_kg_ha = sum(flux, na.rm = TRUE)) %>%
        mutate(flux_annual_kg_ha = flux_annual_kg_ha/area) %>%
        mutate(method = 'real')
    
    if(!dir.exists(glue('data/fluxes/daily/true/{var}'))){
        dir.create(glue('data/fluxes/daily/true/{var}'), recursive = TRUE)
    }
    
    write_feather(true_flux_daily, glue('data/fluxes/daily/true/{var}/{site_code}.feather'))
    
    if(!dir.exists(glue('data/fluxes/annual/true/{var}'))){
        dir.create(glue('data/fluxes/annual/true/{var}'), recursive = TRUE)
    }
    
    write_feather(true_flux_annual, glue('data/fluxes/annual/true/{var}/{site_code}.feather'))
    
    for(t in 1:length(thinning_intervals)){
        
        thinning_interval <- thinning_intervals[t]
        # Prep data for functions 
        conc_data <- try(read_feather(glue('data/thinned/{var}/{thinning_interval}/{site_code}.feather')))
        
        if(inherits(conc_data, 'try-error') || nrow(conc_data) == 0) next
        
        conc_data_prep <- conc_data %>%
            select(date, con = val)
        
        q_data_prep <- q_data %>%
            mutate(date = date(datetime)) %>%
            group_by(site_code, date, var) %>%
            summarise(val = mean(val, na.rm = T)) %>%
            mutate(q_lps = val*28.316847) %>%
            select(date, q_lps)
        
        # Prep data 
        conc_data_prep <- conc_data_prep %>%
            mutate(water_year = water_year(date, origin = "usgs")) %>%
            filter(water_year %in% !!good_years)
        
        q_data_prep <- q_data_prep %>%
            mutate(water_year = water_year(date, origin = "usgs")) %>%
            filter(water_year %in% !!good_years)
        
        ## Flux Methods
        
        # HBEF
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/hbef/')
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
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/hbef/')
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
        directory_raw <- glue('data/fluxes/daily/{thinning_interval}/{var}/egret/')
        if(!dir.exists(directory_raw)){
            dir.create(directory_raw, recursive = TRUE)
        }
        directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/egret/')
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
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/fernow/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        fernow_flux_daily <- try(estimate_flux_fernow_weekly(chem_df = conc_data_prep, 
                                                             q_df = q_data_prep, 
                                                             ws_size = area) %>%
            filter(wy %in% !!good_years))
        
        if(inherits(fernow_flux_daily, 'try-error')){
            
        } else{
            write_feather(fernow_flux_daily, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/fernow/')
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
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/bear/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        conv_q <- q_data_raw %>%
            mutate(site_code = !!site_code,
                   q_lps = Q*28.316847) %>%
            select(site_code, date = datetime, q_lps)
        
        bear_flux_daily <- try(estimate_flux_bear(chem_df = conc_data_prep, 
                                                   q_df = conv_q, 
                                                   ws_size = area) %>%
            mutate(wy = water_year(date, origin = 'usgs')) %>%
            filter(wy %in% !!good_years))
        
        if(inherits(bear_flux_daily, 'try-error')){
            
        } else{
            write_feather(bear_flux_daily, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/bear/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        
        bear_flux_annual <- try(estimate_flux_bear_annual(chem_df = conc_data_prep, 
                                                          q_df = conv_q, 
                                                          ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(bear_flux_annual, 'try-error')){
        } else{
            write_feather(bear_flux_annual, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # Santee
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/santee/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        daily_santee_flux <- try(estimate_flux_santee(chem_df = conc_data_prep, 
                                                     q_df = q_data_prep, 
                                                     ws_size = area) %>%
                                     filter(wy %in% !!good_years))
        
        if(inherits(daily_santee_flux, 'try-error')){
        } else{
            write_feather(daily_santee_flux, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/santee/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        annual_santee_flux <- try(estimate_flux_santee_annual(chem_df = conc_data_prep, 
                                                              q_df = q_data_prep, 
                                                              ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(annual_santee_flux, 'try-error')){
        } else{
            write_feather(annual_santee_flux, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # Beale
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/beale/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        monthly_beale_flux <- try(estimate_flux_beale_monthly(chem_df = conc_data_prep,
                                     q_df = q_data_raw,
                                     ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(monthly_beale_flux, 'try-error')){
        } else{
            write_feather(monthly_beale_flux, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/beale/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        annual_beale_flux <- try(estimate_flux_beale_annual(chem_df = conc_data_prep, 
                                                              q_df = q_data_raw, 
                                                              ws_size = area) %>%
                                     filter(wy %in% !!good_years)) 
            
        
        if(inherits(annual_beale_flux, 'try-error')){
        } else{
            write_feather(annual_beale_flux, glue('{annual_directory}/{site_code}.feather'))
        }
        
        # Rating 
        daily_directory <- glue('data/fluxes/daily/{thinning_interval}/{var}/rating/')
        if(!dir.exists(daily_directory)){
            dir.create(daily_directory, recursive = TRUE)
        }
        
        daily_rating_flux <- try(estimate_flux_rating_daily(chem_df = conc_data_prep,
                                                              q_df = q_data_prep,
                                                              ws_size = area)) %>%
            filter(wy %in% !!good_years)
        
        if(inherits(daily_rating_flux, 'try-error')){
        } else{
            write_feather(daily_rating_flux, glue('{daily_directory}/{site_code}.feather'))
        }
        
        annual_directory <- glue('data/fluxes/annual/{thinning_interval}/{var}/rating/')
        if(!dir.exists(annual_directory)){
            dir.create(annual_directory, recursive = TRUE)
        }
        annual_rating_flux <- try(estimate_flux_rating_annual(chem_df = conc_data_prep, 
                                                            q_df = q_data_raw, 
                                                            ws_size = area) %>%
                                     filter(wy %in% !!good_years)) 
        
        
        if(inherits(annual_beale_flux, 'try-error')){
        } else{
            write_feather(annual_beale_flux, glue('{annual_directory}/{site_code}.feather'))
        }
        
        
    }
    
    print(paste(site_code, ' done!'))
}
 


