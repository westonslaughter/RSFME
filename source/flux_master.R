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

thinning_intervals <- c('monthly', 'bi_weekely', 'weekely')
vars <- c('nitrate_nitrite_mgl', 'spcond_uscm')

#### download raw data #### 
# Read in site data 
# (Sites were identified with the identify_usgs_gauges.R script)
site_var_data <- read.csv('data/general/site_var_data.csv', colClasses = 'character') %>%
    filter(parm_cd %in% c('00095', '99133'))
n_sites <- site_var_data %>%
    filter(parm_cd == '99133')

site_var_data <- site_var_data %>%
    filter(site_code %in% !!n_sites$site_code)


site_data <- read.csv('data/general/site_data.csv', colClasses = 'character')
variable_data <- read.csv('data/general/variable_data.csv', colClasses = 'character')

# get site areas, lat, and long
# areas <- dataRetrieval::readNWISsite(site_data$site_code) %>%
#     select(site_no, drain_area_va, dec_lat_va, dec_long_va) %>%
#     mutate(ws_area_ha = drain_area_va*259) %>%
#     select(site_code = site_no, ws_area_ha, lat = dec_lat_va, long = dec_long_va)
# site_data <- left_join(site_data, areas)
# write_csv(site_data, 'data/general/site_data.csv')

# Variables
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
        next
    }
    
    var_data <- var_data %>%
        select(site_code = site_no, datetime = dateTime, val = !!paste0('X_', parm_cd, '_00000')) %>%
        mutate(val = na.approx(val, rule = 2)) %>%
        mutate(var = !!var)
    
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
for(i in 1:length(sites)){
    # Download Q data 
    site_code <- sites[i]
    q_data <- readNWISuv(site_code, parameterCd = '00060')
    
    # Incase a site has a varibale but not Q
    if(nrow(q_data) == 0) next
    
    q_data <- q_data %>%
        select(site_code = site_no, datetime = dateTime, val = X_00060_00000) %>%
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

# get sites that have both Q and chem 
q_sites <- list.files('data/raw/q_cfs/')
chem_sites <- list.files('data/raw/nitrate_nitrite_mgl/')
spcond_sites <- list.files('data/raw/spcond_uscm/')
common_sites <- q_sites[q_sites %in% chem_sites]
common_sites <- common_sites[common_sites %in% spcond_sites]
common_sites <- str_split_fixed(common_sites, '\\.', n = Inf)[,1] 

# Filter out sites that failed to download for all parameters (Q, nitrate, and spcond)
site_var_data <- site_var_data %>%
    filter(site_code %in% !!common_sites)


#### Thin data to desired intervals ####
# You have more built out code for this from what I remember but I am just 
# doing a simple thinning here as an example 

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
    for(p in 1:length(thinning_intervals)){
        
        if(thinning_intervals[p] == 'monthly'){
            
            chem_data_thin <- chem_data %>%
                mutate(date = date(datetime)) %>%
                group_by(site_code, var, date) %>%
                summarise(val = mean(val, na.rm = T)) %>%
                filter(day(date) == 1)
            
            directory <- glue('data/thinned/{var}/{t}',
                              t = thinning_intervals[p])
            
            if(!dir.exists(directory)){
                dir.create(directory, recursive = TRUE)
            }
            
            write_feather(chem_data_thin, glue('{directory}/{site}.feather'))
            
        }
        
        if(thinning_intervals[p] == 'bi_weekely'){
            # add code here 
        }
        
        if(thinning_intervals[p] == 'weekely'){
            # add code here 
    }
    }
}


#### Calculate flues with various methods ####

for(i in 50:nrow(site_var_data)){
    
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
    
    for(t in 1:length(thinning_intervals)){
        
        thinning_interval <- thinning_intervals[t]
        # Prep data for functions 
        conc_data <- try(read_feather(glue('data/thinned/{var}/{thinning_interval}/{site_code}.feather')))
        
        if(inherits(conc_data, 'try-error') || nrow(conc_data) == 0) next
        
        conc_data_prep <- conc_data %>%
            select(date, con = val)
        
        q_data <- try(read_feather(glue('data/raw/q_cfs/{site_code}.feather')))
        
        if(inherits(q_data, 'try-error') || nrow(q_data) == 0) next
        
        q_data_prep <- q_data %>%
            mutate(date = date(datetime)) %>%
            group_by(site_code, date, var) %>%
            summarise(val = mean(val, na.rm = T)) %>%
            mutate(q_lps = val*28.316847) %>%
            select(date, q_lps)
        
        # Add methods here 
        
        # HBEF
        directory <- glue('data/fluxes/{thinning_interval}/{var}/hbef/')
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
        }
        
        source('source/flux_method_hbef_daily.R')
        hbef_flux <- try(estimate_flux_hbef_daily(chem_df = conc_data_prep, 
                                                  q_df = q_data_prep, 
                                                  ws_size = area))
        
        if(inherits(hbef_flux, 'try-error')){
        
        } else{
            write_feather(hbef_flux, glue('{directory}/{site_code}.feather'))
        }
        
        # EGRET (WRTDS) 
        directory <- glue('data/fluxes/{thinning_interval}/{var}/egret/')
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
        }
        
        source('source/flux_method_egret_daily.R')
        egret_flux <- try(adapt_ms_egret(chem_df = conc_data_prep, 
                                         q_df = q_data_prep, 
                                         ws_size = area,
                                         lat = lat, 
                                         long = long))
        
        if(inherits(egret_flux, 'try-error')){
            
        } else{
            write_rds(egret_flux, glue('{directory}/{site_code}.rds'))
        }
        
        # Next Method 
    }
    
}
 


