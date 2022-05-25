# Build flux master 

# Load in packages 
library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(glue)
library(feather)
library(zoo)
library(here)

#### download raw data #### 
# Read in site data 
site_data <- read.csv('data/general/site_data.csv', colClasses = 'character')
variable_data <- read.csv('data/general/variable_data.csv', colClasses = 'character')

for(i in 1:nrow(site_data)){
    
    # Load in site info 
    site_info <- site_data[i,]
    site_code <- site_info$site_code
    parm_cd <- site_info$parm_cd
    
    var <- variable_data %>%
        filter(usgs_parm_cd == !!parm_cd) %>%
        pull(var)
    
    # Download Q data 
    q_data <- readNWISuv(site_code, parameterCd = '00060')
    
    # Incase a site has a varibale but not Q
    if(nrow(q_data) == 0) next
    
    q_data <- q_data %>%
        select(site_code = site_no, datetime = dateTime, val = X_00060_00000) %>%
        mutate(val = na.approx(val, rule = 2)) %>%
        mutate(var = 'q_cfs')
    
    q_directory <- glue('data/raw/q_cfs')
    
    if(!dir.exists(q_directory)){
        dir.create(q_directory)
    }
    
    file_path <- glue('{q_directory}/{s}.feather',
                      s = site_code)
    
    write_feather(q_data, file_path)
    
    # Download variable  
    no3_data <- readNWISuv(site_code, parameterCd = parm_cd)
    
    if(nrow(no3_data) == 0) next
    
    if(!paste0('X_', parm_cd, '_00000') %in% names(no3_data)) next
    
    no3_data <- no3_data %>%
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
    
    write_feather(no3_data, file_path)
    
}

# get sites that have both Q and chem 
q_sites <- list.files('data/raw/q_cfs/')
chem_sites <- list.files('data/raw/nitrate_nitrite_mgl/')
common_sites <- q_sites[q_sites %in% chem_sites]






#### Thin data to desired intervals ####
# You have more built out code for this from what I remember but I am just 
# doing a simple thinning here as an example 

thinning_intervals <- c('monthly', 'bi_weekely', 'weekely')

for(i in 1:length(common_sites)) {
    
    site <- common_sites[i]
    
    q_data <- read_feather(glue('data/raw/q_cfs/{site}'))
    chem_data <- read_feather(glue('data/raw/nitrate_nitrite_mgl/{site}'))
    
    var <- 'nitrate_nitrite_mgl'
    
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
            
            write_feather(chem_data_thin, glue('{directory}/{site}'))
            
        }
        
        if(thinning_intervals[p] == 'bi_weekely'){
            # add code here 
        }
        
        if(thinning_intervals[p] == 'weekely'){
            # add code here 
        }
    }
    
}


#### Calcaulte flues with various methods ####
flux_methods <- c('hbef', 'fernow')

for(i in 1:length(flux_methods)){
    
}


