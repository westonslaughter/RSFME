# Streamlined USGS data retrieval for sites w continuous Nitrate
# load sources ####
library(dataRetrieval)
library(here)
library(dplyr)
library(feather)
library(zoo)
library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(glue)
library(feather)
library(zoo)
library(lfstat)
library(RiverLoad)
# read in df of USGS sites w continous Nitrate
usgs <- read.csv("streamlined/data/site/usgs_nitrate_sites.csv",
                 colClasses = "character")

# Parameter
# Nitrate
parm_cd <- "99133"
var <- "nitrate_nitrite_mgl"

failed_sites <- c()
for(i in 1:nrow(usgs)){

    # Load in site info
    site_info <- usgs[i,]
    site_code <- site_info$site_code

    # Download variable
    var_data <- readNWISuv(site_code,
                           parameterCd = parm_cd)

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

    directory <- glue('streamlined/data/raw/{v}',
                      v = var)

    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }

    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)

    write_feather(var_data, file_path)
    print(paste(site_code, 'done'))

}

# Discharge
sites <- unique(usgs$site_code)
sites <- sites[!sites %in% failed_sites]

for(i in 1:length(sites)){
    # Download Q data
    site_code <- sites[i]
    q_data <- readNWISdv(site_code,
                         parameterCd = '00060')

        date_name <- 'Date'


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

    q_directory <- glue('streamlined/data/raw/q_cfs')

    if(!dir.exists(q_directory)){
        dir.create(q_directory, recursive = TRUE)
    }

    file_path <- glue('{q_directory}/{s}.feather',
                      s = site_code)

    write_feather(q_data, file_path)
}
