#Identify usgs gauges with vars of intest

library(tidyverse)
library(lubridate)
library(USAboundaries)
library(feather)
library(glue)

states <- USAboundaries::state_codes$state_abbr

all_sensors <- tibble()
# all_sensors <- tibble()
for(i in 1:length(states)){
    sites_99136 <- dataRetrieval::readNWISdata(parameterCd = '99136',
                                               service = 'iv',
                                               stateCd = states[i])

    if(nrow(sites_99136) == 0){
        sites_99136 <- tibble()
    } else {
        sites_99136 <- sites_99136 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = 99136,
                   state = !!states[i])
    }
    
    sites_51290 <- dataRetrieval::readNWISdata(parameterCd = '51290',
                                               service = 'iv',
                                               stateCd = states[i])
    if(nrow(sites_51290) == 0){
        sites_51290 <- tibble()
    } else {
        sites_51290 <- sites_51290 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = 51290,
                   state = !!states[i])
    }
    
    sites_83561 <- dataRetrieval::readNWISdata(parameterCd = '83561',
                                               service = 'iv',
                                               stateCd = states[i])
    if(nrow(sites_83561) == 0){
        sites_83561 <- tibble()
    } else {
        sites_83561 <- sites_83561 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = 83561,
                   state = !!states[i])
    }
    
    sites_99133 <- dataRetrieval::readNWISdata(parameterCd = '99133',
                                               service = 'iv',
                                               stateCd = states[i])
    if(nrow(sites_99133) == 0){
        sites_99133 <- tibble()
    } else {
        sites_99133 <- sites_99133 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = 99133,
                   state = !!states[i])
    }
    
    sites_99137 <- dataRetrieval::readNWISdata(parameterCd = '99137',
                                               service = 'iv',
                                               stateCd = states[i])
    if(nrow(sites_99137) == 0){
        sites_99137 <- tibble()
    } else {
        sites_99137 <- sites_99137 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = 99137,
                   state = !!states[i])
    }

    sites_00095 <- dataRetrieval::readNWISdata(parameterCd = '00095',
                                               service = 'iv',
                                               stateCd = states[i])

    if(nrow(sites_00095) == 0){
        sites_00095 <- tibble()
    } else {
        sites_00095 <- sites_00095 %>%
            select(site_no, dateTime) %>%
            mutate(parm_cd = '00095',
                   state = !!states[i])
    }
    
    all_sensors <- rbind(all_sensors, sites_99136, sites_51290, sites_83561, sites_99133, sites_99137, sites_00095)
    
    
}

site_info <- dataRetrieval::readNWISsite(all_sensors$site_no) %>%
    filter(site_tp_cd == 'ST') %>%
    filter(!is.na(drain_area_va))

all_sensors <- all_sensors %>%
    filter(site_no %in% site_info$site_no)

write_csv(all_sensors, 'data/general/usgs_sensors.csv')
# all_sensors <- read.csv('data/general/usgs_sensors.csv', colClasses = 'character')
codes <- unique(all_sensors$site_no)

sites_with_q <- dataRetrieval::whatNWISdata(siteNumber = codes, parameterCd="00060")


site_var_data <- all_sensors %>%
    filter(site_no %in% !!sites_with_q$site_no) %>%
    select(site_code = site_no, start_date = dateTime, parm_cd)

write_csv(site_var_data, 'data/general/site_var_data.csv')


# site_areas <- site_info %>%
#     select(site_no, station_nm, drain_area_va)
# 
# site_info_codes <- left_join(all_n_sensors, site_areas)
# 
# dir.create('flux_comparison/data/usgs', showWarnings = FALSE)
# for(i in 1:nrow(site_info_codes)){
#     nitrogen <- dataRetrieval::readNWISuv(site_info_codes[[i,1]], site_info_codes[[i,3]])
#     
#     
#     write_feather(nitrogen, glue('flux_comparison/data/usgs/{s}_{p}.feather',
#                                  s = site_info_codes[[i,1]],
#                                  p = site_info_codes[[i,3]]))
# }
# 
# all_fils <-list.files('flux_comparison/data/usgs/', full.names = T)
# nitrogen <- read_feather(grep('12212050', all_fils, value = T))
# nitrogen <- dataRetrieval::readNWISuv('06899900', '99137')
# 
# nitrogen %>%
#     filter()
# ggplot(aes(dateTime, X_99133_00000)) +
#     geom_line()
# 
