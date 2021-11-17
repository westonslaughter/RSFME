library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(feather)
library(zoo)

setwd('C:/Users/gubbins/Dropbox/macrosheds/flux estimates/timing_analysis')

###### prep data ####
# read in eagle creek
eagle <- readNWISuv('03353200', parameterCd = c('00060', '99137'),
                    startDate = '2019-10-01', endDate = '2020-09-29')

eagle_ws_area_ha <- 27453.9

eagle_chem <- eagle %>%
    select(date = dateTime, nitrate_mgL = X_99137_00000) %>%
    mutate(nitrate_mgL = na.approx(nitrate_mgL, rule = 2))

eagle_q <- eagle %>%
    select(date = dateTime, q_cfs = X_00060_00000) %>%
    mutate(q_cfs = na.approx(q_cfs, rule = 2))

ggplot(eagle_chem, aes(x = date, y = nitrate_mgL))+
    geom_line()

ggplot(eagle_q, aes(x = date, y = q_cfs))+
    geom_line()+
    scale_y_log10()

# read in oconaluftee
ocon <- readNWISuv('03512000', parameterCd = c('00060', '99137'),
                   startDate = '2019-10-01', endDate = '2020-09-29')

ocon_ws_area_ha <- 47655.8

ocon_chem <- ocon %>%
    select(date = dateTime, nitrate_mgL = X_99137_00000) %>%
    mutate(nitrate_mgL = na.approx(nitrate_mgL, rule = 2))

ocon_q <- ocon %>%
    select(date = dateTime, q_cfs = X_00060_00000) %>%
    mutate(q_cfs = na.approx(q_cfs, rule = 2))

ggplot(ocon_chem, aes(x = date, y = nitrate_mgL))+
    geom_line()

ggplot(ocon_q, aes(x = date, y = q_cfs))+
    geom_line()+
    scale_y_log10()