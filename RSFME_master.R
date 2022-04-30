library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(feather)
library(zoo)
library(here)

###### prep data ####
# read in eagle creek
eagle <- readNWISuv('03353200', parameterCd = c('00060', '99137'),
                    startDate = '2019-10-01', endDate = '2020-09-30',
                    tz = 'EST')

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

# create weekly data for testing
chem_df <- eagle_chem %>%
    mutate(week = floor_date(date, unit = 'week')) %>%
    group_by(week) %>%
    arrange(date) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(date, con = nitrate_mgL)%>%
    mutate(date = floor_date(date, unit = 'days'))

q_df <- eagle_q %>%
    mutate(q_lps = q_cfs*28.316847) %>%
    select(date, q_lps)

###### test completed functions
source('source/flux_method_hbef_daily.R')
test_hbef <- estimate_flux_hbef_daily(chem_df = chem_df, q_df = q_df, ws_size = eagle_ws_area_ha)
test_hbef

ggplot(test_hbef, aes(x = date, y = flux_daily_kg_ha))+
  geom_point()+
  geom_line()

source('source/flux_method_hbef_annual.R')
estimate_flux_hbef_annual(chem_df = chem_df, q_df = q_df, ws_size = eagle_ws_area_ha)

source('source/flux_method_fernow_weekly.R')
estimate_flux_fernow_weekly(chem_df = chem_df, q_df = q_df, ws_size = eagle_ws_area_ha)

source('source/flux_method_fernow_annual.R')
estimate_flux_fernow_annual(chem_df = chem_df, q_df = q_df, ws_size = eagle_ws_area_ha)


