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


source('streamlined/source/usgs_helpers.R')
#source('streamlined/source/usgs_retrieval.R')
source('source/helper_functions.R')

# read in USGS sites w continuous Nitrate ####
usgs_n <- read_csv("streamlined/data/site/usgs_nitrate_sites.csv")

good_sites <- c(1,3,8,9,12,13,14,17,18)
for(i in good_sites){
# select site #####
### prep data ####
#### read in raw ######
site_no <- usgs_n$site_code[i]
area <- usgs_n$ws_area_ha[i]
lat <- usgs_n$lat[i]
long <- usgs_n$long[i]

raw_data_n <- read_feather(glue('streamlined/data/raw/nitrate_nitrite_mgl/', site_no, '.feather'))
raw_data_q <- read_feather(glue('streamlined/data/raw/q_cfs/', site_no, '.feather')) %>%
    group_by(datetime = date(datetime)) %>%
    summarize(val = mean(val)) %>%
    mutate(var = 'q_cfs',
           site_code = site_no) %>%
    select(site_code, datetime, var, val)

#### isolate to full years######

q_check <- raw_data_q %>%
    mutate(date = date(datetime)) %>%
    distinct(., date, .keep_all = TRUE) %>%
    mutate(water_year = water_year(datetime, origin = "usgs")) %>%
    group_by(water_year) %>%
    summarise(n = n()) %>%
    filter(n >= 311)

conc_check <- raw_data_n %>%
    mutate(date = date(datetime)) %>%
    distinct(., date, .keep_all = TRUE) %>%
    mutate(water_year = water_year(date, origin = "usgs")) %>%
    group_by(water_year) %>%
    summarise(n = n()) %>%
    filter(n >= 311)

q_good_years <- q_check$water_year
conc_good_years <- conc_check$water_year

good_years <- q_good_years[q_good_years %in% conc_good_years]
# if(length(good_years) != 0){
#     print(i)
# }else{print('fail')}


#### join data and cut to good years ######
# at this site, q is daily and nitrate is high frequency, so averaging n to daily
daily_data_n <- raw_data_n %>%
    mutate(date = date(datetime)) %>%
    group_by(date) %>%
    summarize(val = mean(val)) %>%
    mutate(site_code = site_no, var = 'nitrate_nitrite_mgl') %>%
    select(site_code, datetime = date, var, val)

raw_data_full <- rbind(daily_data_n, raw_data_q) %>%
    pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
    mutate(wy = water_year(datetime, origin = 'usgs'),
           q_lps = q_cfs*28.316847) %>%
    filter(wy %in% good_years)



# Loop through good years ####
# needs DAILY q and any chem
for(t in length(good_years)){
target_year <- as.numeric(as.character(good_years[t]))
raw_data_full <- raw_data_full %>%
    mutate(wy = as.numeric(as.character(wy))) %>%
    filter(wy == target_year)
### TRUTH (via composite) ######
# first make c:q rating
paired_df <- raw_data_full%>%
    rename(con = nitrate_nitrite_mgl) %>%
    na.omit() %>%
    filter(q_lps > 0,
           is.finite(q_lps))

q_log <- log10(paired_df$q_lps)
c_log <- log10(paired_df$con)
model_data <- tibble(c_log, q_log) %>%
    filter(is.finite(c_log),
           is.finite(q_log))%>%
    na.omit()

rating <- summary(lm(model_data$c_log ~ model_data$q_log))

# extract model info
intercept <- rating$coefficients[1]
slope <- rating$coefficients[2]

# create modeled c, calc residuals, adjust modeled c by interpolated residuals
rating_filled_df <- raw_data_full %>%
    mutate(con_reg = 10^(intercept+(slope*log10(q_lps)))) %>%
    select(datetime, con_reg) %>%
    full_join(., raw_data_full, by = 'datetime') %>%
    select(site_code, datetime, con = nitrate_nitrite_mgl, con_reg, q_lps, wy) %>%
    mutate(res = con_reg-con,
           res = imputeTS::na_interpolation(res),
           con_com = con_reg-res)

# calculate true annual flux
true_flux <- rating_filled_df %>%
    mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
    group_by(wy) %>%
    summarize(flux = sum(flux)) %>%
    mutate(site_code = site_no,
           method = 'true',
           thin = 'none')

###### prep q data#####
prep_data_q <- raw_data_q %>%
    mutate(wy = water_year(datetime, origin = 'usgs'),
           q_lps = val*28.316847) %>%
    filter(wy %in% good_years) %>%
    select(site_code, date = datetime, q_lps, wy)


### THIN data to selected intervals #######
thinned_daily_c <- raw_data_n %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    mutate(date = lubridate::date(datetime),
           wy = water_year(datetime, origin = 'usgs')) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val) %>%
    filter(wy %in% good_years)

thinned_weekly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    filter(lubridate::wday(datetime) == 3) %>%
    mutate(date = lubridate::date(datetime)) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)

thinned_biweekly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
    mutate(date = lubridate::date(datetime)) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)

thinned_monthly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    mutate(date = date(datetime)) %>%
    filter(day(date) == 1) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)

nmonth <- nrow(thinned_monthly_c)

thinned_quarterly_c <- rbind(thinned_monthly_c[1,], 
                             thinned_monthly_c[as.integer(.5*nmonth),],
                             thinned_monthly_c[as.integer(.75*nmonth),],
                             thinned_monthly_c[nmonth,])

source('source/flux_method_egret_daily.R')
### Riverload conversion function #####
prep_raw_for_riverload <- function(chem_df, q_df){
    conv_q <- q_df %>%
        mutate(datetime = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
        mutate(flow = q_lps*0.001) %>% # convert lps to cubic meters per second)
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
### DAILY ######
chem_df <- thinned_daily_c
###### calculate period weighted#########
calculate_pw <- function(chem_df, q_df){
rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)

flux_from_pw <- method1(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
return(flux_from_pw)
}
flux_from_daily_pw <- calculate_pw(chem_df, prep_data_q)

###### calculate beale ######
calculate_beale <- function(chem_df, q_df){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
flux_from_beale <- beale.ratio(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
return(flux_from_beale)
}
flux_from_daily_beale <- calculate_beale(chem_df, prep_data_q)
##### calculate rating #####
calculate_rating <- function(chem_df, q_df){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
    flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1) %>%
        sum(.)/(1000*area)
    return(flux_from_reg)
}
flux_from_daily_rating <- calculate_rating(chem_df, prep_data_q)

###### calculate wrtds ######
#Currently not working
# egret_q <- raw_data_q %>%
#     mutate(q_lps = val*28.316847) %>%
#     select(date = datetime, q_lps)


# egret_flux <- try(adapt_ms_egret(chem_df = thinned_daily_c,
#                                  q_df = prep_data_q,
#                                  ws_size = area,
#                                  lat = lat,
#                                  long = long))


###### calculate composite ######
generate_residual_corrected_con <- function(chem_df, q_df){
# first make c:q rating
paired_df <- q_df %>%
    full_join(chem_df, by = c('date', 'site_code', 'wy')) %>%
    na.omit() %>%
    filter(q_lps > 0,
           is.finite(q_lps))

q_log <- log10(paired_df$q_lps)
c_log <- log10(paired_df$con)
model_data <- tibble(c_log, q_log) %>%
    filter(is.finite(c_log),
           is.finite(q_log))%>%
    na.omit()

rating <- summary(lm(model_data$c_log ~ model_data$q_log))

# extract model info
intercept <- rating$coefficients[1]
slope <- rating$coefficients[2]

# create modeled c, calc residuals, adjust modeled c by interpolated residuals
rating_filled_df <- q_df %>%
    mutate(con_reg = 10^(intercept+(slope*log10(q_lps)))) %>%
    select(date, con_reg, q_lps) %>%
    full_join(., chem_df, by = 'date') %>%
    select(site_code, date, con, con_reg, q_lps, wy) %>%
    mutate(res = con_reg-con,
           res = imputeTS::na_interpolation(res),
           con_com = con_reg-res,
           site_code = site_no,
           wy = water_year(date, origin = 'usgs'))

return(rating_filled_df)
}

rating_filled_df <- generate_residual_corrected_con(chem_df = thinned_daily_c,
                                q_df = prep_data_q)

# calculate annual flux from composite
calculate_composite_from_rating_filled_df <- function(rating_filled_df){
flux_from__comp <- rating_filled_df %>%
    mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
    group_by(wy) %>%
    summarize(flux = sum(flux)) %>%
    mutate(site_code = site_no)
}
flux_from_daily_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

##### congeal daily ####
daily_out <- tibble(wy = flux_from_daily_comp$wy[1], 
                    flux = c(flux_from_daily_pw, flux_from_daily_beale, 
                                                flux_from_daily_rating, flux_from_daily_comp$flux[1]),
                    site_code = flux_from_daily_comp$site_code[1], 
                    method = c('pw', 'beale', 'rating', 'composite'),
                    thin = 'daily')


## WEEKLY ######
chem_df <- thinned_weekly_c
###### calculate period weighted#########
flux_from_weekly_pw <- calculate_pw(chem_df, prep_data_q)

###### calculate beale ######
flux_from_weekly_beale <- calculate_beale(chem_df, prep_data_q)

##### calculate rating #####
flux_from_weekly_rating <- calculate_beale(chem_df, prep_data_q)

#### calculate wrtds ######

#### calculate composite ######
rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                    q_df = prep_data_q)
# calculate annual flux from composite
flux_from_weekly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

##### congeal weekly ####
weekly_out <- tibble(wy = flux_from_weekly_comp$wy[1], 
                    flux = c(flux_from_weekly_pw, flux_from_weekly_beale, 
                             flux_from_weekly_rating, flux_from_weekly_comp$flux[1]),
                    site_code = flux_from_weekly_comp$site_code[1], 
                    method = c('pw', 'beale', 'rating', 'composite'),
                    thin = 'weekly')

## BIWEEKLY ######
chem_df <- thinned_biweekly_c
###### calculate period weighted#########
flux_from_biweekly_pw <- calculate_pw(chem_df, prep_data_q)

###### calculate beale ######
flux_from_biweekly_beale <- calculate_beale(chem_df, prep_data_q)

##### calculate rating #####
flux_from_biweekly_rating <- calculate_rating(chem_df, prep_data_q)

#### calculate wrtds ######

#### calculate composite ######
rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                    q_df = prep_data_q)
# calculate annual flux from composite
flux_from_biweekly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

##### congeal biweekly ####
biweekly_out <- tibble(wy = flux_from_biweekly_comp$wy[1], 
                     flux = c(flux_from_biweekly_pw, flux_from_biweekly_beale, 
                              flux_from_biweekly_rating, flux_from_biweekly_comp$flux[1]),
                     site_code = flux_from_biweekly_comp$site_code[1], 
                     method = c('pw', 'beale', 'rating', 'composite'),
                     thin = 'biweekly')

## MONTHLY ######
chem_df <- thinned_weekly_c
###### calculate period weighted#########
flux_from_monthly_pw <- calculate_pw(chem_df, prep_data_q)

###### calculate beale ######
flux_from_monthly_beale <- calculate_beale(chem_df, prep_data_q)

##### calculate rating #####
flux_from_monthly_rating <- calculate_rating(chem_df, prep_data_q)

#### calculate wrtds ######

#### calculate composite ######
rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                    q_df = prep_data_q)
# calculate annual flux from composite
flux_from_monthly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

##### congeal monthly ####
monthly_out <- tibble(wy = flux_from_monthly_comp$wy[1], 
                        flux = c(flux_from_monthly_pw, flux_from_monthly_beale, 
                                 flux_from_monthly_rating, flux_from_monthly_comp$flux[1]),
                        site_code = flux_from_monthly_comp$site_code[1], 
                        method = c('pw', 'beale', 'rating', 'composite'),
                        thin = 'monthly')

## QUARTERLY ######
chem_df <- thinned_quarterly_c
###### calculate period weighted#########
flux_from_quarterly_pw <- calculate_pw(chem_df, prep_data_q)

###### calculate beale ######
flux_from_quarterly_beale <- calculate_beale(chem_df, prep_data_q)

##### calculate rating #####
flux_from_quarterly_rating <- calculate_beale(chem_df, prep_data_q)

#### calculate wrtds ######

#### calculate composite ######
rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                    q_df = prep_data_q)
# calculate annual flux from composite
flux_from_quarterly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

##### congeal quarterly ####
quarterly_out <- tibble(wy = flux_from_quarterly_comp$wy[1], 
                       flux = c(flux_from_quarterly_pw, flux_from_quarterly_beale, 
                                flux_from_quarterly_rating, flux_from_quarterly_comp$flux[1]),
                       site_code = flux_from_quarterly_comp$site_code[1], 
                       method = c('pw', 'beale', 'rating', 'composite'),
                       thin = 'quarterly')

## congeal results ####
out_frame <- rbind(true_flux, daily_out, weekly_out, biweekly_out, monthly_out, quarterly_out)

## save out of loop ####
directory <- glue('streamlined/out/{wy}',
                  wy = target_year)
if(!dir.exists(directory)){
    dir.create(directory, recursive = TRUE)
            }
file_path <- glue('{directory}/{s}.feather',
                  s = site_no)

write_feather(out_frame, file_path)
    }
}
