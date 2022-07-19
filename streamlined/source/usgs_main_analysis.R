# load sources ####
#remotes::install_github("leppott/ContDataQC")
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
library(ContDataQC)


source('streamlined/source/usgs_helpers.R')
#source('streamlined/source/usgs_retrieval.R')
source('source/helper_functions.R')

# read in USGS sites w continuous Nitrate ####
usgs_n <- read_csv("streamlined/data/site/usgs_nitrate_sites.csv") %>%
    filter(site_code != '03275500',
           site_code != '03381495',
           site_code != '01646500')

# find good sites #####
# Q and N present
# 85% data coverage by day
good_list <- tibble(site_code = as.character(), 
                    index = as.integer())

for(i in 1:nrow(usgs_n)){ #check for good sites
    #for(i in 1:length(good_sites)){
    # select site #####
    ### prep data ####
    #### read in raw ######
    site_no <- usgs_n$site_code[i]
    area <- usgs_n$ws_area_ha[i]
    lat <- usgs_n$lat[i]
    long <- usgs_n$long[i]

  tryCatch(
    expr = {
      print(paste("---- attempting read:", site_no))
      raw_data_n <- read_feather(glue('streamlined/data/raw/nitrate_nitrite_mgl/', site_no, '.feather'))
      raw_data_q <- read_feather(glue('streamlined/data/raw/q_cfs/', site_no, '.feather')) %>%
          group_by(datetime = date(datetime)) %>%
          summarize(val = mean(val)) %>%
          mutate(var = 'q_cfs',
                 site_code = site_no) %>%
          select(site_code, datetime, var, val)
    },
    error = function(e) {
      print('ERROR', site_no)
    }
    )
    
    #### isolate to full years######
    
    q_check <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(datetime, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)
## +
##   labs(caption = "(Pauloo, et al. 2017)")
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
    if(length(good_years) != 0){
        
        append <- tibble(site_code = site_no, index = i)
        good_list <- rbind(good_list, append)
        print(paste(' SUCCESS', site_no, '  years:', length(good_years)))
    }else{
      print('fail')
    }
}

i <- 32

for(i in 1:nrow(good_list)){
# select site #####
### prep data ####
#### read in raw ######
index <- good_list$index[i]
site_no <- usgs_n$site_code[index]
area <- usgs_n$ws_area_ha[index]
lat <- usgs_n$lat[index]
long <- usgs_n$long[index]

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

#### join data and cut to good years ######
# at this site, q is daily and nitrate is high frequency, so averaging n to daily
daily_data_n <- raw_data_n %>%
    mutate(date = date(datetime)) %>%
    group_by(date) %>%
    summarize(val = mean(val)) %>%
    mutate(site_code = site_no, var = 'nitrate_nitrite_mgl') %>%
    select(site_code, datetime = date, var, val)

raw_data_full_pre <- rbind(daily_data_n, raw_data_q) %>%
    pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
    mutate(wy = water_year(datetime, origin = 'usgs'),
           q_lps = q_cfs*28.316847) %>%
    filter(wy %in% good_years)



# Loop through good years ####
## t <- 1
# needs DAILY q and any chem
for(t in 1:length(good_years)){
target_year <- as.numeric(as.character(good_years[t]))
raw_data_full <- raw_data_full_pre %>%
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

rating <- summary(lm(model_data$c_log ~ model_data$q_log, singular.ok = T))

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
rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0

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
    filter(wy == target_year) %>%
    select(site_code, date = datetime, q_lps, wy)


### THIN data to selected intervals #######
thinned_daily_c <- raw_data_n %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    mutate(date = lubridate::date(datetime),
           wy = water_year(datetime, origin = 'usgs')) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val) %>%
    filter(wy == target_year)

thinned_weekly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    filter(lubridate::wday(datetime) == 3) %>%
    mutate(date = lubridate::date(datetime)) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)%>%
    filter(wy == target_year)

thinned_biweekly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
    mutate(date = lubridate::date(datetime)) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)%>%
    filter(wy == target_year)

thinned_monthly_c <- raw_data_n %>%
    mutate(wy = water_year(datetime, origin = 'usgs')) %>%
    filter(wy %in% good_years) %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    mutate(date = date(datetime)) %>%
    filter(day(date) == 1) %>%
    distinct(date, .keep_all = T) %>%
    select(-datetime, -var) %>%
    rename(con = val)%>%
    filter(wy == target_year)

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
library(EGRET)
source('streamlined/source/egret_overwrites.R')

calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long) {

  tryCatch(
    expr = {
      egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long)

      flux_from_egret <- egret_results$Daily$FluxDay %>%
        sum(.)/(1000 * area)
    },
    error = function(e) {
      print('ERROR, EGRET FAILED')
    })
  return(flux_from_egret)
}

flux_from_daily_wrtds <- calculate_wrtds(
  chem_df = thinned_daily_c,
                                 q_df = prep_data_q,
                                 ws_size = area,
                                 lat = lat,
                                 long = long)

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

rating <- summary(lm(model_data$c_log ~ model_data$q_log, singular.ok = T))

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
rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0
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
                    flux = c(flux_from_daily_pw,
                             flux_from_daily_beale,
                             flux_from_daily_rating,
                             flux_from_daily_wrtds$Daily,
                             flux_from_daily_comp$flux[1]),
                    site_code = flux_from_daily_comp$site_code[1], 
                    method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
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

## save flux out of loop ####
directory <- glue('streamlined/out/{wy}',
                  wy = target_year)
if(!dir.exists(directory)){
    dir.create(directory, recursive = TRUE)
            }
file_path <- glue('{directory}/{s}.feather',
                  s = site_no)

write_feather(out_frame, file_path)
## take meta and save out of loop ####
meta_n <- raw_data_n %>%
    filter(wy == target_year)

out_meta <- tibble(max_q_lps = max(prep_data_q$q_lps, na.rm = T),
                   min_q_lps = min(prep_data_q$q_lps, na.rm = T),
                   mean_q_lps = mean(prep_data_q$q_lps, na.rm = T),
                   sd_q_lps = sd(prep_data_q$q_lps, na.rm = T),
                   max_n = max(meta_n$val, na.rm = T),
                   min_n = min(meta_n$val, na.rm = T),
                   mean_n = mean(meta_n$val, na.rm = T),
                   sd_n = sd(meta_n$val, na.rm = T))


directory <- glue('streamlined/data/meta/{wy}',
                  wy = target_year)
if(!dir.exists(directory)){
    dir.create(directory, recursive = TRUE)
}
file_path <- glue('{directory}/{s}.feather',
                  s = site_no)

write_feather(out_meta, file_path)
    } # end year loop
} # end site loop
