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
library(EGRET)

## Local Imports
# flux helpers
source('source/egret_overwrites.R')
source('source/flux_methods.R')

# general helpers
source('source/usgs_helpers.R')
source('source/helper_functions.R')

# read in USGS sites w continuous Nitrate ####
# TODO: replace with dynamic call to USGS site info, or, pulling from dynamic source at least
usgs_n <- read_csv("data/site/usgs_nitrate_sites.csv") %>%
    filter(site_code != '03275500',
           site_code != '03381495',
           site_code != '01646500')

# find good sites #####
# Q and N present
# 85% data coverage by day
good_list <- tibble(site_code = as.character(),
                    index = as.integer())

for(i in 1:nrow(usgs_n)){ #check for good sites
    # site information
    site_no <- usgs_n$site_code[i]
    area <- usgs_n$ws_area_ha[i]
    lat <- usgs_n$lat[i]
    long <- usgs_n$long[i]

    tryCatch(
      expr = {
          # read in Nitrate and Q data
          raw_data_n <- read_feather(glue('data/raw/nitrate_nitrite_mgl/', site_no, '.feather'))
          raw_data_q <- read_feather(glue('data/raw/q_cfs/', site_no, '.feather')) %>%
              group_by(datetime = date(datetime)) %>%
              summarize(val = mean(val)) %>%
              mutate(var = 'q_cfs',
                     site_code = site_no) %>%
              select(site_code, datetime, var, val)
        },
      error = function(e) {
        # TODO: fix printing error on digits 1477050
        print('ERROR: ', site_no, ' could not read in raw data')
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
        print(paste('---- SUCCESS', site_no, '  years:', length(good_years)))
    } else{
        print(paste('---- FAIL', site_no, 'site has no years with sufficient data'))
    }
}

## write.csv(good_list, 'data/good_sites.csv')
## write.csv(good_years, 'data/good_years.csv')

# failed WRTDS flux calcs
## failed.df <- data.frame(matrix(ncol = 5, nrow= 1))
## colnames(failed.df) <- c('site_no', 'year', 'method', 'thinning')

# df to populate with annual flux values by method

out_frame_main <- tibble(wy = as.integer(),
                    flux = as.numeric(),
                    site_code = as.character(),
                    method = as.character(),
                    thin = as.character())

# loop thru 'good' sites and perform thinnings and flux calcs
for(i in 1:nrow(good_list)){
    # site information
    index <- good_list$index[i]
    site_no <- usgs_n$site_code[index]
    area <- usgs_n$ws_area_ha[index]
    lat <- usgs_n$lat[index]
    long <- usgs_n$long[index]


    # read in raw data
    raw_data_n <- read_feather(glue('data/raw/nitrate_nitrite_mgl/', site_no, '.feather'))
    raw_data_q <- read_feather(glue('data/raw/q_cfs/', site_no, '.feather')) %>%
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
    # at some sites, q is daily and nitrate is high frequency, so averaging n to daily
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

    # print info
    writeLines(paste('\n\n',
                   'START OF FLUX CALCS, SITE:', site_no, '\n',
                   'YEARS:', good_years, '\n',
                   'site ', i, 'of', nrow(good_list), '\n'))

    print('-------------------------------------')

    # Loop through good years ####
    # needs DAILY q and any chem
    for(t in 1:length(good_years)){
        target_year <- as.numeric(as.character(good_years[t]))

        # print info
        writeLines(paste('--------- year:', target_year,
                         '| 1 of ', length(good_years)))

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
            # NOTE: some sites throw NAs in record, ask NG if cool to replace?
            summarize(flux = warn_sum(flux)) %>%
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

        ### DAILY ######
        print('DAILY')
        chem_df <- thinned_daily_c

        ###### calculate period weighted#########
        flux_from_daily_pw <- calculate_pw(chem_df, prep_data_q)

        ###### calculate beale ######
        flux_from_daily_beale <- calculate_beale(chem_df, prep_data_q)

        ##### calculate rating #####
        flux_from_daily_rating <- calculate_rating(chem_df, prep_data_q)

        ###### calculate wrtds ######
        flux_from_daily_wrtds <- calculate_wrtds(
          chem_df = chem_df,
          q_df = prep_data_q,
          ws_size = area,
          lat = lat,
          long = long)

        # plot TS of true flux and wrtds
        ## true_flux_d <- rating_filled_df %>%
        ##   mutate(flux = con_com*q_lps*86400*(1/area)*1e-6)
        ## egret_flux <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long)
        ## plot(true_flux_d$datetime, true_flux_d$flux)
        ## plot(egret_flux$Daily$Date, egret_flux$Daily$FluxDay)

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
        flux_from_daily_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

        ##### congeal daily ####
        daily_out <- tibble(wy = flux_from_daily_comp$wy[1],
                            flux = c(flux_from_daily_pw,
                                     flux_from_daily_beale,
                                     flux_from_daily_rating,
                                     flux_from_daily_wrtds,
                                     flux_from_daily_comp$flux[1]),
                            site_code = site_no,
                            method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
                            thin = 'daily')


        ## WEEKLY ######
        print('WEEKLY')
        chem_df <- thinned_weekly_c

        ###### calculate period weighted#########
        flux_from_weekly_pw <- calculate_pw(chem_df, prep_data_q)

        ###### calculate beale ######
        flux_from_weekly_beale <- calculate_beale(chem_df, prep_data_q)

        ##### calculate rating #####
        flux_from_weekly_rating <- calculate_beale(chem_df, prep_data_q)

        #### calculate wrtds ######
        flux_from_weekly_wrtds <- calculate_wrtds(chem_df,
                                                  prep_data_q,
                                                  ws_size = area,
                                                  lat = lat,
                                                  long = long)

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = prep_data_q)
        # calculate annual flux from composite
        flux_from_weekly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

        ##### congeal weekly ####
        weekly_out <- tibble(wy = flux_from_weekly_comp$wy[1],
                             flux = c(flux_from_weekly_pw,
                                      flux_from_weekly_beale,
                                      flux_from_weekly_rating,
                                      flux_from_weekly_wrtds,
                                      flux_from_weekly_comp$flux[1]),
                            site_code = site_no,
                            method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
                            thin = 'weekly')

        ## BIWEEKLY ######
        print('BIWEEKLY')
        chem_df <- thinned_biweekly_c

        ###### calculate period weighted#########
        flux_from_biweekly_pw <- calculate_pw(chem_df, prep_data_q)

        ###### calculate beale ######
        flux_from_biweekly_beale <- calculate_beale(chem_df, prep_data_q)

        ##### calculate rating #####
        flux_from_biweekly_rating <- calculate_rating(chem_df, prep_data_q)

        #### calculate wrtds ######
        flux_from_biweekly_wrtds <- calculate_wrtds(chem_df,
                                                  prep_data_q,
                                                  ws_size = area,
                                                  lat = lat,
                                                  long = long)

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = prep_data_q)
        # calculate annual flux from composite
        flux_from_biweekly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

        ##### congeal biweekly ####
        biweekly_out <- tibble(wy = flux_from_biweekly_comp$wy[1],
                               flux = c(flux_from_biweekly_pw,
                                        flux_from_biweekly_beale,
                                        flux_from_biweekly_rating,
                                        flux_from_biweekly_wrtds,
                                        flux_from_biweekly_comp$flux[1]),
                             site_code = site_no,
                             method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
                             thin = 'biweekly')

        ## MONTHLY ######
        print('MONTHLY')
        chem_df <- thinned_monthly_c

        ###### calculate period weighted#########
        flux_from_monthly_pw <- calculate_pw(chem_df, prep_data_q)

        ###### calculate beale ######
        flux_from_monthly_beale <- calculate_beale(chem_df, prep_data_q)

        ##### calculate rating #####
        flux_from_monthly_rating <- calculate_rating(chem_df, prep_data_q)

        #### calculate wrtds ######
        flux_from_monthly_wrtds <- calculate_wrtds(chem_df,
                                                  prep_data_q,
                                                  ws_size = area,
                                                  lat = lat,
                                                      long = long)

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = prep_data_q)
        # calculate annual flux from composite
        flux_from_monthly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

        ##### congeal monthly ####
        monthly_out <- tibble(wy = flux_from_monthly_comp$wy[1],
                              flux = c(flux_from_monthly_pw,
                                       flux_from_monthly_beale,
                                       flux_from_monthly_rating,
                                       flux_from_monthly_wrtds,
                                       flux_from_monthly_comp$flux[1]),
                                site_code = site_no,
                                method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
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
        # TODO: make WRTDS run on quarterly data
        ## flux_from_quarterly_wrtds <- calculate_wrtds(chem_df,
        ##                                           prep_data_q,
        ##                                           ws_size = area,
        ##                                           lat = lat,
        ##                                           long = long)

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = prep_data_q)
        # calculate annual flux from composite
        flux_from_quarterly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df, site_no = site_no)

        ##### congeal quarterly ####
        quarterly_out <- tibble(wy = flux_from_quarterly_comp$wy[1],
                                flux = c(flux_from_quarterly_pw,
                                         flux_from_quarterly_beale,
                                         flux_from_quarterly_rating,
                                         ## flux_from_quarterly_wrtds,
                                         NA,
                                         flux_from_quarterly_comp$flux[1]),
                               site_code = site_no,
                               method = c('pw', 'beale', 'rating', 'wrtds', 'composite'),
                               thin = 'quarterly')

        ## congeal results ####
        out_frame <- rbind(true_flux, daily_out, weekly_out, biweekly_out, monthly_out, quarterly_out)

        # all results #
        out_frame_main <- rbind(out_frame_main, out_frame)

        print('----- head main outframe ------')
        print(head(out_frame_main))
        print('NROW:')
        print(nrow(out_frame_main))
        print('-------------------------------')

        ## save flux out of loop ####
        directory <- 'out/'
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
                    }
        file_path <- 'out/usgs_annual_flux.feather'

        write_feather(out_frame_main, file_path)

        ## save flux out of loop ####
        directory <- glue('out/{wy}',
                          wy = target_year)
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
                    }
        file_path <- glue('{directory}/{s}.feather',
                          s = site_no)

        write_feather(out_frame, file_path)

        ## take meta and save out of loop ####
        meta_n <- raw_data_n %>%
          mutate(wy = water_year(datetime)) %>%
            filter(wy == target_year)

        out_meta <- tibble(max_q_lps = max(prep_data_q$q_lps, na.rm = T),
                           min_q_lps = min(prep_data_q$q_lps, na.rm = T),
                           mean_q_lps = mean(prep_data_q$q_lps, na.rm = T),
                           sd_q_lps = sd(prep_data_q$q_lps, na.rm = T),
                           max_n = max(meta_n$val, na.rm = T),
                           min_n = min(meta_n$val, na.rm = T),
                           mean_n = mean(meta_n$val, na.rm = T),
                           sd_n = sd(meta_n$val, na.rm = T))


        directory <- glue('data/meta/{wy}',
                          wy = target_year)
        if(!dir.exists(directory)){
            dir.create(directory, recursive = TRUE)
        }

        file_path <- glue('{directory}/{s}.feather',
                          s = site_no)

        write_feather(out_meta, file_path)


    } # end year loop
} # end site loop
## fail on :
##  0167300055
