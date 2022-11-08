library(tidyverse)
library(feather)
library(here)
library(glue)
library(lubridate)
library(EGRET)
library(macrosheds)
library(bootstrap)

source('source/helper_functions.R')
source('source/egret_overwrites.R')
source('source/flux_methods.R')
source('source/usgs_helpers.R')

data_dir <- here('streamlined/data/ms/hbef/')
site_files  <- list.files('streamlined/data/ms/hbef/discharge', recursive = F)
site_info  <- read_csv(here('streamlined/data/site/ms_site_info.csv'))

# non-fluxable solutes
non_fluxable <- c('anionCharge')

# df to populate with annual flux values by method
out_frame <- tibble(wy = as.integer(),
                    site_code = as.character(),
                    var = as.character(),
                    val = as.numeric(),
                    method = as.character(),
                    ms_reccomended = as.integer())
## i = 2
# Loop through sites #####
for(i in 1:length(site_files)){

    site_file <- site_files[i]
    site_code <- strsplit(site_file, split = '.feather')[[1]]

    area <- site_info %>%
        filter(site_code == !!site_code) %>%
        pull(ws_area_ha)

    lat <- site_info %>%
        filter(site_code == !!site_code) %>%
        pull(Y)

    long <- site_info %>%
        filter(site_code == !!site_code) %>%
        pull(X)

    # read in chemistry data
    raw_data_con_in <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0)

    # read in discharge data
    raw_data_q <- read_feather(here(glue(data_dir, '/discharge/', site_code, '.feather')))

    # initialize next loop
    solutes <- raw_data_con_in %>%
        ## filter(!str_detect(var, 'Charge'),
        ##        !str_detect(var, 'temp'),
        ##        !str_detect(var, 'pH')) %>%
        select(var) %>%
        unique() %>%
        pull(var)

  writeLines(paste("FLUX CALCS:", site_code))

  ## Loop through solutes at site #####
  ## j = 1
  for(j in 1:length(solutes)){
    writeLines(paste("site:", site_code,
                     "var:", solutes[j]))

    #set to target solute
    target_solute <- solutes[j]

    if(str_split_fixed(target_solute, pattern = '_', n = 2)[1,2] %in% non_fluxable){next}

    raw_data_con <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0,
               val > 0) %>%
        filter(var == target_solute) %>%
        select(datetime, val) %>%
        na.omit()

    # find acceptable years
    q_check <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        filter(ms_interp == 0) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(datetime, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)

    conc_check <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(date, origin = "usgs"),
               quart = quarter(date)) %>%
        group_by(water_year) %>%
        summarise(count = n_distinct(quart),
                  n = n()) %>%
        filter(n >= 4,
               count > 3)

    q_good_years <- q_check$water_year
    conc_good_years <- conc_check$water_year

    # 'good years' where Q and Chem data both meet min requirements
    good_years <- q_good_years[q_good_years %in% conc_good_years]
    n_yrs <- length(good_years)

    #join data and cut to good years
    daily_data_con <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean(val)) %>%
        mutate(site_code = !!site_code, var = 'con') %>%
        select(site_code, datetime = date, var, val)

    daily_data_q <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean(val)) %>%
        mutate(site_code = !!site_code, var = 'q_lps') %>%
        select(site_code, datetime = date, var, val)

    raw_data_full <- rbind(daily_data_con, daily_data_q) %>%
        pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
        mutate(wy = water_year(datetime, origin = 'usgs')) %>%
        filter(wy %in% good_years)

    ## k = 16
    ### Loop through good years #####
    for(k in 1:length(good_years)){

      writeLines(paste("site:", site_code,
                       'year:', good_years[k]))

        target_year <- as.numeric(as.character(good_years[k]))

        # calculate flag ratios to carry forward
        flag_df <- carry_flags(raw_q_df = raw_data_q,
                               raw_con_df = raw_data_con_in,
                               target_year = target_year,
                               target_solute = target_solute,
                               period = 'month')

        raw_data_target_year <- raw_data_full %>%
            mutate(wy = as.numeric(as.character(wy))) %>%
            filter(wy == target_year)

        # find months with at least 3 samples in them
        good_months <- raw_data_target_year %>%
            mutate(month = month(datetime)) %>%
            group_by(month) %>%
            tally(!is.na(con)) %>%
            filter(n > 2) %>%
            pull(month)

        # isolate to target year
        q_target_year <- raw_data_target_year %>%
            select(site_code, datetime, q_lps, wy)%>%
            filter(wy == target_year) %>%
            na.omit()

        con_target_year <- raw_data_target_year %>%
            select(site_code, datetime, con, wy) %>%
            na.omit()

        ### calculate monthly flux ######
        chem_df <- con_target_year
        q_df <- q_target_year

        #### calculate average ####
        flux_monthly_average <- raw_data_target_year %>%
            mutate(month = month(datetime)) %>%
            group_by(wy, month) %>%
            summarize(q_lps = mean(q_lps, na.rm = TRUE),
                      con = mean(con, na.rm = TRUE),
                      date = max(datetime)) %>%
            # multiply by seconds in a year, and divide my mg to kg conversion (1M)
            mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6,
                   date = substr(as.character(date),1,nchar(as.character(date))-3)) %>%
            ungroup() %>%
            select(date, flux)

        ##### error estimation ######
        theta <- function(x, raw_data_target_month, n_index){
            flux_monthly_average <- raw_data_target_month[n_index[x],] %>%
                summarize(q_lps = mean(q_lps, na.rm = TRUE),
                          con = mean(con, na.rm = TRUE)) %>%
                # multiply by seconds in a year, and divide my mg to kg conversion (1M)
                mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
                pull(flux)
            return(flux_monthly_average)
        }

        ###### initialize output for jackknife #####
        flux_monthly_average_jack <- tibble(month = as.integer(),
                                  se = as.numeric(),
                                  bias = as.numeric())

        ###### apply jackknife ######
        for(l in good_months){
            target_month <- good_months[l]

            raw_data_target_month <- raw_data_target_year %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            n_index <- raw_data_target_month %>%
                arrange(datetime) %>%
                mutate(index = row_number(datetime)) %>%
                select(con, index) %>%
                na.omit() %>%
                pull(index)

            n = length(n_index)

            month_jack <- jackknife(1:n,theta, raw_data_target_year, n_index)

            out_jack <- tibble(month = target_month,
                               se = month_jack$jack.se,
                               bias = month_jack$jack.bias)

            flux_monthly_average_jack <- rbind(flux_monthly_average_jack, out_jack)
        }

        #### calculate period weighted #####
        pw_con_df <- chem_df %>%
            mutate(month = month(datetime)) %>%
            filter(month %in% good_months) %>%
            select(-month)

        flux_monthly_pw <- calculate_pw(chem_df = pw_con_df, q_df,
                                       datecol = 'datetime', period = 'month')

        ##### error estimation ######
        theta <- function(x, pw_con_df_month, q_df, n_index){
            chem_df_jack <- pw_con_df_month[n_index[x],]

            flux_monthly_pw <- calculate_pw(chem_df = chem_df_jack, q_df = q_df,
                                            datecol = 'datetime')

            return(flux_monthly_pw)
        }

        ###### initialize output for jackknife #####
        flux_monthly_pw_jack <- tibble(month = as.integer(),
                                            se = as.numeric(),
                                            bias = as.numeric())

        ###### apply jackknife ######
        for(l in good_months){
            target_month <- good_months[l]

            pw_con_df_month <- pw_con_df %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            n_index <- pw_con_df_month  %>%
                arrange(datetime) %>%
                mutate(index = row_number(datetime)) %>%
                select(con, index) %>%
                na.omit() %>%
                pull(index)

            n = length(n_index)

            month_jack <- jackknife(1:n,theta, pw_con_df_month, q_df, n_index)

            out_jack <- tibble(month = target_month,
                               se = month_jack$jack.se,
                               bias = month_jack$jack.bias)

            flux_monthly_pw_jack <- rbind(flux_monthly_pw_jack, out_jack)
        }

        #### calculate beale ######
        beale_df <- chem_df %>%
            mutate(month = month(datetime)) %>%
            filter(month %in% good_months)

        flux_monthly_beale <- calculate_beale(chem_df = beale_df, q_df, datecol = 'datetime',
                                              period = 'month')

        ##### error estimation ######
        theta <- function(x, beale_df_month, q_df_month, n_index){
            beale_df_jack <- beale_df_month[n_index[x],]

            flux_monthly_beale <- calculate_pw(chem_df = beale_df_jack, q_df = q_df_month,
                                            datecol = 'datetime', period = 'month')

            return(flux_monthly_beale)
        }

        ###### initialize output for jackknife #####
        flux_monthly_beale_jack <- tibble(month = as.integer(),
                                       se = as.numeric(),
                                       bias = as.numeric())

        ###### apply jackknife ######
        for(l in good_months){
            target_month <- good_months[l]

            beale_df_month <- beale_df %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            q_df_month <- q_df %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            n_index <- beale_df_month  %>%
                arrange(datetime) %>%
                mutate(index = row_number(datetime)) %>%
                select(con, index) %>%
                na.omit() %>%
                pull(index)

            n = length(n_index)

            month_jack <- jackknife(1:n, theta, beale_df_month, q_df_month, n_index)

            out_jack <- tibble(month = target_month,
                               se = month_jack$jack.se,
                               bias = month_jack$jack.bias)

            flux_monthly_beale_jack <- rbind(flux_monthly_beale_jack, out_jack)
        }

        #### calculate rating #####
        flux_monthly_rating <- calculate_rating(chem_df, q_df, datecol = 'datetime',
                                                period = 'month')

        ##### error estimation ######
        theta <- function(x, chem_df_exclude, q_df_month, n_index){
            rating_df_jack <- chem_df_month[-n_index[-x],]

            flux_monthly_rating <- calculate_pw(chem_df = chem_df_month, q_df = q_df_month,
                                               datecol = 'datetime')

            return(flux_monthly_rating)
        }

        ###### initialize output for jackknife #####
        flux_monthly_rating_jack <- tibble(month = as.integer(),
                                          se = as.numeric(),
                                          bias = as.numeric())

        ###### apply jackknife ######
        for(l in good_months){
            target_month <- good_months[l]

            chem_df_month <- chem_df %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            q_df_month <- q_df %>%
                mutate(month = month(datetime)) %>%
                filter(month == target_month)

            n_index <- chem_df_month  %>%
                arrange(datetime) %>%
                mutate(index = row_number(datetime)) %>%
                select(con, index) %>%
                na.omit() %>%
                pull(index)

            n = length(n_index)

            month_jack <- jackknife(1:n, theta, chem_df_month, q_df_month, n_index)

            out_jack <- tibble(month = target_month,
                               se = month_jack$jack.se,
                               bias = month_jack$jack.bias)

            flux_monthly_rating_jack <- rbind(flux_monthly_rating_jack, out_jack)
        }

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = q_df,
                                                            datecol = 'datetime',
                                                            sitecol = 'site_code')

        # calculate annual flux from composite
        flux_monthly_comp <- calculate_composite_from_rating_filled_df(rating_filled_df, period = 'month') %>%#,
                                                                    #sitecol = 'site_code')
            mutate(date = substr(as.character(date),1,nchar(as.character(date))-3)) %>%
            ungroup() %>%
            select(date, flux)

        #### select MS favored ####
        paired_df <- q_df %>%
            full_join(chem_df, by = c('datetime', 'site_code', 'wy')) %>%
            na.omit() %>%
            filter(q_lps > 0,
                   is.finite(q_lps))

        q_log <- log10(paired_df$q_lps)
        c_log <- log10(paired_df$con)
        model_data <- tibble(c_log, q_log) %>%
            filter(is.finite(c_log),
                   is.finite(q_log))%>%
            na.omit()

        rating <- summary(lm(model_data$c_log ~ model_data$q_log, singular.ok = TRUE))

        r_squared <- rating$r.squared

        resid_acf <- abs(acf(rating$residuals, lag.max = 1, plot = FALSE)$acf[2])

        con_acf <- abs(acf(paired_df$con, lag.max = 1, plot = FALSE)$acf[2])

        # modified from figure 10 of Aulenbach et al 2016
        if(r_squared > 0.3){
            if(resid_acf > 0.2){
                ideal_method <- 'composite'
            }else{
                ideal_method <- 'rating'
            }
        }else{
            if(con_acf > 0.20){
                ideal_method <- 'pw'
            }else{
                ideal_method <- 'average'
            }
        }

        #### congeal fluxes ####
        target_year_out <- mutate(flux_monthly_average, method = 'average') %>%
            rbind(mutate(flux_monthly_pw, method = 'pw'),
                  mutate(flux_monthly_beale, method = 'beale'),
                  mutate(flux_monthly_rating, method = 'rating'),
                  mutate(flux_monthly_comp, method = 'composite')) %>%
            mutate(site_cod = !! site_code,
                   var = !!target_solute,
                   ms_recommended = ifelse(method == !!ideal_method, 1, 0),
                   good = ifelse(month(as_date(date)) %in% good_months, 1, 0)) %>%
            filter(good == 1) %>%
            select(-good)

        out_frame <- rbind(out_frame, target_year_out)

        } # end year loop
    } # end solute loop

    directory <- glue(data_dir,'stream_flux/')
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)
write_feather(out_frame, file_path)
} # end site loop

## w4df <- read_feather('data/ms/hbef/stream_flux/w4.feather')
