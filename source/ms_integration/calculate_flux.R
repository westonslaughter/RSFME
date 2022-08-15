library(tidyverse)
library(feather)
library(here)
library(glue)
library(lubridate)
library(EGRET)
library(macrosheds)

source('source/helper_functions.R')
source('source/egret_overwrites.R')
source('source/flux_methods.R')
source('source/usgs_helpers.R')

data_dir <- 'data/ms/hbef/'
site_files  <- list.files('data/ms/hbef/discharge', recursive = F)
site_info  <- read_csv('data/site/ms_site_info.csv')

# df to populate with annual flux values by method
out_frame <- tibble(wy = as.integer(), 
                    site_code = as.character(),
                    var = as.character(),
                    val = as.numeric(), 
                    method = as.character(),
                    ms_reccomended = as.integer())
## i = 3
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
    raw_data_con <- read_feather(here(glue(data_dir, 'stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0)

    # read in discharge data
    raw_data_q <- read_feather(here(glue(data_dir, 'discharge/', site_code, '.feather')))

    # initialize next loop
    solutes <- raw_data_con %>%
        ## switch to only Nitrate for now
        filter(var == 'GN_NO3_N') %>%
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

    raw_data_con <- read_feather(here(glue(data_dir, 'stream_chemistry/', site_code, '.feather'))) %>%
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
        mutate(water_year = water_year(date, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 4)
    
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

    ## write_feather(raw_data_full, "data/ms/hbef/true/w3_chem_samples.feather")
    ## k = 1
    ### Loop through good years #####
    ## for(k in 42:47) {
    for(k in 1:length(good_years)){

      writeLines(paste("site:", site_code,
                       'year:', good_years[k]))

        target_year <- as.numeric(as.character(good_years[k]))
        
        raw_data_target_year <- raw_data_full %>%
            mutate(wy = as.numeric(as.character(wy))) %>%
            filter(wy == target_year)
        
        q_target_year <- raw_data_target_year %>%
            select(site_code, datetime, q_lps, wy)%>%
            na.omit()
        
        con_target_year <- raw_data_target_year %>%
            select(site_code, datetime, con, wy) %>%
            na.omit()
    
        ### calculate annual flux ######
        chem_df <- con_target_year
        q_df <- q_target_year
        
        #### calculate average ####
        flux_annual_average <- raw_data_target_year %>%
            group_by(wy) %>%
            summarize(q_lps = mean(q_lps, na.rm = TRUE),
                      con = mean(con, na.rm = TRUE)) %>%
            # multiply by seconds in a year, and divide my mg to kg conversion (1M)
            mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
            pull(flux)
            
        
        #### calculate period weighted #####
        flux_annual_pw <- calculate_pw(chem_df, q_df, datecol = 'datetime')
        
        #### calculate beale ######
        flux_annual_beale <- calculate_beale(chem_df, q_df, datecol = 'datetime')

        #### calculate rating #####

        flux_annual_rating <- calculate_rating(chem_df, q_df, datecol = 'datetime')

        #### calculate WRTDS ######
        flux_annual_wrtds <- calculate_wrtds(
          chem_df = chem_df,
          q_df = q_df,
          ws_size = area,
          lat = lat,
          long = long,
          ## datamode = 'ms',
          datecol = 'datetime')

        # get daily WRTDS
        ## flux_daily <- adapt_ms_egret(
        ##   chem_df = chem_df,
        ##   q_df = q_df,
        ##   ws_size = area,
        ##   lat = lat,
        ##   long = long,
        ##   datecol = 'datetime')

        filepath = paste0('data/ms/hbef/true/w3_dailyWRTDS_', good_years[k], '.feather')
        write_feather(flux_daily$Daily, filepath)

        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = q_df,
                                                            datecol = 'datetime',
                                                            sitecol = 'site_code')
        
        # calculate annual flux from composite
        flux_annual_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)
        
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
        target_year_out <- tibble(wy = as.character(target_year), 
                                  val = c(flux_annual_average,
                                          flux_annual_pw,
                                          flux_annual_beale,
                                          flux_annual_rating,
                                          flux_annual_wrtds,
                                          flux_annual_comp$flux[1]),
                            site_code = !!site_code, 
                            var = !!target_solute,
                            method = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite')) %>%
            mutate(ms_recommended = ifelse(method == !!ideal_method, 1, 0))
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
