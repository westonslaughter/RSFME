library(tidyverse)
library(feather)
library(here)
library(glue)
library(lubridate)
library(EGRET)
library(macrosheds)
library(foreach)
library(doParallel)
library(bootstrap)
library(lfstat)

ncores <- detectCores()

if(.Platform$OS.type == "windows"){
    cl <- makeCluster(ncores, type = 'PSOCK')
} else {
    cl <- makeCluster(ncores, type = 'FORK')
}

registerDoParallel(cl)

source('source/helper_functions.R')
source('source/egret_overwrites.R')
source('ms_overwrites.R')
source('source/flux_methods.R')
source('source/usgs_helpers.R')

# ng

# site_files  <- list.files('streamlined/data/ms/hbef/discharge', recursive = F)
# site_info  <- read_csv(here('streamlined/data/site/ms_site_info.csv'))
## var_info <- nick/file/path

# ws
# data_dir <- here('/home/ws184/science/macrosheds/portal/data/')
# ms_root <- here('/home/ws184/science/macrosheds/portal/data/')
data_dir <- here('/home/ws184/science/macrosheds/RSFME/data/ms')
ms_root <- here('/home/ws184/science/macrosheds/RSFME/data/ms')
## run below if you do not already have macrosheds core data and catalogs
## set path to ms data
# ms_root <- 'data/lter/'
# library(remotes)
# remotes::install_github('https://github.com/MacroSHEDS/macrosheds.git')
# data_dir <-  ms_download_core_data(ms_root, 'all')

# site_files  <- list.files(here(data_dir, 'discharge'), recursive = F)
site_info <- ms_download_site_data()
var_info <-  ms_download_variables()

logfile <- '~/log.txt' #log to a file (necessary for receiving output from parallel processes)
#logfile = stdout() #log to console as usual

# set wether you are running with or without uncertainty estimation (LOOCV)
error_estimation = TRUE

## initializing loop output doesn't work in a parallel framework, since the initialized object
## is only accessible to the parent process
#   # df to populate with annual flux values by method
#   out_frame <- tibble(wy = as.integer(),
#                       site_code = as.character(),
#                       var = as.character(),
#                       val = as.numeric(),
#                       method = as.character(),
#                       ms_reccomended = as.integer(),
#                       ms_interp_ratio = as.numeric(),
#                       ms_status_ratio = as.numeric(),
#                       ms_missing_ratio = as.numeric())

## i = 2
## i = 3
# Loop through sites #####
#for(i in 1:length(site_files)){
domains <- list.files(file.path(data_dir))
data_dirs <- file.path(ms_root, 
                       # site_dmn_nwk$network, 
                       domains)
site_dmn_nwk <- site_info %>% 
  # filter(network == 'lter') %>% 
  filter(domain %in% domains) %>%
  select(network, domain, site_code)

site_files <- site_dmn_nwk %>% 
  # filter(domain %in% domains) %>%
  pull(site_code) 

out_frame <- foreach(i = 21:length(site_files), .combine = bind_rows) %do% {

    # site_file <- site_files[i]
    # site_code <- strsplit(site_file, split = '.feather')[[1]]
    site_code <- site_files[i]
    this_site_dmn_nwk <- site_dmn_nwk %>% filter(site_code == !!site_code) 
    # create this network-domain data fp
    # data_dir <- data_dirs[i]
    nwk = this_site_dmn_nwk$network
    dmn = this_site_dmn_nwk$domain
    data_dir <- file.path(ms_root, 
                           # nwk,
                           dmn)
    
    writeLines(glue('running flux calcs for \n network: {nwk} \n domain: {dmn}\n site: {site_code}', 
                    nwk = nwk,
                    dmn = dmn, 
                    site_code = site_code))
    
    area <- site_info %>%
        filter(site_code == !!site_code) %>%
        pull(ws_area_ha)

    lat <- site_info %>%
        filter(site_code == !!site_code) %>%
        select(any_of(c('Y', 'latitude'))) %>%
        pull()

    long <- site_info %>%
        filter(site_code == !!site_code) %>%
        select(any_of(c('X', 'longitude'))) %>%
        pull()

    # read in chemistry data
    raw_data_con_in <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0)

    # read in discharge data
    raw_data_q <- read_feather(here(glue(data_dir, '/discharge/', site_code, '.feather')))

    # errors
    raw_data_q$val = errors::set_errors(raw_data_q$val, raw_data_q$val_err)
    raw_data_q$val_err = NULL

    # initialize next loop
    solutes <- raw_data_con_in %>%
        ## switch to only Nitrate for now
        #filter(var == 'GN_NO3_N') %>%
        ## filter(!str_detect(var, 'Charge'),
        ##        !str_detect(var, 'temp'),
        ##        !str_detect(var, 'pH')) %>%
        select(var) %>%
        unique() %>%
        pull(var)

  writeLines(paste("FLUX CALCS:", site_code), con = logfile)

  ## Loop through solutes at site #####
  ## j = 1
  ## j = 20

  foreach(j = 1:length(solutes), .combine = bind_rows) %do% {

    writeLines(paste("site:", site_code,
                     "var:", solutes[j]),
               con = logfile)

    #set to target solute
    target_solute <- solutes[j]

    # convert all solutes to mg/L
    solute_name <- ms_drop_var_prefix(target_solute)
    solute_default_unit <- var_info[var_info$variable_code == solute_name,] %>%
      filter(variable_type %in% c('chem_discrete', 'chem_mix')) %>%
      ## filter(chem_category == 'stream_conc') %>%
      pull(unit)

    # read out conversions
    writeLines(paste("\n  unit conversion for", target_solute,
                     '\n    converting', solute_default_unit, 'to grams per liter (g/L)\n'),
               con = logfile)

    raw_data_con <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
              filter(ms_interp == 0,
                     val > 0) %>%
              filter(var == target_solute) %>%
              # conver units from macrosheds default to g/L
              ms_conversions(convert_units_from = tolower(solute_default_unit),
                                        convert_units_to = "mg/l",
                             macrosheds_root = ms_root) %>%
              select(datetime, val, val_err) %>%
              tidyr::drop_na(datetime, val)

    # errors
    raw_data_con$val = errors::set_errors(raw_data_con$val, raw_data_con$val_err)
    raw_data_con$val_err = NULL

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

    # TODO: calculate BREAK input to egret

    #join data and cut to good years
    daily_data_con <- raw_data_con %>%
      mutate(date = date(datetime)) %>%
      group_by(date) %>%
      summarize(val = mean_or_x(val)) %>%
      # this is the step where concentration value errors turn to NA
        mutate(site_code = !!site_code, var = 'con') %>%
        select(site_code, datetime = date, var, val)

    daily_data_q <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean_or_x(val)) %>%
      # this is the step where discharge value errors turn to NA
        mutate(site_code = !!site_code, var = 'q_lps') %>%
        select(site_code, datetime = date, var, val)

    q_df <- daily_data_q %>%
      pivot_wider(names_from = var,
                  values_from = val)

    raw_data_full <- rbind(daily_data_con, daily_data_q) %>%
        pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
        mutate(wy = water_year(datetime, origin = 'usgs')) %>%
        filter(wy %in% good_years)

    ## big nope on this i think
    ## q_full <- raw_data_full %>%
    ##       mutate(wy = as.numeric(as.character(wy))) %>%
    ##         select(site_code, datetime, q_lps, wy)%>%
    ##         na.omit()

    con_full <- raw_data_full %>%
          mutate(wy = as.numeric(as.character(wy))) %>%
            select(site_code, datetime, con, wy) %>%
            ## filter(wy < 1975) %>%
            na.omit()

    #### calculate WRTDS ######
    flux_annual_wrtds <- calculate_wrtds(
          chem_df = con_full,
          q_df = q_df,
          ws_size = area,
          lat = lat,
          long = long,
          datecol = 'datetime',
          agg = 'annual',
          minNumObs = 100,
          minNumUncen = 50
         )

    ## write_feather(raw_data_full, "data/ms/hbef/true/w3_chem_samples.feather")
    ## k = 1
    ### Loop through good years #####

    # for(k in 1:length(good_years)){
    foreach(k = 1:length(good_years), .combine = bind_rows) %dopar% {

      writeLines(paste("site:", site_code,
                       'year:', good_years[k]),
                 con = logfile)

        target_year <- as.numeric(as.character(good_years[k]))

        # calculate flag ratios to carry forward
        flag_df <- carry_flags(raw_q_df = raw_data_q,
                               raw_con_df = raw_data_con_in,
                               target_year = target_year,
                               target_solute = target_solute,
                               period = 'annual')

        raw_data_target_year <- raw_data_full %>%
            mutate(wy = as.numeric(as.character(wy))) %>%
            filter(wy == target_year)

        q_target_year <- raw_data_target_year %>%
            select(site_code, datetime, q_lps, wy)%>%
            na.omit()

        con_target_year <- raw_data_target_year %>%
            select(site_code, datetime, con, wy) %>%
            na.omit()

        ##### calculate annual flux ######
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

        ###### error estimation ######

        n_index <- raw_data_target_year %>%
            arrange(datetime) %>%
            mutate(index = row_number(datetime)) %>%
            select(con, index) %>%
            na.omit() %>%
            pull(index)

        n = length(n_index)

        theta <- function(x, n_index, raw_data_target_year){

            flux_annual_average <- raw_data_target_year[n_index[x],] %>%
                group_by(wy) %>%
                summarize(q_lps = mean(q_lps, na.rm = TRUE),
                          con = mean(con, na.rm = TRUE)) %>%
                # multiply by seconds in a year, and divide my mg to kg conversion (1M)
                mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
                pull(flux)
            return(flux_annual_average)
        }
        
        if(error_estimation == TRUE) {
          flux_annual_average_jack <- jackknife(1:n, theta, n_index, raw_data_target_year)
        }


        #### calculate period weighted #####
        flux_annual_pw <- calculate_pw(chem_df, q_df, datecol = 'datetime')

        ###### error estimation ######
        n = nrow(chem_df)
        theta <- function(x, chem_df, q_df){calculate_pw(chem_df[x,], q_df, datecol = 'datetime') }
        
        if(error_estimation == TRUE) {
          flux_annual_pw_jack <- jackknife(0:n, theta, chem_df, q_df)
        }


        #### calculate beale ######
        flux_annual_beale <- calculate_beale(chem_df, q_df, datecol = 'datetime')

        ###### error estimation ######
        n = nrow(chem_df)
        theta <- function(x, chem_df, q_df){calculate_beale(chem_df[x,], q_df, datecol = 'datetime') }
        
        if(error_estimation == TRUE) {
          flux_annual_beale_jack <- jackknife(1:n,theta, chem_df, q_df)
        }

        #### calculate rating #####
        flux_annual_rating <- calculate_rating(chem_df, q_df, datecol = 'datetime')

        ###### error estimation ######
        n = nrow(chem_df)
        theta <- function(x, chem_df, q_df){calculate_rating(chem_df[x,], q_df, datecol = 'datetime') }
        
        if(error_estimation == TRUE) {
          flux_annual_rating_jack <- jackknife(1:n,theta, chem_df, q_df)
        }

        #### calculate WRTDS ######

        # put agg here
        #### calculate composite ######
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = q_df,
                                                            datecol = 'datetime',
                                                            sitecol = 'site_code')

        # calculate annual flux from composite
        flux_annual_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)

        ###### error estimation ######
        n = nrow(chem_df)
        theta <- function(x, chem_df, q_df){
            #### calculate composite ######

            rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df[-x,],
                                                                q_df = q_df,
                                                                datecol = 'datetime',
                                                                sitecol = 'site_code')
            #if(!is.logical(rating_filled_df)){
            out <- calculate_composite_from_rating_filled_df(rating_filled_df)
            return(out$flux[1])
            #}else{return(NA)}

        }
        
        if(error_estimation == TRUE) {
          flux_annual_comp_jack <- jackknife(1:n,theta, chem_df, q_df)
        }

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

        # placeholder for testing
        flux_annual_wrtds = 99999
        
        if(error_estimation == TRUE) {
          flux_annual_wrtds_jack = tibble(jack.se = 9999, jack.bias = 9999)
        }

        #### congeal fluxes ####
        target_year_out <- tibble(wy = as.character(target_year),
                                  val = c(flux_annual_average,
                                          flux_annual_pw,
                                          flux_annual_beale,
                                          flux_annual_rating,
                                          ## flux_annual_wrtds,
                                          flux_annual_comp$flux[1]),
                            site_code = !!site_code,
                            var = !!target_solute,
                            method = c('average', 'pw', 'beale', 'rating', 'composite')) %>%
            mutate(ms_recommended = ifelse(method == !!ideal_method, 1, 0))

        #### congeal SE ####
        if(error_estimation == TRUE) {
          target_year_out_jack <- tibble(wy = as.character(target_year),
                                  se = c(flux_annual_average_jack$jack.se,
                                          flux_annual_pw_jack$jack.se,
                                          flux_annual_beale_jack$jack.se,
                                          flux_annual_rating_jack$jack.se,
                                          flux_annual_wrtds_jack$jack.se,
                                          flux_annual_comp_jack$jack.se),
                                  bias = c(flux_annual_average_jack$jack.bias,
                                         flux_annual_pw_jack$jack.bias,
                                         flux_annual_beale_jack$jack.bias,
                                         flux_annual_rating_jack$jack.bias,
                                         flux_annual_wrtds_jack$jack.bias,
                                         flux_annual_comp_jack$jack.bias),
                                  method = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite'))
        }

        #### OUTPUT ####
        
        if(error_estimation == TRUE) {
            target_year_out_combined <- target_year_out %>%
              left_join(., target_year_out_jack, by = 'method')
            
            # out_frame <- rbind(out_frame, target_year_out_combined)
            return(target_year_out_combined)
        } else {
          return(target_year_out)
        }

        } # end year loop

        wrtds_out <- flux_annual_wrtds %>%
          filter(wy %in% good_years) %>%
          rename(val = flux) %>%
          mutate(site_code = site_code,
                 var = solutes[j],
                 method = 'wrtds',
                 ms_recommended = 0)

        out_frame <- rbind(out_frame, wrtds_out) %>%
          arrange(desc(method), desc(wy))
    } # end solute loop

    directory <- file.path(data_dir,'stream_flux')
    
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    
    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)

  write_feather(out_frame, file_path)
} # end site loop

stopCluster(cl)
