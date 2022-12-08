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
library(RiverLoad)

# clean processing environment
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

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
# data_dir <- here('streamlined/data/ms/hbef/')
## site_files  <- list.files('streamlined/data/ms/hbef/discharge', recursive = F)
# site_info  <- read_csv(here('streamlined/data/site/ms_site_info.csv'))
## var_info <- nick/file/path

## run below if you do not already have macrosheds core data and catalogs
## set path to ms data
# data_dir <-  ms_download_core_data(ms_root)

# ws
ms_root <- 'data/ms'
site_info <- ms_download_site_data()
var_info <-  ms_download_variables()

logfile <- '~/log.txt' #log to a file (necessary for receiving output from parallel processes)
#logfile = stdout() #log to console as usual

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


# Loop through domains #####
domains <- list.dirs(file.path(ms_root), recursive = FALSE)
domain_names <- basename(list.dirs(file.path(ms_root), recursive = FALSE))

# TODO: skipping any domains this run?
# domain_filter = c('baltimore', 'fernow')
all_stream_chem_sites <- c()

for(d in 1:length(domains)) {
  dmn = domains[d]
  directory <- glue(file.path(dmn, 'stream_flux'))

  stream_chem_dir <- glue(file.path(dmn, 'stream_chemistry'))
  stream_chem_files <- list.files(stream_chem_dir, pattern = ".feather")
  sites <- tools::file_path_sans_ext(stream_chem_files)

  all_stream_chem_sites <- c(all_stream_chem_sites, sites)
  if(dir.exists(directory)){
    print(glue('{dmn} stream flux directory found'))
    files_flux <- list.files(directory, pattern = ".feather")
    sites_flux <- tools::file_path_sans_ext(files_flux)
  }

  need_sites <- sites[!sites %in% sites_flux]
}


# clear stream flux?
# for(d in 1:length(domains)) {
#   dmn = domains[d]
#   directory <- glue(file.path(dmn, 'stream_flux'))
#   directory_ftr <- glue(file.path(dmn, 'stream_flux', '*.feather'))
#   directory_txt <- glue(file.path(dmn, 'stream_flux', '*.txt'))
# 
#   if(dir.exists(directory)){
#     print(glue('deleting {dmn} stream flux feather amd text files'))
#     unlink(directory_ftr)
#     unlink(directory_txt)
#   }
# }

main_frame <- site_info %>% 
  filter(domain %in% domain_names,
         site_code %in% all_stream_chem_sites
         ) %>%
  select(domain, site_code) %>%
  mutate(
    flux = FALSE
  )
  
# index=2
# for(index in 4:length(domains)){
foreach(index = 1:length(domains), .packages = c('tidyverse', 'feather', 'glue', 'lubridate',
                                                 'EGRET', 'macrosheds', 'foreach', 'doParallel',
                                                 'bootstrap', 'lfstat', 'RiverLoad'), .verbose = TRUE) %do% {
                                                   
    # want to use LOOCV or not?
    error_estimation = FALSE
    
    # set warnings to print as they occur
    options(warn = 1)
    
    # load list of flux convertible vars
    ms_flux_vars <- ms_download_variables() %>%
      filter(flux_convertible == 1) %>%
      pull(variable_code)
                                                   
    source('source/helper_functions.R')
    source('source/egret_overwrites.R')
    source('ms_overwrites.R')
    source('source/flux_methods.R')
    source('source/usgs_helpers.R')

    dmn = domains[index]
    dmn_nm = domain_names[index]
      
    writeLines(glue('\n\n ----- \n\n FLUX CALCULATION for data in {dmn_nm} \n\n-----\n\n', dmn = dmn))
    data_dir = dmn
    
    site_files = list.files(file.path(dmn, 'discharge'), full.names = FALSE)
    
    # create domain stream flux directory
    directory <- glue(file.path(data_dir, 'stream_flux'))
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    
    # cat(glue(''), file = log_fp, sep = "\n", append = TRUE)
    
      # i=1
    # for(i in 1:length(site_files)){
      # Loop through sites #####
    
    foreach(i = 1:length(site_files), .packages = c('tidyverse', 'feather', 'glue', 'lubridate',
                                                    'EGRET', 'macrosheds', 'foreach', 'doParallel',
                                                    'bootstrap', 'lfstat', 'RiverLoad'), 
            .verbose = TRUE, .errorhandling = "pass") %dopar% {
  
      site_file <- site_files[i]
      site_code <- strsplit(site_file, split = '.feather')[[1]]
  
      # add log file into this directory, named by site
      log_fn <- glue('{dmn_nm}_{site_code}__annual_flux_log.txt')
      log_fp <- file.path(directory, log_fn) 
      
      start.tm <- Sys.time()
      start.tz <- Sys.timezone()
      start.dt <- Sys.Date()
      
      cat(glue('{dmn_nm} annual flux log\n'), file = log_fp, sep = "\n")
      cat(glue('time: {start.tm} {start.tz}\ndate:{start.dt}'), file = log_fp, sep = "\n", append = TRUE)
      # NOTE: logfile standard, I will add two hyphens "--" in every log statement for every for loop level deep we go
      #       starting here, so:
      cat(glue('-- {site_code}'), file = log_fp, sep = "\n", append = TRUE)
      
      cat(glue('-- STEP_1: read in site info and all chem and discharge data'), file = log_fp, sep = "\n", append = TRUE)
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
      raw_data_con_in <- tryCatch(
        expr = {
          raw_data_con_in <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
              filter(ms_interp == 0,
                   # NOTE: dropping all vars that are not marked as 'flux convertible'
                   ms_drop_var_prefix(var) %in% ms_flux_vars)
        },
        error = function(e) {
          cat(glue('-- STEP_1_ERROR: {site_code} chem data read-in failed'), file = log_fp, sep = "\n", append = TRUE)
          NULL
        }
      )
      
      if(is.null(raw_data_con_in)) {
        print('raw chem failed, return NULL')
        stop("SKIP")
      }
      
      # read in discharge data
      raw_data_q <- tryCatch(
        expr = {
          read_feather(here(glue(data_dir, '/discharge/', site_code, '.feather')))
        },
        error = function(e) {
          cat(glue('-- STEP_1_ERROR: {site_code} chemdischarge data read-in failed'), file = log_fp, sep = "\n", append = TRUE)
        }
      )
      
      if(is.null(raw_data_q)) {
        print('raw Q failed')
        stop("SKIP")
      }
      
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
  
    out_frame <- foreach(j = 1:length(solutes), .combine = bind_rows, .errorhandling = "pass") %do% {
  
      
      writeLines(paste("site:", site_code,
                       "var:", solutes[j]))
  
      #set to target solute
      target_solute <- solutes[j]
  
      
      # two more hypehns on every line, one level deeper
      cat(glue('---- {target_solute}'), file = log_fp, sep = "\n", append = TRUE)
      cat(glue('---- STEP_2: solute data prep'), file = log_fp, sep = "\n", append = TRUE)
      
      # convert all solutes to mg/L
      solute_name <- ms_drop_var_prefix(target_solute)
      solute_default_unit <- var_info[var_info$variable_code == solute_name,] %>%
        filter(variable_type %in% c('chem_discrete', 'chem_mix')) %>%
        ## filter(chem_category == 'stream_conc') %>%
        pull(unit)
  
      # if no defaultunit, skip solute
      if(length(solute_default_unit) == 0) {
        cat(glue('---- STEP_2_ERROR: solute failed to read in default unit'), file = log_fp, sep = "\n", append = TRUE)
        stop("skip")
      }
      
      # read out conversions
      writeLines(paste("\n  unit conversion for", target_solute,
                       '\n    converting', solute_default_unit, 'to grams per liter (g/L)\n'))
  
      raw_data_con <- tryCatch(
        expr = {
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
        },
        error = function(e) {
          cat(glue('---- STEP_2_ERROR: solute failed to read in data'), file = log_fp, sep = "\n", append = TRUE)
          NULL
        }
        )
      
      if(is.null(raw_data_con)) {
        print("chemistry data read failed")
        # next
        stop("SKIP")
      }
  
      # errors
      raw_data_con$val = errors::set_errors(raw_data_con$val, raw_data_con$val_err)
      raw_data_con$val_err = NULL
  
      cat(glue('---- STEP_3: data checks'), file = log_fp, sep = "\n", append = TRUE)
      # set calc minimums
      # Q
      q_days_min = 311
      # C
      c_count_quarters_min = 3 # number of distinct seasons
      c_n_sample_min = 4 # number of samples
      
      # find acceptable years
      q_check <- raw_data_q %>%
          mutate(date = date(datetime)) %>%
          filter(ms_interp == 0) %>%
          distinct(., date, .keep_all = TRUE) %>%
          mutate(water_year = water_year(datetime, origin = "usgs")) %>%
          group_by(water_year) %>%
          summarise(n = n()) %>%
          # filter to above minimum stndards
          filter(n >= q_days_min) %>%
          mutate(
           water_year = droplevels(water_year)
          )
  
      conc_check <- raw_data_con %>%
          mutate(date = date(datetime)) %>%
          distinct(., date, .keep_all = TRUE) %>%
          mutate(water_year = water_year(date, origin = "usgs"),
                 quart = quarter(date)) %>%
          group_by(water_year) %>%
          summarise(count = n_distinct(quart),
                    n = n()) %>%
          # filter to above minimum stndards
          filter(n >= c_n_sample_min,
                 count > c_count_quarters_min) %>%
        mutate(
          water_year = droplevels(water_year)
        )
      
      q_good_years <- q_check$water_year
      conc_good_years <- conc_check$water_year
  
      # 'good years' where Q and Chem data both meet min requirements
      good_years <- q_good_years[q_good_years %in% conc_good_years] %>% droplevels()
      n_yrs <- length(good_years)
      
      
      if(is.null(good_years)) {
        print("good years object is NULL")
        cat(glue('---- STEP_3_ERROR: good years is NULL '), file = log_fp, sep = "\n", append = TRUE)
        
        print('good years NULL, SKIP')
        stop()
        # next
      } else if(length(good_years) == 0) {
          cat(glue('---- STEP_3_ERROR: no good years '), file = log_fp, sep = "\n", append = TRUE)
          print(glue('\n----    Alert: data at {dmn}, site: {site_code}, solute: {target_solute},',
                       'data has no "good years" of raw data, where:\n',
                       '            Q has samples in at least {q_days_min} days \n',
                       '            C has samples in at least {c_count_quarters_min} seasons',
                       ' and at least {c_n_sample_min} total samples',
                       '. jumping to next solute'))
          
        print('no good years, SKIP')
        stop()
        # next
      }
  
      # if good years pass checks, log years and number of
      cat(glue('---- (STEP_3) good years: {good_years} '), file = log_fp, sep = "\n", append = TRUE)
      cat(glue('---- (STEP_3) total good years: {n_yrs} '), file = log_fp, sep = "\n", append = TRUE)
      
      # TODO: calculate BREAK input to egret
      #join data and cut to good years
      daily_data_con <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean_or_x(val)) %>%
        # this is the step where concentration value errors turn to NA
          mutate(
            site_code = !!site_code, 
            # why turn all vars to 'con'?
            var = 'con'
            ) %>%
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
  
      # TODO: must add a record of ALL WRTDS arguments to official logfile. 
      #       tricky part is that most important things is if minNumObs or Uncen
      #       changes, and this happens inside func. may have to add a "con ="
      # arg to flux_annual_wrtds() to set the logfile to write this data to
    
      cat(glue('---- STEP_4: calculate flux '), file = log_fp, sep = "\n", append = TRUE)
      
      #### calculate WRTDS ######
      tryCatch(
        expr = {
          cat(glue('---- STEP_4_1: calculate WRTDS '), file = log_fp, sep = "\n", append = TRUE)
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
        },
        error = function(e) {
          cat(glue('---- STEP_4_1_ERROR: calculate WRTDS failed '), file = log_fp, sep = "\n", append = TRUE)
          writeLines(glue("---- warning: WRTDS failed, set to double NA ----"))
          flux_annual_wrtds <- c('wy' = NA, 'flux' = NA)
        }
      )
      writeLines(glue("---- this WRTDS COMPLETE ----"))
      # TODO: LOG stats about WRTDS, and most importantly every single argument value should be 
      #       supplied to log file, esp minNum args, which are dynamic and effect the model
  
      ### Loop through good years #####
      # for(k in 0:length(good_years)){
      year_out_frame <- foreach(k = 1:length(good_years), .combine = bind_rows) %do% {
        this_yr = good_years[k]
        
        # parallel info (using pid as proxy for CPU id within given site-solute flux calc)
        rpid <- Sys.getpid()
        
        cat(glue('------ (STEP_4) {this_yr} flux calcs (paralell) '), file = log_fp, sep = "\n", append = TRUE)
        cat(glue('------ (STEP_4)       session: {rpid} '), file = log_fp, sep = "\n", append = TRUE)
        
        year_fail <- NULL
        tryCatch(
          expr = {
          writeLines(glue("---- year loop ----"))
          
    
            target_year <- as.numeric(as.character(good_years[k]))
            
            writeLines(paste("site:", site_code,
                           "solute:", target_solute,
                           'year:', target_year),
                     con = logfile)
    
            writeLines(paste("site:", site_code,
                           "solute:", target_solute,
                           'year:', target_year))
    
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
    
            # TODO: log each flux calc
            
            #### calculate average ####
            cat(glue('------ STEP_4_2: calculate average '), file = log_fp, sep = "\n", append = TRUE)
            flux_annual_average <- raw_data_target_year %>%
                group_by(wy) %>%
                summarize(q_lps = mean(q_lps, na.rm = TRUE),
                          con = mean(con, na.rm = TRUE)) %>%
                # multiply by seconds in a year, and divide my mg to kg conversion (1M)
                mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
                pull(flux)
    
            ###### error estimation ######
            if(error_estimation == TRUE) {
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
            flux_annual_average_jack <- jackknife(1:n,theta, n_index, raw_data_target_year)
            }
    
    
            #### calculate period weighted #####
            cat(glue('------ STEP_4_3: calculate period weighted '), file = log_fp, sep = "\n", append = TRUE)
            flux_annual_pw <- calculate_pw(chem_df, q_df, datecol = 'datetime')
    
            ###### error estimation ######
            if(error_estimation == TRUE) {
            n = nrow(chem_df)
            theta <- function(x, chem_df, q_df){calculate_pw(chem_df[x,], q_df, datecol = 'datetime') }
            flux_annual_pw_jack <- jackknife(1:n,theta, chem_df, q_df)
            }
    
            #### calculate beale ######
            cat(glue('------ STEP_4_4: calculate beale '), file = log_fp, sep = "\n", append = TRUE)
            flux_annual_beale <- calculate_beale(chem_df, q_df, datecol = 'datetime')
    
            ###### error estimation ######
            if(error_estimation == TRUE) {
            n = nrow(chem_df)
            theta <- function(x, chem_df, q_df){calculate_beale(chem_df[x,], q_df, datecol = 'datetime') }
            flux_annual_beale_jack <- jackknife(1:n,theta, chem_df, q_df)
            }
    
            #### calculate rating #####
            cat(glue('------ STEP_4_5: calculate rating'), file = log_fp, sep = "\n", append = TRUE)
            flux_annual_rating <- calculate_rating(chem_df, q_df, datecol = 'datetime')
    
            ###### error estimation ######
            if(error_estimation == TRUE) {
            n = nrow(chem_df)
            theta <- function(x, chem_df, q_df){calculate_rating(chem_df[x,], q_df, datecol = 'datetime') }
            flux_annual_rating_jack <- jackknife(1:n,theta, chem_df, q_df)
            }
            
            #### calculate composite ######
            cat(glue('------ STEP_4_6: calculate composite'), file = log_fp, sep = "\n", append = TRUE)
            rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                                q_df = q_df,
                                                                datecol = 'datetime',
                                                                sitecol = 'site_code')
    
            # calculate annual flux from composite
            flux_annual_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)
    
            ###### error estimation ######
            if(error_estimation == TRUE) {
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
            flux_annual_comp_jack <- jackknife(1:n,theta, chem_df, q_df)
            }
    
            cat(glue('------ STEP_5: flux meta operations (recomendation) '), file = log_fp, sep = "\n", append = TRUE)
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
            tryCatch(
              expr = {
                if(is.na(r_squared)) {
                  ideal_method <- NA
                  writeLines(glue("---- warning: r_squared failed to compute"), con = logfile)
                  writeLines(glue("---- warning: r_squared failed to compute"))
                  } else {
                    # browser()
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
                }
              },
              error = function(e) {
                ideal_method <- NA
                writeLines(glue("---- warning: Aulenbach conditional error"), con = logfile)
                writeLines(glue("---- warning: Aulenbach conditional error"))
              }
            )
    
            # placeholder for testing
            # if(error_estimation == TRUE) {
                # flux_annual_wrtds = 99998
            #   flux_annual_wrtds_jack = tibble(jack.se = 9999, jack.bias = 9999)
            # }
    
            cat(glue('------ STEP_6: combine fluxes '), file = log_fp, sep = "\n", append = TRUE)
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
            
            writeLines(glue("---- complete: {target_year}"))
            
            #### add to output ####
            if(error_estimation == TRUE) {
                target_year_out_combined <- target_year_out %>%
                    left_join(., target_year_out_jack, by = 'method')
          
                # out_frame <- rbind(out_frame, target_year_out_combined)
                return(target_year_out_combined)
              } else {
               cat(glue('------ STEP_6_SUCCESS: {this_yr} complete'), file = log_fp, sep = "\n", append = TRUE)
               return(target_year_out)
              }
          },
          error = function(e) {
            year_fail <- "error"
          })
        
        if(!is.null(year_fail)) {
          writeLines(glue("year: {this_yr} FAILED, moving on"))
          cat(glue('------ STEP_6_ERROR: {this_yr} failed'), file = log_fp, sep = "\n", append = TRUE)
          next
        }
        
        } # end year loop
  
      cat(glue('---- STEP_7: combine flux calcs'), file = log_fp, sep = "\n", append = TRUE)
      
      wrtds_out <- flux_annual_wrtds %>%
        filter(wy %in% good_years) %>%
        rename(val = flux) %>%
        mutate(site_code = site_code,
               var = solutes[j],
               method = 'wrtds',
               ms_recommended = 0)
  
      # TODO: LOG every solute completion (and some stats perhaps)
      target_solute_out_frame <- rbind(year_out_frame, wrtds_out) %>%
        arrange(desc(method), desc(wy))
          
      writeLines(glue("---- complete: {target_solute}"))
      
      domain = dmn_nm
      site = site_code
      solute = target_solute
      flux_calc = TRUE
      flux_calc_dt = Sys.time()
      flux_calc_tz = Sys.timezone()
            
      ## TRACKER DF BIND
      if(!exists("main_flux_frame")) {
        main_flux_frame <- data.frame(domain, site, solute, flux_calc, flux_calc_dt, flux_calc_tz)
      } else {
        this_flux_frame <- c(domain, site, solute, flux_calc, flux_calc_dt, flux_calc_tz)
        main_flux_frame <- rbind(main_flux_frame, this_flux_frame)
      }
      
      fp <- glue(file.path('data/ms/main_flux_tracker.csv'))
      write_csv(main_flux_frame, fp)
      
      file_path <- glue('{directory}/{s}.feather',
                        s = site_code)
      
      tryCatch(
        expr = {
          cat(glue('---- STEP_7_SUCCESS: {target_solute} complete'), file = log_fp, sep = "\n", append = TRUE)
          return(target_solute_out_frame)
        },
        error = function(e) {
          cat(glue('---- STEP_7_ERROR: {target_solute} failed'), file = log_fp, sep = "\n", append = TRUE)
        }
      )
    
      
    } # end solute loop
      
    directory <- glue(file.path(data_dir,'stream_flux'))
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    
    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)
    
    
    cat(glue('-- STEP_8: {site_code} flux attempting write'), file = log_fp, sep = "\n", append = TRUE)
    tryCatch(
      expr = {
        cat(glue('-- STEP_8_SUCCESS: {site_code} flux written to file'), file = log_fp, sep = "\n", append = TRUE)
        write_feather(outframe, file_path)
      },
      error = function(e) {
        cat(glue('-- STEP_8_ERROR: {site_code} flux failed to write to file'), file = log_fp, sep = "\n", append = TRUE)
      }
    )
    
    end.tm <- Sys.time()
    end.tz <- Sys.timezone()
    end.dt <- Sys.Date()
    
    cat(glue('STEP_8: {site_code} annual flux\n'), file = log_fp, sep = "\n")
    cat(glue('end, time: {end.tm} {end.tz}\ndate:{end.dt}'), file = log_fp, sep = "\n", append = TRUE)
    
  } # end site loop (for loop)
} # end domain loop

stopCluster(cl)
rstudioapi::restartSession()
