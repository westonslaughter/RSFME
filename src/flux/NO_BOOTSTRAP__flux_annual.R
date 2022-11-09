library(tidyverse)
library(feather)
library(here)
library(glue)
library(lubridate)
library(EGRET)
library(macrosheds)

source('source/helper_functions.R')
source('source/egret_overwrites.R')
source('ms_overwrites.R')
source('source/flux_methods.R')
source('source/usgs_helpers.R')

data_dir <- here('data/ms/')
## site_files  <- list.files(data_dir, recursive = F)

hbef_files  <- list.files('data/ms/hbef/discharge', recursive = F)
hj_files  <- list.files('data/ms/hjandrews/discharge', recursive = F)
site_files <- c(hbef_files, hj_files)

site_info  <- read_csv(here('data/site/ms_site_info.csv'))

# df to populate with annual flux values by method
out_frame <- tibble(wy = as.integer(),
                    site_code = as.character(),
                    var = as.character(),
                    val = as.numeric(),
                    method = as.character(),
                    ms_reccomended = as.integer(),
                    ms_interp_ratio = as.numeric(),
                    ms_status_ratio = as.numeric(),
                    ms_missing_ratio = as.numeric())
## i = 1
# Loop through sites #####
## hbef_diTopic: Weston Slaughter's Personal Meeting Room

hbef_dir <- here('data/ms/hbef/')
hj_dir <- here('data/ms/hjandrews/')

choices <- c('hbef')

for(choice in choices) {
  if(choice == 'hbef') {
    data_dir <- hbef_dir
    site_files <- hbef_files
  } else if(choice == 'hjandrews') {
    data_dir <- hj_dir
    site_files <- hj_files
  }
}

# list all networks
## networks <- list.files(data_dir, recursive = F)
## networks <- networks[!networks %in% c("hbef", "hjandrews", "arctic")]
## networks <- "hbef"

# read in variables data
ms_flux_vars <- ms_download_variables() %>%
  filter(flux_convertible == 1) %>%
  pull(variable_code)

for(nwk in networks){

  writeLines(paste("\n\n starting annual flux calculations for netowrk:", nwk, "\n\n"))

  data_dir <- here(glue('data/ms/{network}', network = nwk))
  site_files  <- list.files(glue('data/ms/{network}/discharge', network = nwk), recursive = F)

  for(i in 1:length(site_files)){

    site_file <- site_files[i]
    site_code <- strsplit(site_file, split = '.feather')[[1]]

    area <- site_info %>%
        filter(site_code == !!site_code) %>%
        distinct() %>%
        pull(ws_area_ha)

    lat <- site_info %>%
        filter(site_code == !!site_code) %>%
        distinct() %>%
        pull(Y)

    long <- site_info %>%
        filter(site_code == !!site_code) %>%
        distinct() %>%
        pull(X)

    # read in chemistry data
    raw_data_con_in <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
    # NOTE: is there only temp data at GS Watershed 3?
      filter(ms_interp == 0,
             # dropping all vars that are not marked as 'flux convertible'
             ms_drop_var_prefix(var) %in% ms_flux_vars)

    # read in discharge data
    raw_data_q <- read_feather(here(glue(data_dir, '/discharge/', site_code, '.feather')))

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

  if(length(solutes) == 0) {
    writeLines(
      glue("\n\n no flux convertible solutes in stream chemistry data for {site}", site = site_code),
      "\n   skipping to next site\n"
      )
    next
  }

  writeLines(paste("FLUX CALCS:", site_code))

  ## Loop through solutes at site #####
  ## j = 1
  ## j = 5 # Ca
  ## j = 27 # SpCond
  for(j in 1:length(solutes)){
    writeLines(paste("site:", site_code,
                     "var:", solutes[j]))

    #set to target solute
    target_solute <- solutes[j]

    raw_data_con <- read_feather(here(glue(data_dir, '/stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0,
               val > 0) %>%
        filter(var == target_solute) %>%
        select(datetime, val) %>%
        na.omit()

    # if site has no data for var, skip to next var
    if(nrow(raw_data_con) == 0) next

    # find acceptable years
    q_check <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        # NOTE: should we filter out NAs?
        filter(ms_interp == 0, !is.na(val)) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(datetime, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)

    conc_check <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        # NOTE: should we filter out NAs?
        filter(!is.na(val)) %>%
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

    # NOTE: adding handling if concentration data fails conc check
    if(nrow(conc_check) < 1) {
      writeLines(glue("{site} concentration data insufficient sample size and frequency to warrant flux estimation",
                      "\n   no water years in {site} dataset with minimum standards met", site = site_code))
      next
    } else if(nrow(q_check) < 1) {
      writeLines(glue("{site} discharge data insufficient sample size and frequency to warrant flux estimation",
                      "\n   no water years in {site} dataset with minimum standards met", site = site_code))
      next
    } else if(length(good_years) == 0) {
      writeLines(glue("no water years where q data and concentration data both meet minimum standards",
                      "skipping site: {site}", site = site_code))
      next
    }

    #join data and cut to good years
    daily_data_con <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean_or_x(val)) %>%
        mutate(site_code = !!site_code, var = 'con') %>%
        select(site_code, datetime = date, var, val)

    daily_data_q <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean_or_x(val)) %>%
        mutate(site_code = !!site_code, var = 'q_lps') %>%
        select(site_code, datetime = date, var, val)

    q_df <- daily_data_q %>%
      pivot_wider(names_from = var,
                  values_from = val)

    raw_data_full <- rbind(daily_data_con, daily_data_q) %>%
        pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
        mutate(wy = water_year(datetime, origin = 'usgs')) %>%
        filter(wy %in% good_years)

    con_full <- raw_data_full %>%
          mutate(wy = as.numeric(as.character(wy))) %>%
            select(site_code, datetime, con, wy) %>%
            ## filter(wy < 1975) %>%
            na.omit()

    #### calculate WRTDS ######
    tryCatch(
      expr = {
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
        writeLines(paste('\nWRTDS run failed for \n     site', site_code,
                         '\n     variable', target_solute, '\n WRTDS TRYING AGAIN'))
        tryCatch(
          expr = {
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
            print("WRTDS failed, setting to NA")
            flux_annual_wrtds <- NA
          }
        )

        ## next
      }
    )

    ## write_feather(raw_data_full, "data/ms/hbef/true/w3_chem_samples.feather")
    ## k = 1
    ### Loop through good years #####
    ## for(k in 42:47) {
    for(k in 1:length(good_years)){

      writeLines(paste("site:", site_code,
                       'year:', good_years[k]))

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

        ### calculate annual flux ######
        chem_df_errors <- con_target_year
        q_df_errors <- q_target_year

        ### save and then remove errors attribute for calcs
        chem_df <- errors::drop_errors(chem_df_errors)
        q_df <- errors::drop_errors(q_df_errors)

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

        # ``model_data`` is the site-variable-year dataframe of Q and concentration
        # log-log rating curve of C by Q
        rating <- summary(lm(model_data$c_log ~ model_data$q_log, singular.ok = TRUE))
        # R^2 value of rating curve
        r_squared <- rating$r.squared
        # auto-correlation of residuals of rating curve
        resid_acf <- abs(acf(rating$residuals, lag.max = 1, plot = FALSE)$acf[2])
        # auto-correlation of concentration data
        con_acf <- abs(acf(paired_df$con, lag.max = 1, plot = FALSE)$acf[2])

        # modified from figure 10 of Aulenbach et al 2016
      if(!is.nan(r_squared)) {
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
      } else {
        writeLines("\n\n ideal method error: r_squared value was NaN, ideal method set to NA\n\n")
        ideal_method <- NA
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
        out_frame <- rbind(out_frame, target_year_out)

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

    directory <- glue(data_dir,'stream_flux/')
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }

    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)

  ## out_frame_save <- out_frame
  out_frame <- errors::drop_errors(out_frame)

  write_feather(out_frame, file_path)
} # end site loop

## w3 <- read_feather('data/ms/hbef/stream_flux/w3.feather')
  ## filter(var == 'GN_DOC') %>%
  ## pivot_wider(id_cols = c('site_code', 'wy', 'var'),
  ##             values_from = val,
  ##             names_from = c('method'))
}
