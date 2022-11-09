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

data_dir <- here('data/ms/hbef/')
site_files  <- list.files('data/ms/hbef/discharge', recursive = F)
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
        ## switch to only Nitrate for now
        #filter(var == 'GN_NO3_N') %>%
        ## filter(!str_detect(var, 'Charge'),
        ##        !str_detect(var, 'temp'),
        ##        !str_detect(var, 'pH')) %>%
        select(var) %>%
        unique() %>%
        pull(var)

  writeLines(paste("FLUX CALCS:", site_code))

  # Loop through solutes at site #####
  ## j = 1
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

    ## find acceptable years ####
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

    # Check WRTDS assumptions #####
    if(n_yrs < 5){
        next
    }else{

    # get data ready for egret adaptation ####
    # chem_df <- raw_data_con %>%
    #     mutate(date = date(datetime)) %>%
    #     group_by(date) %>%
    #     summarize(val = mean(val)) %>%
    #     mutate(site_code = !!site_code, var = 'con',
    #            wy = water_year(date, origin = 'usgs')) %>%
    #     select(site_code, datetime = date, con = val, wy) %>%
    #     na.omit()

    q_df <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean(val)) %>%
        mutate(site_code = !!site_code, var = 'q_lps',
               wy = water_year(date, origin = 'usgs')) %>%
        #select(site_code, datetime = date, var, val, wy) %>%
        select(site_code, datetime = date, q_lps = val, wy) %>%
        na.omit()


    con_full <- raw_data_con %>%
        mutate(wy = as.numeric(as.character(water_year(datetime, origin = 'usgs'))),
               site_code = site_code) %>%
        select(site_code, datetime, con = val, wy) %>%
        ## filter(wy < 1975) %>%
        na.omit()

    #### calculate WRTDS ######
    flux <- calculate_wrtds(
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


        } # end good data else statement
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
