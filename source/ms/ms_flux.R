# MacroSheds
library(devtools)
#install_github("https://github.com/MacroSHEDS/macrosheds.git")
library(macrosheds)
library(here)
library(RiverLoad)
library(stringr)
library(tidyverse)
library(lfstat)

source('source/usgs_helpers.R')
source('source/helper_functions.R')
source('source/flux_methods.R')
source('source/egret_overwrites.R')

#declare ms core dir
my_ms_dir <- here('data/ms')

# load in a df of MS sites with
ms <- ms_download_site_data() %>%
  select(site_code, ws_area_ha, latitude, longitude) %>%
  rename(lat = latitude,
         long = longitude)


# find all ms sites with nitrate
ms_vars <- ms_catalog()
ms_sites <- ms_vars[grep("^Nitrate", ms_vars$variable_name),]

# calculate the amount of days in ten years
ms_sites <- ms_sites %>%
    mutate(period = as.numeric(
        difftime(last_record_utc,
                 first_record_utc,
                 units = "days"
                 ))) %>%
    filter(period >= 365)
    # to get records > 10 years
  #filter(period > ten_yrs_days)

# get Q for RBI calc
ms_sites_ls <- unique(ms_sites$site_code)
raw_data_q <- ms_load_product(
    my_ms_dir,
    prodname = "discharge",
    site_codes = ms_sites_ls,
    sort_result = TRUE,
    warn = TRUE
)

actual_files <- list.files('data/ms', full.names = T, recursive = T)
actual_files_ls <- tibble(file = actual_files) %>%
    filter(grepl('stream_chemistry', file),
       !grepl('documentation', file)) %>%
  mutate(site_file = str_split(file, '/', simplify = TRUE)[,5],
         site_code = str_extract(site_file, "[^\\.]*")) %>%
    filter(site_code %in% unique(raw_data_q$site_code)) %>%
    pull(site_code)


# for each MS sites
s <- actual_files_ls[1]
ms_dir <- my_ms_dir

for(s in actual_files_ls) {
  ## s <- actual_files_ls[1]
  # load in site area
  area <- ms[ms$site_code == s,]$ws_area_ha

  # load in site Q
  raw_data_q <- ms_load_product(
      ms_dir,
      prodname = "discharge",
      site_codes = s,
      sort_result = TRUE,
      warn = TRUE
  )

  # load in stream chemistry
  raw_data_chem <- ms_load_product(
      ms_dir,
      prodname = "stream_chemistry",
      site_codes = s,
      sort_result = TRUE,
      warn = TRUE
  )

  ms_q <- raw_data_q %>%
    #filter(ms_status == 0) %>%
    pivot_wider(id_cols = c(datetime, site_code),
                names_from = var,
                values_from = val)

  ms_chem <- raw_data_chem %>%
    filter(ms_status == 0,
           grepl('NO3', var)) %>%
    select(datetime, site_code, nitrate_n_mgl = val)


  ms_df <- full_join(ms_q, ms_chem, by = c('datetime', 'site_code')) %>%
    mutate(water_year = water_year(datetime, origin = "usgs"),
           date = date(datetime)) %>%
    group_by(site_code, water_year, date) %>%
    summarise(nitrate_n_mgl = mean(nitrate_n_mgl, na.rm = T),
              q_lps = mean(IS_discharge, na.rm = T)) %>%
    select(site_code, water_year, date, nitrate_n_mgl, q_lps)

  # loop thru water years
  ## y <- 1989
  for(y in unique(ms_df$water_year)) {
    ms_df_wy <- ms_df %>%
        filter(water_year == y) %>%
        rename(wy = water_year)
    
    chem_df <- ms_df_wy %>%
        ungroup()%>%
        select(date, con = nitrate_n_mgl, site_code, wy) %>%
        na.omit()
    
    q_df <- ms_df_wy %>%
        ungroup() %>%
        select(date, q_lps, site_code, wy)
    
    chem_check <- nrow(ms_df_wy[!is.na(ms_df_wy$nitrate_n_mgl),]) > 311
    q_check <- nrow(ms_df_wy[!is.na(ms_df_wy$q_lps),]) > 311
    com_check <- chem_check*q_check

    if(com_check == 0 ) {
      print(paste(s, y, 'insufficient data, skipped'))
      next
    } else {
    
      pw_wy <- calculate_pw(chem_df = chem_df, q_df = q_df)
    
      beale_wy <- calculate_beale(chem_df = chem_df, q_df = q_df)
    
      rating_wy <- calculate_rating(chem_df = chem_df, q_df = q_df)
    
      composite_wy <- calculate_composite_from_rating_filled_df(
        generate_residual_corrected_con(chem_df = chem_df, q_df = q_df)
    )
    
    out_frame <- tibble(wy = y,
                        flux = c(pw_wy[[1]], beale_wy[[1]], 
                                 rating_wy [[1]], composite_wy$flux[[1]]), 
                        site_code = s,
                        method = c('pw', 'beale', 'rating', 'composite'))
    
    ## save flux out of loop ####
    directory <- glue('streamlined/out_ms/{wy}',
                      wy = y)
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    file_path <- glue('{directory}/{s}.feather',
                      s = s)
    
    write_feather(out_frame, file_path)

  }

  }
}


