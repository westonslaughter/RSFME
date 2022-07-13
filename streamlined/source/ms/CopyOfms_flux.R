# MacroSheds
library(devtools)
#install_github("https://github.com/MacroSHEDS/macrosheds.git")
library(macrosheds)
library(here)
library(RiverLoad)

source('streamlined/source/usgs_helpers.R')
source('source/helper_functions.R')

#declare ms core dir
my_ms_dir <- here('streamlined/data/ms')

# load in a df of MS sites with
ms <- ms_download_site_data() %>%
  select(site_code, ws_area_ha, latitude, longitude) %>%
  rename(lat = latitude,
         long = longitude)


# find all ms sites with nitrate
ms_vars <- ms_catalog()
ms_sites <- ms_vars[grep("^Nitrate", ms_vars$variable_name),]

# calculate the amount of days in ten years
#ten_yrs_days <- 365 * 10

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

#acutal_ms_site_ls <- 
actual_files <- list.files('streamlined/data/ms', full.names = T, recursive = T)
actual_files_ls <- tibble(file = actual_files) %>%
    filter(grepl('stream_chemistry', file),
       !grepl('documentation', file)) %>%
    mutate(site_file = str_split_fixed(file, '/', n = Inf)[,6],
           site_code =  str_split_fixed(site_file, '[.]', n = Inf)[,1]) %>%
    filter(site_code %in% unique(raw_data_q$site_code)) %>%
    pull(site_code)


# for each MS sites
s <- actual_files_ls[1]
ms_dir <- my_ms_dir

for(s in actual_files_ls) {
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
    select(datetime, site_code, nitrate_n_mgl = val) %>%
      na.omit()

  #### isolate to full years######
  
  q_check <- ms_q %>%
      mutate(date = date(datetime)) %>%
      na.omit() %>%
      distinct(., date, .keep_all = TRUE) %>%
      mutate(water_year = water_year(datetime, origin = "usgs")) %>%
      group_by(water_year) %>%
      summarise(n = n()) %>%
      filter(n >= 311)
  
  conc_check <- ms_chem %>%
      na.omit()
  if(nrow(conc_check != 0)){
    conc_check <- conc_check %>%
      mutate(date = date(datetime)) %>%
      distinct(., date, .keep_all = TRUE) %>%
      mutate(water_year = water_year(date, origin = "usgs")) %>%
      group_by(water_year) %>%
      summarise(n = n()) %>%
      filter(n >= 311)
  
  conc_good_years <- conc_check$water_year}
  if(nrow(conc_check) == 0){conc_good_years = NA}
  
  q_good_years <- q_check$water_year
  good_years <- q_good_years[q_good_years %in% conc_good_years]
  
  if(length(good_years) == 0){print(paste0('no good years at ', s))}
  if(length(good_years) != 0){
  
  #### join data and cut to good years ######
  daily_data_n <- ms_chem
  
  daily_data_q <- ms_q %>%
      group_by(date = as_date(datetime)) %>%
      summarize(q_lps = mean(IS_discharge)) %>%
      select(datetime = date, q_lps)
  
  raw_data_full_pre <- full_join(daily_data_n, daily_data_q, by = c('datetime')) %>%
      mutate(wy = water_year(datetime, origin = 'usgs')) %>%
      filter(wy %in% good_years)
  
  
  
  # Loop through good years ####
  # needs DAILY q and any chem
  for(t in 1:length(good_years)){
      target_year <- as.numeric(as.character(good_years[t]))
      raw_data_full <- raw_data_full_pre %>%
          mutate(wy = as.numeric(as.character(wy))) %>%
          filter(wy == target_year)
      ### TRUTH (via composite) ######
      # first make c:q rating
      paired_df <- raw_data_full%>%
          rename(con = nitrate_n_mgl) %>%
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
          select(site_code, datetime, con = nitrate_n_mgl, con_reg, q_lps, wy) %>%
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
      prep_data_q <- ms_q %>%
          mutate(wy = water_year(datetime, origin = 'usgs')) %>%
          filter(wy == target_year) %>%
          select(site_code, date = datetime, q_lps = IS_discharge, wy)
      
      
      ### THIN data to selected intervals #######
      thinned_daily_c <- ms_chem %>%
          #filter(hour(datetime) %in% c(13:18)) %>%
          mutate(date = lubridate::date(datetime),
                 wy = water_year(datetime, origin = 'usgs')) %>%
          distinct(date, .keep_all = T) %>%
          select(-datetime) %>%
          rename(con = nitrate_n_mgl) %>%
          filter(wy == target_year)
      
      thinned_weekly_c <- thinned_daily_c %>%
          #mutate(wy = water_year(datetime, origin = 'usgs')) %>%
          filter(wy %in% good_years) %>%
          #filter(hour(datetime) %in% c(13:18)) %>%
          filter(lubridate::wday(date) == 3) %>%
          #mutate(date = lubridate::date(datetime)) %>%
          distinct(date, .keep_all = T) %>%
          filter(wy == target_year)
      
      thinned_biweekly_c <- thinned_daily_c %>%
          #mutate(wy = water_year(datetime, origin = 'usgs')) %>%
          #filter(wy %in% good_years) %>%
          #filter(hour(datetime) %in% c(13:18)) %>%
          filter(lubridate::mday(date) %in% c(1, 15)) %>%
          #mutate(date = lubridate::date(datetime)) %>%
          distinct(date, .keep_all = T) %>%
          filter(wy == target_year)
      
      thinned_monthly_c <- thinned_daily_c %>%
          #mutate(wy = water_year(datetime, origin = 'usgs')) %>%
          #filter(wy %in% good_years) %>%
          #filter(hour(datetime) %in% c(13:18)) %>%
          #mutate(date = date(datetime)) %>%
          filter(day(date) == 1) %>%
          distinct(date, .keep_all = T) %>%
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
      
      ## save flux out of loop ####
      directory <- glue('streamlined/out_ms/{wy}',
                        wy = target_year)
      if(!dir.exists(directory)){
          dir.create(directory, recursive = TRUE)
      }
      file_path <- glue('{directory}/{s}.feather',
                        s = s)
      
      write_feather(out_frame, file_path)
      # ## take meta and save out of loop ####
      # meta_n <- raw_data_n %>%
      #     filter(wy == target_year)
      # 
      # out_meta <- tibble(max_q_lps = max(prep_data_q$q_lps, na.rm = T),
      #                    min_q_lps = min(prep_data_q$q_lps, na.rm = T),
      #                    mean_q_lps = mean(prep_data_q$q_lps, na.rm = T),
      #                    sd_q_lps = sd(prep_data_q$q_lps, na.rm = T),
      #                    max_n = max(meta_n$val, na.rm = T),
      #                    min_n = min(meta_n$val, na.rm = T),
      #                    mean_n = mean(meta_n$val, na.rm = T),
      #                    sd_n = sd(meta_n$val, na.rm = T))
      # 
      # 
      # directory <- glue('streamlined/data/meta/{wy}',
      #                   wy = target_year)
      # if(!dir.exists(directory)){
      #     dir.create(directory, recursive = TRUE)
      # }
      # file_path <- glue('{directory}/{s}.feather',
      #                   s = site_no)
      # 
      # write_feather(out_meta, file_path)
  }# end year loop
} # end site loop

}
