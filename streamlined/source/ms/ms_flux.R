# MacroSheds
library(macrosheds)
source('streamlined/source/usgs_helpers.R')
source('source/helper_functions.R')

# load in a df of MS sites with
ms <- ms_download_site_data() %>%
  select(site_code, ws_area_ha, latitude, longitude) %>%
  rename(lat = latitude,
         long = longitude)

# find all ms sites with nitrate
ms_vars <- ms_catalog()
ms_sites <- ms_vars[grep("^Nitrate", ms_vars$variable_name),]

# calculate the amount of days in ten years
ten_yrs_days <- 365 * 10

ms_sites <- ms_sites %>%
    mutate(period = as.numeric(
        difftime(last_record_utc,
                 first_record_utc,
                 units = "days"
                 ))) %>%
    # to get records > 60 years
  filter(period > ten_yrs_days)

# get Q for RBI calc
ms_sites_ls <- unique(ms_sites$site_code)
raw_data_q <- ms_load_product(
    my_ms_dir,
    prodname = "discharge",
    site_codes = ms_sites_ls,
    sort_result = TRUE,
    warn = TRUE
)

# for each MS sites
s <- ms_sites_ls[1]
ms_dir <- 'streamlined/data/ms/'

for(s in ms_sites_ls) {
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
    filter(ms_status == 0) %>%
    pivot_wider(id_cols = c(datetime, site_code),
                names_from = var,
                values_from = val)

  ms_chem <- raw_data_chem %>%
    filter(ms_status == 0,
           grepl('NO3', var)) %>%
    pivot_wider(id_cols = c(datetime, site_code),
                names_from = var,
                values_from = val)

  ms_df <- merge(ms_q, ms_chem, by = c('datetime', 'site_code')) %>%
    mutate(water_year = water_year(datetime, origin = "usgs"),
           date = date(datetime)) %>%
    group_by(site_code, water_year, date) %>%
    summarise(nitrate_n_mgl = mean(GN_NO3_N),
              q_lps = mean(IS_discharge)) %>%
    select(site_code, water_year, date, nitrate_n_mgl, q_lps)

  # loop thru water years
  ## y <- 1989
  for(y in unique(ms_df$water_year)) {
    ms_df_wy <- ms_df %>%
      filter(water_year == y)

    if(nrow(ms_df_wy[!is.na(ms_df_wy$nitrate_n_mgl),]) < 311) {
      print(paste(s, y, 'insufficient data, skipped'))
      next
    } else {
      ### TRUTH (via composite) ######
        # first make c:q rating
        paired_df <- ms_df_wy %>%
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
        rating_filled_df <- ms_df_wy %>%
            mutate(con_reg = 10^(intercept+(slope*log10(q_lps)))) %>%
            select(date, con_reg) %>%
          full_join(., ms_df_wy, by = 'date') %>%
          select(site_code = site_code.x, water_year = water_year.x, date, con = nitrate_n_mgl, con_reg, q_lps) %>%
          mutate(res = con_reg-con,
                   res = imputeTS::na_interpolation(res),
                   con_com = con_reg-res)

        rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0

        # calculate true annual flux
        true_flux <- rating_filled_df %>%
            mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
            group_by(water_year) %>%
            summarize(flux = sum(flux)) %>%
            mutate(method = 'true',
                   thin = 'none')

        ### THIN data to selected intervals #######
        ms_raw <- merge(ms_q, ms_chem, by = c('datetime', 'site_code'), .keep_all = TRUE) %>%
          mutate(nitrate_n_mgl = GN_NO3_N,
                 q_lps = IS_discharge,
                 wy = water_year(datetime, origin = "usgs")) %>%
          select(site_code, wy, datetime, nitrate_n_mgl, q_lps)

       # MacroSheds is already thinned to daily
        thinned_weekly_c <- ms_raw %>%
          filter(wy == y) %>%
          filter(!is.na(nitrate_n_mgl)) %>%
          filter(lubridate::wday(datetime) == 3) %>%
            mutate(date = lubridate::date(datetime)) %>%
            distinct(date, .keep_all = T) %>%
            select(-datetime) %>%
            rename(con = nitrate_n_mgl)

        thinned_biweekly_c <- ms_raw %>%
          filter(wy == y) %>%
          filter(!is.na(nitrate_n_mgl)) %>%
          filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
          mutate(date = lubridate::date(datetime)) %>%
          distinct(date, .keep_all = T) %>%
          select(-datetime) %>%
          rename(con = nitrate_n_mgl)

        thinned_monthly_c <- ms_raw %>%
          filter(wy == y) %>%
          filter(!is.na(nitrate_n_mgl)) %>%
            mutate(date = date(datetime)) %>%
            filter(day(date) == 1) %>%
            distinct(date, .keep_all = T) %>%
            select(-datetime) %>%
            rename(con = nitrate_n_mgl)

        nmonth <- nrow(thinned_monthly_c)

        thinned_quarterly_c <- rbind(thinned_monthly_c[1,],
                             thinned_monthly_c[as.integer(.5*nmonth),],
                             thinned_monthly_c[as.integer(.75*nmonth),],
                             thinned_monthly_c[nmonth,])

    }

  }

}

# get only if > 85% days coevred
q_check <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        distinct(., date, .keep_all = TRUE) %>%
        mutate(water_year = water_year(datetime, origin = "usgs")) %>%
        group_by(site_code, water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 311)
# all fine!
                                        #

# calculate flux
