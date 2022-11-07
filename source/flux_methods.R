# Flux Functions
# PREP and HELPERS
### Riverload conversion function #####
prep_raw_for_riverload <- function(chem_df, q_df, datecol = 'date'){
    conv_q <- q_df %>%
                mutate(datetime = as.POSIXct(get(datecol), format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
                mutate(flow = q_lps*0.001) %>% # convert lps to cubic meters per second)
                select(datetime, flow) %>%
                arrange(datetime) %>%
                data.frame()
    if(month(conv_q$datetime[1]) == 9 & day(conv_q$datetime[1]) == 30){
        conv_q = conv_q[-1,]
    }

    conv_c <- chem_df %>%
                mutate(datetime = as.POSIXct(get(datecol), format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
                select(datetime, con) %>%
                data.frame()

    db <- full_join(conv_q, conv_c, by = "datetime") %>%
                #filter(!is.na(flow)) %>%
                arrange(datetime)

    return(db)
}


# FLUX CALCS
###### calculate period weighted#########
calculate_pw <- function(chem_df, q_df, datecol = 'date', period = NULL){
  rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)

  if(is.na(rl_data[1,2])){
      rl_data <- rl_data[-1,]
  }

  if(is.null(period)){
  flux_from_pw <- method6(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
  }else{

  if(period == 'month'){

      method6_month <- function (db, ncomp, period){
          if (requireNamespace("imputeTS")) {
              n <- nrow(db)
              interpolation <- data.frame(imputeTS::na_interpolation(db,
                                                                     "linear"))
              load <- data.frame(interpolation[, 2] * interpolation[,
                                                                    -c(1:2)])
              difference <- matrix(nrow = (nrow(db) - 1), ncol = 1)
                for (i in 1:(nrow(db) - 1)) {
                    difference[i] <- difftime(db[i + 1, 1], db[i, 1],
                                            units = "days")
            }

              loadtot <- cbind.data.frame(interpolation$datetime[-nrow(interpolation)],
                                          flux)
              colnames(loadtot)[1] <- c("datetime")
              loadtot[, 1] <- format(as.POSIXct(loadtot[, 1]), format = "%Y-%m")
              forrow <- aggregate(loadtot[, 2] ~ datetime, loadtot,
                                  sum)
              agg.dataC <- matrix(nrow = nrow(forrow), ncol = (ncomp))
              for (i in 1:ncomp) {
                  agg.data <- aggregate(loadtot[, i + 1] ~ datetime,
                                        loadtot, sum)
                  agg.dataC[, i] <- as.matrix(agg.data[, 2])
              }
              colnames(agg.dataC) <- c(names(db)[3:(ncomp + 2)])
              rownames(agg.dataC) <- forrow$datetime
              return(agg.dataC)
          }

              colnames(agg.dataC) <- c(names(db)[3:(ncomp + 2)])
              rownames(agg.dataC) <- forrow$datetime
              return(agg.dataC)
      }

      ##### apply #####

      flux_from_pw <- method6_month(rl_data, ncomp = 1, period = period)

      flux_from_pw <- tibble(date = rownames(flux_from_pw),
                          flux = (flux_from_pw[,1]/(1000*area)))
}
}
  return(flux_from_pw)
}

###### calculate beale ######
calculate_beale <- function(chem_df, q_df, datecol = 'date', period = NULL){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)

    if(is.na(rl_data[1,2])){
        rl_data <- rl_data[-1,]
    }

    if(is.null(period)){
    flux_from_beale <- beale.ratio(rl_data, ncomp = 1) %>%
      sum(.)/(1000*area)
    }else{

    if(period == 'month'){
        flux_from_beale <- beale.ratio(rl_data, ncomp = 1, period = period)

        flux_from_beale <- tibble(date = rownames(flux_from_beale),
                               flux = (flux_from_beale[,1]/(1000*area)))
    }
    }

    return(flux_from_beale)
}

##### calculate rating #####
calculate_rating <- function(chem_df, q_df, datecol = 'date', period = NULL){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)

    if(is.null(period)){
    flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1) %>%
        sum(.)/(1000*area)
    }else{

    if(period == 'month'){
        flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1, period = period)

        flux_from_reg <- tibble(date = rownames(flux_from_reg),
                                  flux = (flux_from_reg[,1]/(1000*area)))
    }
    }

    return(flux_from_reg)
}


##### calculate WRTDS #####
calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long, datecol = 'date', agg = 'default', minNumObs = 2, minNumUncen =2, gap_period = 730) {
  tryCatch(
    expr = {
      # default sums all daily flux values in df
      egret_results <- adapt_ms_egret(chem_df, q_df, ws_size,
                                      lat, long,
                                      datecol = datecol,
                                      minNumObs = minNumObs,
                                      minNumUncen = minNumUncen,
                                      gap_period = gap_period)

      if(agg == 'default') {
        flux_from_egret <- egret_results$Daily$FluxDay %>%
          warn_sum(.)/(area)
      } else if(agg == 'annual') {
        flux_from_egret <- egret_results$Daily %>%
                    mutate(
                      wy = water_year(Date)
                    ) %>%
          group_by(wy) %>%
          summarize(
            flux = warn_sum(FluxDay)/(area)
          )
      } else if(agg == 'monthly') {
        flux_from_egret <- egret_results$Daily %>%
                    mutate(
                      wy = water_year(Date),
                      month = lubridate::month(Date)
                    ) %>%
          group_by(wy, month) %>%
          summarize(
            flux = warn_sum(FluxDay)/(area)
          )
      }
        },
    error = function(e) {
            print('ERROR: WRTDS failed to run')
            return(NA)
        })
    return(flux_from_egret)
}

generate_residual_corrected_con <- function(chem_df, q_df, datecol = 'date', sitecol = 'site_no'){
        # first make c:q rating
        paired_df <- q_df %>%
            full_join(chem_df, by = c(datecol, sitecol, 'wy')) %>%
            na.omit() %>%
            filter(q_lps > 0,
                   is.finite(q_lps))

        if(nrow(paired_df) <= 2){return(NA)}else{

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
        site_code <- paired_df$site_code[[1]]

        rating_filled_df <- q_df %>%
          mutate(con_reg = 10^(intercept+(slope*log10(q_lps)))) %>%
          select(all_of(datecol), con_reg, q_lps) %>%
          full_join(., chem_df, by = datecol) %>%
          select(site_code, all_of(datecol), con, con_reg, q_lps, wy)  %>%
            mutate(res = con_reg-con,
                   res = imputeTS::na_interpolation(res),
                   con_com = con_reg-res,
                   site_code = !!site_code,
                   wy = water_year(get(datecol), origin = 'usgs'))

        rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0
        return(rating_filled_df)
        }
        }

##### calculate monthly flux from composite ####
calculate_composite_from_rating_filled_df <- function(rating_filled_df, site_no = 'site_no', period = NULL){

        if(is.null(period)){
        flux_from_comp <- rating_filled_df %>%
            select(datetime, con_com, q_lps, wy) %>%
            na.omit() %>%
            mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
            group_by(wy) %>%
          summarize(flux = sum(flux)) %>%
            mutate(site_code = site_no)
        }else{

        if(period == 'month'){
            flux_from_comp <- rating_filled_df %>%
                mutate(month = month(datetime),
                       flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
                group_by(wy, month) %>%
                summarize(date = max(datetime),
                          flux = sum(flux))
        }
        }

        return(flux_from_comp)
        }

# functions for adapt_ms_egret
detect_record_break <- function(data) {
    data_time <- data %>%
      select(Date, Julian) %>%
      # how many days between this day and the next
      mutate(days_gap = lead(Julian, 1) - Julian)
    return(data_time)
  }

get_break_dates <- function(data, gap_period = 730) {
  start = list()
  end = list()

  index = 1
  # default 730 bc USGS says breaks below 2 years probably fine
  for(i in 1:nrow(data)) {
    gap <- data$days_gap[i]

    if(!is.na(gap)) {

    if(gap > gap_period) {
      start[[index]] = c(data$Date[i])
      end[[index]] = c(data$Date[i+1])
      index = index + 1
    }
    }
  }

    break_dates = list('start' = start, 'end'= end)
    return(break_dates)
  }

# wrapper for WRTDS in EGRET
adapt_ms_egret <- function(chem_df, q_df, ws_size, lat, long,
                           site_data = NULL, kalman = FALSE,
                           datecol = 'date', minNumObs = 2, minNumUncen = 2, gap_period = 730){
  # TODO:  reorder site data args to make fully optional

    get_MonthSeq <- function(dates){
        ## dates <- Sample_file$Date
        years <- lubridate::year(dates)-1850

        MonthSeq <- years*12

        MonthSeq <- MonthSeq + lubridate::month(dates)

        return(MonthSeq)
    }


    decimalDate <- function(rawData){

        dateTime <- as.POSIXlt(rawData)
        year <- dateTime$year + 1900

        startYear <- as.POSIXct(paste0(year,"-01-01 00:00"))
        endYear <- as.POSIXct(paste0(year+1,"-01-01 00:00"))

        DecYear <- year + as.numeric(difftime(dateTime, startYear, units = "secs"))/as.numeric(difftime(endYear, startYear, units = "secs"))
        return(DecYear)
    }

    wyday <- function(dates, wy_type = 'usgs') {
        # get earliest and latest date from series
        start_date <- min(dates)
        end_date <- max(dates)

        # get the normal year, and make the start of WY date from that year
        year_firsthalf <- year(start_date)
        wy_start <- paste0(year_firsthalf, '-10-01')
        wy_end <- paste0(year_firsthalf+1, '-09-30')

        # gap between start WY and end of normal year, always 91 days
        correction_val <- yday('2021-12-31') - yday('2021-10-01')


        # adjust yday from start of gregorian year to (usgs) WY year
        wy_firsthalf <- yday(dates[month(dates) %in% c(10, 11, 12)]) - yday(wy_start) + 1
        wy_secondhalf <- yday(dates[!month(dates) %in% c(10, 11, 12)]) + correction_val + 1

        wydays <- c(wy_firsthalf, wy_secondhalf)

      return(wydays)
    }

    decimalDateWY <- function(dates, wy_type = 'usgs') {
        # get earliest and latest date from series
        start_date <- min(dates)
        end_date <- max(dates)

        # get the normal year, and make the start of WY date from that year
        year_firsthalf <- year(start_date)
        wy_start <- paste0(year_firsthalf, '-10-01')
        wy_end <- paste0(year_firsthalf+1, '-09-30')
        wy_seq <- seq(as.Date(wy_start), as.Date(wy_end), by = 1)

        wydays <- wyday(dates) - 1
        wylength <- length(wy_seq)

        ydec <- as.integer(as.character(water_year(dates, wy_type))) + wydays/wylength

      return(ydec)
    }

    # TODO: rename with real descriptive name
    enlightened_yday <- function(dates, wy_type = 'usgs') {
        # get earliest and latest date from series
        start_date <- min(dates)
        end_date <- max(dates)

        # get the normal year, and make the start of WY date from that year
        wy_first <- year(start_date)
        wy_start <- paste0(wy_first, '-10-01')
        wy_end <- paste0(wy_first+1, '-09-30')

        # unenlightened yday of date vector
        wy_firsthalf <- yday(dates[month(dates) %in% c(10, 11, 12)])
        wy_secondhalf <- yday(dates[!month(dates) %in% c(10, 11, 12)])

        # if first half of WY is a leap year,
        if(366 %in% wy_firsthalf) {
          ## print('first half of WY is leap year, adjusting')
          wy_firsthalf = wy_firsthalf - 1
        } else if(366 %in% wy_secondhalf) {
          ## print('second half of WY is leap year, no action')
        } else {
          ## print('neither side of water year is a leap year')
          invisible()
        }

        # adjust yday from start of gregorian year to (usgs) WY year
        wy_ydays <- c(wy_firsthalf, wy_secondhalf)

      return(wy_ydays)
    }

    get_start_end <- function(d){
        start_date <- min(d$datetime)
        start_year <- lubridate::year(start_date)
        start_wy <- ifelse(lubridate::month(start_date) %in% c(10, 11, 12), start_year+1, start_year)
        filter_start_date <- lubridate::ymd(paste0(start_wy-1, '-10-01'))

        end_date <- max(d$datetime)
        end_year <- lubridate::year(end_date)
        end_wy <- ifelse(lubridate::month(end_date) %in% c(10, 11, 12), end_year+1, end_year)
        filter_end_date <- lubridate::ymd(paste0(end_wy, '-10-01'))

        fin_dates <- c(filter_start_date, filter_end_date)
        return(fin_dates)

    }


    ms_run_egret_adapt <- function(stream_chemistry, discharge, prep_data = TRUE,
             run_egret = TRUE, kalman = FALSE, quiet = FALSE,
             site_data = NULL, min_q_method = 'USGS', minNumObs = 2, minNumUncen = 2, gap_period = 730){

        # Checks
        if(any(! c('site_code', 'var', 'val', 'datetime') %in% names(stream_chemistry))){
            stop('stream_chemistry must be a data.frame in MacroSheds format with the columns site_code,
             datetime, var, and, val')
        }

        if(any(! c('site_code', 'var', 'val', 'datetime') %in% names(discharge))){
            stop('discharge must be a data.frame in MacroSheds format with the columns site_code,
             datetime, var, and, val')
        }

        if(! length(unique(macrosheds::ms_drop_var_prefix(stream_chemistry$var))) == 1){
            stop('Only one chemistry variable can be run at a time.')
        }

        if((! length(unique(stream_chemistry$site_code)) == 1) || (! length(unique(discharge$site_code)) == 1)){
            stop('Only one site can be run in EGRET at a time')
        }

        if(! unique(stream_chemistry$site_code) == unique(discharge$site_code)){
            stop('stream_chemistry and discharge must contain the same site_code')
        }

        # Get var and site info
        if(is.null(site_data)){
            site_data <- macrosheds::ms_download_site_data()

            if(! unique(stream_chemistry$site_code) %in% site_data$site_code){
                stop('This site is not in the MacroSheds dataset, provide a site_data table with the names: site_code, ws_area_ha, latitude, longitude')
            }
        } else{
            if(!all(names(site_data) %in% c('site_code', 'ws_area_ha', 'latitude', 'longitude', 'site_type'))){
                stop('If you are not using a macrosheds site, you must supply site_data with a tibble with the names: site_code, ws_area_ha, latitude, longitude')
            }
            site_data <- site_data %>%
                mutate(site_type = 'stream_gauge')
        }

        ## ms_vars <- read.csv('data/ms/macrosheds_vardata.csv')
        ms_vars <- macrosheds::ms_download_variables()
        site_code <- unique(stream_chemistry$site_code)


        #### Prep Files ####

        if(prep_data){
            # EGRET does not like NAs
            stream_chemistry <- stream_chemistry %>%
                filter(!is.na(val))
            discharge <- discharge %>%
                filter(!is.na(val))

            # Egret can't accept 0s in the column for min val (either hack egret or do this or
            # look for detection limits)
            # TODO: use DL system to replace min vals, get_hdl()
            # TODO: AND/OR use EGRET built in uncertainty system, ConcLow = NA, ConcHigh = DL
            ## min_chem <- stream_chemistry %>%
            ##   filter(val > 0,
            ##          !is.na(val),
            ##          !is.infinite(val),
            ##          !is.null(val)) %>%
            ##   pull(val) %>%
            ##     # NOTE: upsettingly, use of errors package makes min() now return NA instead of Inf
            ##   min()

            min_chem <- min(as.numeric(stream_chemistry[!is.infinite(stream_chemistry$val) & !is.na(stream_chemistry$val),]$val))

            # NOTE: changed so now we are forced to have a real number minimum value
            # NOTE: change to no conditional- min chem should never be infinite
            ## if(!is.infinite(min_chem)){
            ## }
            stream_chemistry <- stream_chemistry %>%
                mutate(val = ifelse(val == 0, !!min_chem, val))

            # Filter so there is only Q going into the model that also has chem
            stream_chemistry <- stream_chemistry %>%
                mutate(year = lubridate::year(datetime),
                       month = lubridate::month(datetime)) %>%
                mutate(waterYear = ifelse(month %in% c(10, 11, 12), year+1, year))

            # TODO: filter isn't even active... and why 6?
            # Get years with at least 6 chemistry samples (bi-monthly sampling is a
            # reasonable requirement?)
            years_with_data <- stream_chemistry %>%
                group_by(waterYear) %>%
                summarise(n = n()) %>%
                # filter(n >= 6) %>%
                pull(waterYear)

            # Filter Q and chem to overlap in a water-year
            chem_dates <- get_start_end(stream_chemistry)
            q_dates <- get_start_end(discharge)

            start_date <- if(chem_dates[1] > q_dates[1]) { chem_dates[1] } else{ q_dates[1] }
            end_date <- if(chem_dates[2] < q_dates[2]) { chem_dates[2] } else{ q_dates[2] }

            discharge <- discharge %>%
                filter(datetime >= !!start_date,
                       datetime <= !!end_date)
            stream_chemistry <- stream_chemistry %>%
                filter(datetime >= !!start_date,
                       datetime <= !!end_date)

            # Filter discharge to only include water years with chemistry sampling
            discharge <- discharge %>%
                mutate(year = lubridate::year(datetime),
                       month = lubridate::month(datetime)) %>%
                mutate(waterYear = ifelse(month %in% c(10, 11, 12), year+1, year)) %>%
                filter(waterYear %in% !!years_with_data)

            # Remove times when there is a chem sample and no Q reported
            samples_to_remove <- left_join(stream_chemistry, discharge, by = 'datetime') %>%
                filter(is.na(val.y)) %>%
                pull(datetime)

            stream_chemistry <- stream_chemistry %>%
                filter(!datetime %in% samples_to_remove)
        }

      # last redundant assurance of no NA, Inf, or 0 values
      stream_chemistry <- stream_chemistry %>%
        filter(!is.na(val),
               !is.infinite(val),
               !is.null(val),
               val > 0)

      n_records <- length(stream_chemistry$val)
      n_years <- length(unique(stream_chemistry$wy))

      # if n_records < 100, change minNumObs accordingly
      if(n_records < 100) {
        minNumObs = n_records - ceiling(n_records*0.1)
        minNumUncen = ceiling(minNumObs/2)
        warning(paste('number of samples less than 100, modifying WRTDS arguments',
                      '\n     minNumObs:', minNumObs,
                      '\n     minNumUncen', minNumUncen))
      }


      writeLines(paste('stream chemistry dataframe being passed into EGRET sample, \nfile has:',
                       n_records, 'samples over ', n_years, 'water years'))

        # Set up EGRET Sample file
        Sample_file <- tibble(Name = site_code,
                              Date = as.Date(stream_chemistry$datetime),
                              ConcLow = stream_chemistry$val,
                              ConcHigh = stream_chemistry$val,
                              Uncen = 1,
                              ConcAve = stream_chemistry$val,
                              Julian = as.numeric(julian(lubridate::ymd(stream_chemistry$datetime),origin=as.Date("1850-01-01"))),
                              Month = lubridate::month(stream_chemistry$datetime),
                              Day = enlightened_yday(stream_chemistry$datetime),
                              DecYear = decimalDate(stream_chemistry$datetime),
                              MonthSeq = get_MonthSeq(stream_chemistry$datetime)) %>%
            mutate(SinDY = sin(2*pi*DecYear),
                   CosDY = cos(2*pi*DecYear))  %>%
            mutate(waterYear = ifelse(Month %in% c(10, 11, 12), lubridate::year(Date) + 1, lubridate::year(Date))) %>%
            select(Name, Date, ConcLow, ConcHigh, Uncen, ConcAve, Julian, Month, Day,
                   DecYear, MonthSeq, waterYear, SinDY, CosDY)

        # Set up EGRET Daily file
        # make sure dischagre is daily aggregated
        discharge_daily <- discharge %>%
          mutate(datetime = date(datetime)) %>%
          group_by(datetime, site_code, var, ms_status, ms_interp, year, month) %>%
          summarise(val = mean(val))

        Daily_file <- tibble(Name = site_code,
                             Date = as.Date(discharge_daily$datetime),
                             # converting lps to m^3/s
                             Q = discharge_daily$val/1000,
                             Julian = as.numeric(julian(lubridate::ymd(discharge_daily$datetime),origin=as.Date("1850-01-01"))),
                             Month = lubridate::month(discharge_daily$datetime),
                             Day = enlightened_yday(discharge_daily$datetime),
                             DecYear = decimalDate(discharge_daily$datetime),
                             MonthSeq = get_MonthSeq(discharge_daily$datetime),
                             Qualifier = discharge_daily$ms_status)

        if(prep_data){
            # Egret can't handle 0 in Q, setting 0 to the minimum Q ever reported seem reasonable
            no_flow_days <- Daily_file %>%
                filter(Q == 0) %>%
                pull(Date)

            n_nfd <- length(no_flow_days)
            n_record <- length(Daily_file$Date)
            percent_no_flow <- n_nfd/n_record

            writeLines(paste('days with no flow, percent of record:', percent_no_flow))


          # trying USGS WRTDS flow min method,
          # also NOTE: USGS WRTDS manual says no flow > %0.2 of days is an issue
          if(min_q_method == 'USGS'){
            print('using USGS reccomended no flow replacement method')
            mean_flow <- mean(Daily_file$Q[Daily_file$Q > 0], na.rm = TRUE)

            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!mean_flow, Q))
          } else {
            # NOTE: could this be where Inf shows up too? like in min_chem?
            min_flow <- min(Daily_file$Q[Daily_file$Q > 0], na.rm = TRUE)

            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!min_flow, Q))
          }

          # TODO: record zero flow days, and set flux for those days to zero
          # TODO: USGS WRTDS manual says method for replacing 0 and negative flow
          # is to set all 0 and neg to 0, then replace with 0.1% of mean flow
          # and they say final results should have this small flow increment
          # subtracted from Q and flux results (pg 6, manual)
        }

        Daily_file <- Daily_file %>%
            # 'extend' rollmean NOTE: cpuld rollmean args cause mischief
            # tho right align i believe is correct based off of egret docs``
            mutate(Q7 = zoo::rollmean(Q, 7, fill = 'extend', align = 'right'),
                   Q30 = zoo::rollmean(Q, 30, fill = 'extend', align = 'right'),
                   LogQ = log(Q))

        Daily_file <- tibble::rowid_to_column(Daily_file, 'i') %>%
            select(Name, Date, Q, Julian, Month, Day, DecYear, MonthSeq, Qualifier,
                   i, LogQ, Q7, Q30)

        # Set up INFO table
        var <- macrosheds::ms_drop_var_prefix(unique(stream_chemistry$var))
        var_unit <- ms_vars %>%
            filter(variable_code == !!var) %>%
            pull(unit)
        site_lat <- site_data %>%
            filter(site_code == !!site_code) %>%
            pull('latitude')
        site_lon <- site_data %>%
            filter(site_code == !!site_code) %>%
            pull('longitude')
        site_ws_area <- site_data %>%
            filter(site_code == !!site_code,
                   site_type == 'stream_gauge') %>%
            pull('ws_area_ha')

        # ws area change?
        site_ws_area <- site_ws_area / 100

        new_point <- sf::st_sfc(sf::st_point(c(site_lon, site_lat)), crs = 4326) %>%
            sf::st_transform(., crs = 4267)

        # NOTE: question- should we allow interpolating back to Oct 1st from a
        # record where first ever sample is Oct 26th? (as example of issue)
        INFO_file <- tibble(agency_cd = 'macrosheds',
                            site_no = site_code,
                            station_nm = site_code,
                            site_tp_code = 'ST',
                            lat_va = site_lat,
                            long_va = site_lon,
                            dec_lat_va = site_lat,
                            dec_long_va = site_lon,
                            coord_meth_cd = NA,
                            coord_acy_cd = NA,
                            coord_datum_cd = NA,
                            dec_coord_datum_cd = NA,
                            district_cd = NA,
                            state_cd = NA,
                            county_cd = NA,
                            country_cd = NA,
                            land_net_ds = NA,
                            map_nm = NA,
                            map_scale_fc = NA,
                            alt_va = NA,
                            alt_meth_cd = NA,
                            alt_acy_va = NA,
                            alt_datum_cd = NA,
                            huc_cd = NA,
                            basin_cd = NA,
                            topo_cd = NA,
                            instruments_cd = NA,
                            construction_dt = NA,
                            inventory_dt = NA,
                            # ws area change (hectares to square miles)
                            drain_area_va = area / 2.59,
                            contrib_drain_area_va = NA,
                            tz_cd = 'UTC',
                            local_time_fg = 'Y',
                            reliability_cd = NA,
                            gw_file_cd = NA,
                            nat_aqfr_cd = NA,
                            aqfr_cd = NA,
                            aqfr_type_cd = NA,
                            well_depth_va = NA,
                            hole_depth_va = NA,
                            depth_src_cd = NA,
                            project_no = NA,
                            shortName = site_code,
                            drainSqKm = site_ws_area,
                            staAbbrev = site_code,
                            param.nm = var,
                            param.units = var_unit,
                            paramShortName = var,
                            paramNumber = NA,
                            constitAbbrev = var,
                            paStart = 10,
                            paLong = 12)

        eList <- EGRET::mergeReport(INFO_file,
                                  Daily_file,
                                  Sample_file,
                                  verbose = TRUE)

        if(! run_egret){
            return(eList)
        }

        # TODO: print stats on the record being calculated,
        # particularly what the nuber of obs truly is (helps contextualize absurd results)
        # keeping in mind EGRET says that anything less than n=60 is "dangerous"
        eList <- try(modelEstimation(eList,
                                     minNumObs = minNumObs,
                                     minNumUncen = minNumUncen,
                                     windowY = 7,
                                     windowQ = 2,
                                     windowS = 0.5,
                                     verbose = TRUE))

        if(inherits(eList, 'try-error')){
            stop('EGRET failed while running WRTDS. See https://github.com/USGS-R/EGRET for reasons data may not be compatible with the WRTDS model.')
        }

        if(kalman){
            eList <- EGRET::WRTDSKalman(eList, verbose = !quiet)

            if(prep_data){
                # Set flux and Q to 0 and conc to NA on no flow days
                eList$Daily <- eList$Daily %>%
                    mutate(Q = ifelse(Date %in% !!no_flow_days, 0, Q),
                           LogQ = ifelse(Date %in% !!no_flow_days, 0, LogQ),
                           Q7 = ifelse(Date %in% !!no_flow_days, 0, Q7),
                           Q30 = ifelse(Date %in% !!no_flow_days, 0, Q30),
                           ConcDay = ifelse(Date %in% !!no_flow_days, NA, ConcDay),
                           FluxDay = ifelse(Date %in% !!no_flow_days, 0, FluxDay),
                           FNConc = ifelse(Date %in% !!no_flow_days, NA, FNConc),
                           FNFlux = ifelse(Date %in% !!no_flow_days, 0, FNFlux),
                           GenFlux = ifelse(Date %in% !!no_flow_days, 0, GenFlux),
                           GenConc = ifelse(Date %in% !!no_flow_days, 0, GenConc))
            }
        } else if(prep_data){
            # Set flux and Q to 0 and conc to NA on no flow days
            eList$Daily <- eList$Daily %>%
                mutate(Q = ifelse(Date %in% !!no_flow_days, 0, Q),
                       LogQ = ifelse(Date %in% !!no_flow_days, 0, LogQ),
                       Q7 = ifelse(Date %in% !!no_flow_days, 0, Q7),
                       Q30 = ifelse(Date %in% !!no_flow_days, 0, Q30),
                       ConcDay = ifelse(Date %in% !!no_flow_days, NA, ConcDay),
                       FluxDay = ifelse(Date %in% !!no_flow_days, 0, FluxDay),
                       FNConc = ifelse(Date %in% !!no_flow_days, NA, FNConc),
                       FNFlux = ifelse(Date %in% !!no_flow_days, 0, FNFlux))
        }

        # find any 'breaks' in record
        sample_rec <- detect_record_break(Sample_file)
        sample_breaks <- get_break_dates(sample_rec)

        # TODO: make dynamic in case of multiple long breaks
      if(length(sample_breaks['start'][[1]]) < 1) {
        writeLines(paste('no gap in record detected larger than', gap_period, 'days'))
      } else {

        for(i in length(sample_breaks['start'])) {
          # set period
          startBlank = sample_breaks['start'][[1]][[i]]
          endBlank = sample_breaks['end'][[1]][[i]]

          writeLines(paste('\n\nfound gap in record greater than', gap_period, 'days\n',
                           'between', startBlank, 'and', endBlank, 'masking this period',
                           'with EGRET blankTime()'))

          eList <- blankTime(eList,
                           startBlank = startBlank,
                           endBlank = endBlank
                           )
        }
      }

        return(eList)
    }

    ms_chem <- chem_df %>%
      mutate(site_code = 'none',
        # TODO: variable handling
           var = 'variable_chem',
           ms_status = 0,
           ms_interp = 0) %>%
      rename(val = con,
             datetime = all_of(datecol))

    ms_q <- q_df %>%
      mutate(site_code = 'none',
             var = 'variable_q',
             ms_status = 0,
             ms_interp = 0) %>%
      rename(val = q_lps,
             datetime = all_of(datecol))

    site_data <- tibble(site_code = 'none',
                        ws_area_ha = ws_size,
                        latitude = lat,
                        longitude = long,
                        site_type = 'stream_gauge')

  egret_results <- ms_run_egret_adapt(stream_chemistry = ms_chem,
                                      discharge = ms_q,
                                      prep_data = TRUE,
                                      site_data = site_data,
                                      kalman = kalman,
                                      run_egret = TRUE,
                                      minNumObs = minNumObs,
                                      minNumUncen = minNumUncen,
                                      gap_period = gap_period)

    return(egret_results)
}


# WRTDS DeBugging ShortCut argument defs
# run_output <- adapt_ms_egret(chem_df = ocon_chem, q_df = ocon_q, ws_size = ocon_ws_area_ha)
## EGRET::plotConcHist(egret_results)
## EGRET::plotConcTime(egret_results)
## EGRET::plotFluxQ(egret_results)
## EGRET::plotConcPred(egret_results)

# re-run defnitions cheatsheet

## # mas adapt egret
## chem_df = chem_df
## q_df = q_df
## ws_size = area
## lat = lat
## long = long

## ## pre egret
##   ms_chem <- chem_df %>%
##     mutate(site_code = 'none',
##            var = 'IS_NO3',
##            ms_status = 0,
##            ms_interp = 0) %>%
##       rename(val = con,
##              datetime = all_of(datecol))

##     ms_q <- q_df %>%
##       mutate(site_code = 'none',
##              var = 'IS_discharge',
##              ms_status = 0,
##              ms_interp = 0) %>%
##       rename(val = q_lps,
##              datetime = all_of(datecol))

## site_data <- tibble(site_code = 'none',
##   ws_area_ha = ws_size,
##   latitude = lat,
##   longitude = long,
##   site_type = 'stream_gauge')

## stream_chemistry = ms_chem
## discharge = ms_q
## prep_data = TRUE
## site_data = site_data
## kalman = FALSE
## run_egret = TRUE
## quiet = FALSE

## minNumObs = 100
## minNumUncen = 50
## gap_period = 730

## verbose = TRUE

## windowY = 7
## windowQ = 2
## windowS = 0.5
## edgeAdjust = TRUE
## verbose = TRUE
## run.parallel = FALSE

# surfaces
## surfaceStart=NA
## surfaceEnd=NA
## localSample=NA

report_on_df <- function(data) {
  water_years <- unique(data$wy)
  for(watyr in water_years) {
    data_wy <- data %>%
      filter(wy == watyr)
    nsamples <- nrow(data_wy)

    writeLines(paste('Water Year:', watyr,
                     '\n     number of samples:', nsamples))
  }
}

flux_na_reporter <- function(data) {
      print('rows with NA flux day values:')
      print(data[is.na(data$FluxDay),])
      print('rows with Inf flux day values:')
      print(data[is.infinite(data$FluxDay),])

      print('returning DF of all Inf and NA rows')
      data_error <- data[is.infinite(data$FluxDay) | is.na(data$FluxDay),]
      return(data_error)
}
