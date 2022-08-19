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

  if(is.null(period)){
  flux_from_pw <- method1(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
  }else{

  if(period == 'month'){

      method1_month <- function (db, ncomp, period)
      {
          notNA <- na.omit(db)
          if (missing(ncomp)) {
              print("ncomp is missing.")
          }

          if (period == "month") {
              new <- notNA
              new$newdate <- format(as.POSIXct(new$datetime), format = "%Y-%m")
              maximum <- length(unique(new$newdate))
              index <- vector(length = maximum)
              for (i in 1:(maximum)) {
                  index[i] <- length(which(new$newdate == (unique(new$newdate)[i])))
              }
              result <- vector("list", maximum)
              for (i in 1:(maximum)) {
                  seldata <- subset(new, new$newdate == unique(new$newdate)[i])
                  result[[i]] <- seldata[, -which(names(seldata) %in%
                                                      c("datetime", "flow", "newdate"))]/index[i]
              }
              if (is.null(nrow(result[[1]])) == T) {
                  mat.conc <- as.matrix(unlist(result))
              }
              if (is.null(nrow(result[[1]])) == F) {
                  mat.conc <- do.call(rbind, result)
              }
              concdate <- cbind.data.frame(notNA$datetime, mat.conc)
              colnames(concdate)[1] <- c("datetime")
              concdate$newdate <- format(as.POSIXct(concdate$datetime),
                                         format = "%Y-%m")
              aggrg.data <- matrix(nrow = length(unique(concdate$newdate)),
                                   ncol = (ncomp))
              for (i in 1:(ncomp)) {
                  agg.init <- aggregate(concdate[, i + 1] ~ newdate,
                                        concdate, sum)
                  aggrg.data[, i] <- agg.init[, 2]
                  colnames(aggrg.data) <- c(names(db)[3:(ncomp + 2)])
                  if (i == (ncomp + 1))
                      break
              }
              resultF <- vector("list", maximum)
              for (i in 1:(maximum)) {
                  seldata <- subset(new, new$newdate == unique(new$newdate)[i])
                  resultF[[i]] <- seldata$flow/index[i]
              }
              flow.n <- matrix(unlist(resultF), ncol = 1, byrow = TRUE)
              colnames(flow.n) <- c("flow")
              flowdate <- cbind.data.frame(notNA$datetime, flow.n)
              colnames(flowdate) <- c("datetime", "flow")
              flowdate$newdate <- format(as.POSIXct(flowdate$datetime),
                                         format = "%Y-%m")
              aggrg.data1 <- (aggregate(flow ~ newdate, flowdate,
                                        sum))
              aggrg.data2 <- (aggregate(flow ~ newdate, flowdate,
                                        sum))
              aggrg.data2$newdate <- as.Date(paste(aggrg.data2$newdate,
                                                   "-01", sep = ""))
              aggrg.data2$newdate <- as.POSIXct(aggrg.data2$newdate,
                                                format = "%Y-%m")
              aggrg.flow <- (aggrg.data2[, -which(names(aggrg.data2) %in%
                                                      c("newdate"))])
              load <- (aggrg.flow * aggrg.data)
              n <- nrow(db)
              datemonth <- seq(as.POSIXct(db$datetime[1], tz = "CET"),
                               as.POSIXct(db$datetime[n], tz = "CET"), "days")
              b <- length(datemonth)
              dateplus <- as.Date(datemonth[b]) + 2
              daymonths <- as.numeric(round(diff(seq(as.POSIXct(datemonth[1],
                                                                tz = "CET"), as.POSIXct(dateplus, tz = "CET"), "month")),
                                            digits = 0))
              if(length(good_months) < 12){
                  daymonths <- daymonths[good_months]
              }
              method1 <- (load * (daymonths) * 86400)
              rownames(method1) <- ((aggrg.data1)[, 1])
              return(method1)
          }
          else if (period == "year") {
              new <- notNA
              new$newdate <- format(as.POSIXct(new$datetime), format = "%Y")
              maximum <- length(unique(new$newdate))
              index <- vector(length = maximum)
              for (i in 1:(maximum)) {
                  index[i] <- length(which(new$newdate == (unique(new$newdate)[i])))
              }
              result <- vector("list", maximum)
              for (i in 1:(maximum)) {
                  seldata <- subset(new, new$newdate == unique(new$newdate)[i])
                  result[[i]] <- seldata[, -which(names(seldata) %in%
                                                      c("datetime", "flow", "newdate"))]/index[i]
              }
              if (is.null(nrow(result[[1]])) == T) {
                  mat.conc <- as.matrix(unlist(result))
              }
              if (is.null(nrow(result[[1]])) == F) {
                  mat.conc <- do.call(rbind, result)
              }
              concdate <- cbind.data.frame(notNA$datetime, mat.conc)
              colnames(concdate)[1] <- c("datetime")
              concdate$newdate <- format(as.POSIXct(concdate$datetime),
                                         format = "%Y")
              aggrg.data <- matrix(nrow = length(unique(concdate$newdate)),
                                   ncol = (ncomp))
              for (i in 1:(ncomp)) {
                  agg.init <- aggregate(concdate[, i + 1] ~ newdate,
                                        concdate, sum)
                  aggrg.data[, i] <- agg.init[, 2]
                  colnames(aggrg.data) <- c(names(db)[3:(ncomp + 2)])
                  if (i == (ncomp + 1))
                      break
              }
              resultF <- vector("list", maximum)
              for (i in 1:(maximum)) {
                  seldata <- subset(new, new$newdate == unique(new$newdate)[i])
                  resultF[[i]] <- seldata$flow/index[i]
              }
              flow.n <- matrix(unlist(resultF), ncol = 1, byrow = TRUE)
              colnames(flow.n) <- c("flow")
              flowdate <- cbind.data.frame(notNA$datetime, flow.n)
              colnames(flowdate) <- c("datetime", "flow")
              flowdate$newdate <- format(as.POSIXct(flowdate$datetime),
                                         format = "%Y")
              aggrg.data2 <- (aggregate(flow ~ newdate, flowdate,
                                        sum))
              aggrg.flow <- (aggrg.data2[, -which(names(aggrg.data2) %in%
                                                      c("newdate"))])
              load <- (aggrg.flow * aggrg.data)
              is.leapyear = function(year) {
                  return(((year%%4 == 0) & (year%%100 != 0)) | (year%%400 ==
                                                                    0))
              }
              for (i in 1:(nrow(aggrg.data2))) {
                  if (is.leapyear(as.numeric(aggrg.data2$newdate[i])) ==
                      T) {
                      method1 <- (load * (366) * 86400)
                  }
                  else {
                      method1 <- (load * (365) * 86400)
                  }
              }
              rownames(method1) <- aggrg.data2[, 1]
              colnames(method1) <- c(names(db)[3:(ncomp + 2)])
              return(method1)
          }
      }

      ##### apply #####

      flux_from_pw <- method1_month(rl_data, ncomp = 1, period = period)

      flux_from_pw <- tibble(date = rownames(flux_from_pw),
                          flux = (flux_from_pw[,1]/(1000*area)))


  }
  }

  return(flux_from_pw)
}

###### calculate beale ######
calculate_beale <- function(chem_df, q_df, datecol = 'date', period = NULL){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)

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
calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long, datecol = 'date', agg = 'default') {
  tryCatch(
    expr = {
      # default sums all daily flux values in df
      if(agg == 'default') {
        egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long, datecol = datecol)

        flux_from_egret <- egret_results$Daily$FluxDay %>%
          warn_sum(.)/(area)
      } else if(agg == 'annual') {
        egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long, datecol = datecol)
        flux_from_egret <- egret_results$Daily %>%
                    mutate(
                      wy = water_year(Date)
                    ) %>%
          group_by(wy) %>%
          summarize(
            flux = warn_sum(FluxDay)/(area)
          )
      } else if(agg == 'monthly') {
        egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long, datecol = datecol)
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
          select(all_of(datecol), con_reg, q_lps) %>%
          full_join(., chem_df, by = datecol) %>%
          select(site_code, all_of(datecol), con, con_reg, q_lps, wy)  %>%
            mutate(res = con_reg-con,
                   res = imputeTS::na_interpolation(res),
                   con_com = con_reg-res,
                   site_code = !!get(sitecol),
                   wy = water_year(get(datecol), origin = 'usgs'))

        rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0
        return(rating_filled_df)
        }

##### calculate monthly flux from composite ####
calculate_composite_from_rating_filled_df <- function(rating_filled_df, site_no = 'site_no', period = NULL){

        if(is.null(period)){
        flux_from__comp <- rating_filled_df %>%
            mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
            group_by(wy) %>%
          summarize(flux = sum(flux)) %>%
            mutate(site_code = site_no)
        }else{

        if(period == 'month'){
            flux_from__comp <- rating_filled_df %>%
                mutate(month = month(datetime),
                       flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
                group_by(wy, month) %>%
                summarize(date = max(datetime),
                          flux = sum(flux))
        }
        }

        return(flux_from__comp)
        }

# wrapper for WRTDS in EGRET
adapt_ms_egret <- function(chem_df, q_df, ws_size, lat, long,
                           site_data = NULL, kalman = FALSE,
                           datecol = 'date', minNumObs = 2, minNumUncen = 2){
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
             site_data = NULL, min_q_method = 'USGS'){

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
            mean_flow <- mean(Daily_file$Q[Daily_file$Q > 0], na.rm = TRUE)

            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!mean_flow, Q))

          } else {
            # NOTE: could this be where Inf shows up too? like in min_chem?
            min_flow <- min(Daily_file$Q[Daily_file$Q > 0], na.rm = TRUE)

            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!min_flow, Q))
          }



            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!min_flow, Q))
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
                                      run_egret = TRUE)

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
## q_df = prep_data_q
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

## minNumObs = 2
## minNumUncen = 2
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
