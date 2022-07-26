# Flux Functions
# PREP and HELPERS
### Riverload conversion function #####
prep_raw_for_riverload <- function(chem_df, q_df, datecol = 'date'){
    conv_q <- q_df %>%
                mutate(datetime = as.POSIXct(get(datecol), format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
                mutate(flow = q_lps*0.001) %>% # convert lps to cubic meters per second)
                select(datetime, flow) %>%
                data.frame()

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
calculate_pw <- function(chem_df, q_df, datecol = 'date'){
  rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)
  flux_from_pw <- method1(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
  return(flux_from_pw)
}

###### calculate beale ######
calculate_beale <- function(chem_df, q_df, datecol = 'date'){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)
    flux_from_beale <- beale.ratio(rl_data, ncomp = 1) %>%
      sum(.)/(1000*area)
    return(flux_from_beale)
}

##### calculate rating #####
calculate_rating <- function(chem_df, q_df, datecol = 'date'){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df, datecol = datecol)
    flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1) %>%
        sum(.)/(1000*area)
    return(flux_from_reg)
}


##### calculate WRTDS #####
calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long, datecol = 'date') {
  tryCatch(
    expr = {
        egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long, datecol = datecol)

        # still looking for reason why wrtds is 1K higher than others
        flux_from_egret <- egret_results$Daily$FluxDay %>%
          warn_sum(.)/(area)
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

# calculate annual flux from composite
calculate_composite_from_rating_filled_df <- function(rating_filled_df, sitecol = 'site_no'){
        flux_from__comp <- rating_filled_df %>%
            mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
            group_by(wy) %>%
            summarize(flux = sum(flux)) %>%
            mutate(site_code = all_of(sitecol))

        return(flux_from__comp)
        }

# wrapper for WRTDS in EGRET
adapt_ms_egret <- function(chem_df, q_df, ws_size, lat, long, site_data = NULL, kalman = FALSE, datecol = 'date'){
  # TODO:  reorder site data args to make fully optyional

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
             site_data = NULL){

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
            min_chem <- stream_chemistry %>%
                filter(val > 0) %>%
                pull(val) %>%
                min()

            if(!min_chem == Inf){
                stream_chemistry <- stream_chemistry %>%
                    mutate(val = ifelse(val == 0, !!min_chem, val))
            }

            # Filter so there is only Q going into the model that also has chem
            stream_chemistry <- stream_chemistry %>%
                mutate(year = lubridate::year(datetime),
                       month = lubridate::month(datetime)) %>%
                mutate(waterYear = ifelse(month %in% c(10, 11, 12), year+1, year))

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

      ## if(datamode == 'ms') {
      ##   stream_chemistry <- stream_chemistry %>%
      ##     # convert g to kg
      ##       mutate(val = val / 1000000)
      ## }

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
            min_flow <- min(Daily_file$Q[Daily_file$Q > 0], na.rm = TRUE)

            no_flow_days <- Daily_file %>%
                filter(Q == 0) %>%
                pull(Date)

            Daily_file <- Daily_file %>%
                mutate(Q = ifelse(Q <= 0, !!min_flow, Q))
          # TODO: record zero flow days, and set flux for those days to zero

        }

        Daily_file <- Daily_file %>%
            # 'extend' rollmean
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

        eList <- EGRET::mergeReport(INFO_file, Daily_file, Sample_file,
                                    verbose = TRUE)

        if(! run_egret){
            return(eList)
        }

        eList <- try(modelEstimation(eList,
                                        minNumObs = 2,
                                     minNumUncen = 2,
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
           var = 'IS_NO3',
           ms_status = 0,
           ms_interp = 0) %>%
      rename(val = con,
             datetime = all_of(datecol))

    ms_q <- q_df %>%
      mutate(site_code = 'none',
             var = 'IS_discharge',
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
## ## q_df = prep_data_q
## q_df = q_df
## ws_size = area
## lat = lat
## long = long

## pre egret
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

## prep_data = TRUE
## run_egret = TRUE
## kalman = FALSE
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
