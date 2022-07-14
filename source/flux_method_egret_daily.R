adapt_ms_egret <- function(chem_df, q_df, ws_size, lat, long, site_data = NULL, kalman = FALSE){
  
    get_MonthSeq <- function(dates){
        
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
        
        # Set up EGRET Sample file 
        Sample_file <- tibble(Name = site_code,
                              Date = as.Date(stream_chemistry$datetime),
                              ConcLow = stream_chemistry$val,
                              ConcHigh = stream_chemistry$val,
                              Uncen = 1,
                              ConcAve = stream_chemistry$val,
                              Julian = as.numeric(julian(lubridate::ymd(stream_chemistry$datetime),origin=as.Date("1850-01-01"))),
                              Month = lubridate::month(stream_chemistry$datetime),
                              Day = lubridate::yday(stream_chemistry$datetime),
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
                             Q = discharge_daily$val/1000,
                             Julian = as.numeric(julian(lubridate::ymd(discharge_daily$datetime),origin=as.Date("1850-01-01"))),
                             Month = lubridate::month(discharge_daily$datetime),
                             Day = lubridate::yday(discharge_daily$datetime),
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
            
        }
        
        Daily_file <- Daily_file %>%
            mutate(Q7 = zoo::rollmean(Q, 7, fill = NA, align = 'right'),
                   Q30 = zoo::rollmean(Q, 30, fill = NA, align = 'right'),
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
                            drain_area_va = site_ws_area*2.59,
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
                                    verbose = !quiet)
        
        if(! run_egret){
            return(eList)
        }
        
        eList <- try(EGRET::modelEstimation(eList,
                                        minNumObs = 4,
                                        minNumUncen = 4, verbose = !quiet))
        
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
           var = 'IS_NO3',
           ms_status = 0,
           ms_interp = 0) %>%
      rename(val = con,
             datetime = date)
    
    ms_q <- q_df %>%
      mutate(site_code = 'none',
             var = 'IS_discharge',
             ms_status = 0,
             ms_interp = 0) %>%
      rename(val = q_lps,
             datetime = date)
    
    site_data <- tibble(site_code = 'none',
                        ws_area_ha = ws_size,
                        latitude = lat,
                        longitude = long,
                        site_type = 'stream_gauge')
    
    egret_results <- ms_run_egret_adapt(stream_chemistry = ms_chem, discharge = ms_q,
                                        prep_data = TRUE, site_data = site_data,
                                        kalman = kalman, run_egret = FALSE)
    
    eList <- try(EGRET::modelEstimation(egret_results, verbose = TRUE,
                                        minNumObs = 4,
                                        minNumUncen = 4))
    
    return(eList)
}

# run_output <- adapt_ms_egret(chem_df = ocon_chem, q_df = ocon_q, ws_size = ocon_ws_area_ha)
# EGRET::plotConcHist(run_output)
# EGRET::plotConcTime(run_output)
# EGRET::plotFluxQ(run_output)
# EGRET::plotConcPred(run_output)
