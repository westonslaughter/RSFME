data_dir <- 'streamlined/data/ms/hbef/'
site_files  <- list.files('streamlined/data/ms/hbef/discharge', recursive = F)
site_info  <- read_csv('streamlined/data/site/ms_site_info.csv')


out_frame <- tibble(wy = as.integer(), 
                    site_code = as.character(),
                    var = as.character(),
                    val = as.numeric(), 
                    method = as.character(),
                    ms_reccomended = as.integer())

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
    
    raw_data_con <- read_feather(here(glue(data_dir, 'stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0)
    
    raw_data_q <- read_feather(here(glue(data_dir, 'discharge/', site_code, '.feather')))
    
    # initialize next loop
    solutes <- raw_data_con %>%
        filter(!str_detect(var, 'Charge'),
               !str_detect(var, 'temp'),
               !str_detect(var, 'pH')) %>%
        select(var) %>%
        unique() %>%
        pull(var)

## Loop through solutes at site #####
for(j in 1:length(solutes)){
    
    #set to target solute
    target_solute <- solutes[j]
    
    raw_data_con <- read_feather(here(glue(data_dir, 'stream_chemistry/', site_code, '.feather'))) %>%
        filter(ms_interp == 0,
               val > 0)
    
    raw_data_con <- raw_data_con %>%
        filter(var == target_solute) %>%
        select(datetime, val) %>%
        na.omit()
    
    # find acceptable years
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
        mutate(water_year = water_year(date, origin = "usgs")) %>%
        group_by(water_year) %>%
        summarise(n = n()) %>%
        filter(n >= 4)
    
    q_good_years <- q_check$water_year
    conc_good_years <- conc_check$water_year
    
    good_years <- q_good_years[q_good_years %in% conc_good_years]
    
    #join data and cut to good years
    daily_data_con <- raw_data_con %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean(val)) %>%
        mutate(site_code = !!site_code, var = 'con') %>%
        select(site_code, datetime = date, var, val)
    
    daily_data_q <- raw_data_q %>%
        mutate(date = date(datetime)) %>%
        group_by(date) %>%
        summarize(val = mean(val)) %>%
        mutate(site_code = !!site_code, var = 'q_lps') %>%
        select(site_code, datetime = date, var, val)
    
    raw_data_full <- rbind(daily_data_con, daily_data_q) %>%
        pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
        mutate(wy = water_year(datetime, origin = 'usgs')) %>%
        filter(wy %in% good_years)
    
    ### Loop through good years #####
    for(k in 1:length(good_years)){
        target_year <- as.numeric(as.character(good_years[k]))
        
        raw_data_target_year <- raw_data_full %>%
            mutate(wy = as.numeric(as.character(wy))) %>%
            filter(wy == target_year)
        
        q_target_year <- raw_data_target_year %>%
            select(site_code, datetime, q_lps, wy)%>%
            na.omit()
        
        con_target_year <- raw_data_target_year %>%
            select(site_code, datetime, con, wy) %>%
            na.omit()
    
        #### Riverload conversion function #####
        prep_raw_for_riverload <- function(chem_df, q_df){
            conv_q <- q_df %>%
                mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
                mutate(flow = q_lps*0.001) %>% # convert lps to cubic meters per second)
                select(datetime, flow) %>%
                data.frame()
            
            conv_c <- chem_df %>%
                mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')) %>%
                select(datetime, con) %>%
                data.frame()
            
            db <- full_join(conv_q, conv_c, by = "datetime") %>%
                #filter(!is.na(flow)) %>%
                arrange(datetime)
            
            return(db)
        }
        
        ### calculate annual flux ######
        chem_df <- con_target_year
        q_df <- q_target_year
        
        #### calculate average ####
        flux_annual_average <- raw_data_target_year %>%
            group_by(wy) %>%
            summarize(q_lps = mean(q_lps, na.rm = TRUE),
                      con = mean(con, na.rm = TRUE)) %>%
            mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
            pull(flux)
            
        
        #### calculate period weighted #####
        calculate_pw <- function(chem_df, q_df){
            rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
            
            flux_from_pw <- method1(rl_data, ncomp = 1) %>%
                sum(.)/(1000*area)
            return(flux_from_pw)
        }
        flux_annual_pw <- calculate_pw(chem_df, q_df)
        
        #### calculate beale ######
        calculate_beale <- function(chem_df, q_df){
            rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
            flux_from_beale <- beale.ratio(rl_data, ncomp = 1) %>%
                sum(.)/(1000*area)
            return(flux_from_beale)
        }
        flux_annual_beale <- calculate_beale(chem_df, q_df)

        #### calculate rating #####
        calculate_rating <- function(chem_df, q_df){
            rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
            flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1) %>%
                sum(.)/(1000*area)
            return(flux_from_reg)
        }
        flux_annual_rating <- calculate_rating(chem_df, q_df)
        
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
        
        
        #### calculate composite ######
        generate_residual_corrected_con <- function(chem_df, q_df){
            # first make c:q rating
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
            
            # extract model info
            intercept <- rating$coefficients[1]
            slope <- rating$coefficients[2]
            
            # create modeled c, calc residuals, adjust modeled c by interpolated residuals
            rating_filled_df <- q_df %>%
                mutate(con_reg = 10^(intercept+(slope*log10(q_lps)))) %>%
                select(datetime, con_reg, q_lps) %>%
                full_join(., chem_df, by = 'datetime') %>%
                select(site_code, datetime, con, con_reg, q_lps, wy) %>%
                mutate(res = con_reg-con,
                       res = imputeTS::na_interpolation(res),
                       con_com = con_reg-res,
                       site_code = !!site_code,
                       wy = water_year(datetime, origin = 'usgs'))
            rating_filled_df$con_com[!is.finite(rating_filled_df$con_com)] <- 0
            return(rating_filled_df)
        }
        
        rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                            q_df = q_df)
        
        # calculate annual flux from composite
        calculate_composite_from_rating_filled_df <- function(rating_filled_df){
            flux_from__comp <- rating_filled_df %>%
                mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
                group_by(wy) %>%
                summarize(flux = sum(flux)) %>%
                mutate(site_code = !!site_code)
        }
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
        
        rating <- summary(lm(model_data$c_log ~ model_data$q_log, singular.ok = TRUE))
        
        r_squared <- rating$r.squared
        
        resid_acf <- abs(acf(rating$residuals, lag.max = 1, plot = FALSE)$acf[2])
        
        con_acf <- abs(acf(paired_df$con, lag.max = 1, plot = FALSE)$acf[2])
        
        # modified from figure 10 of Aulenbach et al 2016
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
        
        #### congeal fluxes ####
        target_year_out <- tibble(wy = as.character(target_year), 
                            val = c(flux_annual_average, flux_annual_pw, flux_annual_beale, 
                                     flux_annual_rating, flux_annual_comp$flux[1]),
                            site_code = !!site_code, 
                            var = !!target_solute,
                            method = c('average', 'pw', 'beale', 'rating', 'composite')) %>%
            mutate(ms_recommended = ifelse(method == !!ideal_method, 1, 0))
        out_frame <- rbind(out_frame, target_year_out)
        
        } # end year loop
    } # end solute loop
    directory <- glue(data_dir,'stream_flux/')
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
    }
    file_path <- glue('{directory}/{s}.feather',
                      s = site_code)
write_feather(out_frame, file_path)
} # end site loop
