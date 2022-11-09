# simple linear interpolation between con samples at frequency of q
# in:
# out:
interpolate_con_to_q <- function(chem_df, q_df){
 # Join chem and q
 join_df <- full_join(chem_df, q_df, by = 'date') %>%
   arrange(date) %>%
   mutate(interp = 0) %>%
    select(date, con, q_lps, interp)
 join_df$interp[is.na(join_df$con)] <- 1

 #interpolate chem to match q
 interp_df <- na_interpolation(join_df, option = "linear", maxgap = Inf) %>%
     mutate(wy = water_year(date, origin = 'usgs'))
 return(interp_df)
}

# take usgs data and make it ready for riverload functions
# in:
# out:
prep_usgs_for_riverload <- function(chem_df, q_df){
    conv_q <- q_df %>%
        mutate(site_code = !!site_code,
               flow = Q/35.3147) %>% # convert cfs to cubic meters per second)
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


# datetime to water year quarters
dt_to_wy_quarter <- function(datetime) {
  tryCatch(
    expr = {
      date <- date(as.Date(datetime))

      # divide WY into quarters
      m <- as.integer(format(date, "%m"))

      if(m %in% c(10, 11, 12)) {
        return("Q1")
      } else if(m %in% c(1, 2, 3)) {
        return("Q2")
      } else if(m %in% c(4, 5, 6)) {
        return("Q3")
      } else if(m %in% c(7, 8, 9)) {
        return("Q4")
      }
    },
    error = function(e) {
      print("something went wrong in dt_to_wy_quarters() conversion")
    }
  )
}


 old_warn_sum <- function(x) {
   # use to sitll get data with NAs and Infs
   # length of record
   n_x <- length(x)

   if(TRUE %in% is.na(x)) {
     # number of NAs in record
     n_na <- length(x[is.na(x)])

    if(TRUE %in% is.infinite(x)) {
      # number of inf values in record
      n_inf <- length(x[is.infinite(x)])

      # warn user
      writeLines(paste("WARNING: infinite values found in flux record, ignoring during SUM",
                       "\n count NA:", n_na,
                       "\n count Inf:", n_inf))


      xsum <- sum(x[!is.infinite(x)], na.rm = TRUE)
      return(xsum)
    }

    xsum <- sum(x, na.rm = TRUE)

    writeLines(paste('NAs found in flux record, ignoring during SUM. \n count NA:',
                 n_na,
                 '\n percent NA:', (n_na/n_x) * 100))
   } else {
     xsum <- sum(x)
   }
   return(xsum)
 }

 warn_sum <- function(x) {
   # length of record
   n_x <- length(x)

     # number of NAs in record
     n_na <- length(x[is.na(x)])
      # number of inf values in record
      n_inf <- length(x[is.infinite(x)])

      # warn user
      writeLines(paste("WARNING: infinite values found in flux record, ignoring during SUM",
                       "\n count NA:", n_na,
                       "\n count Inf:", n_inf))


      xsum <- sum(x[!is.infinite(x)], na.rm = TRUE)
      return(xsum)
 }

# run flux by site, year(s), method

## ms_run_flux <- function(site, )

##### calculate ms_interp #####

 carry_flags <- function(raw_q_df, raw_con_df, target_solute = NULL, target_year = NULL, period = NULL){
    #### set up ratio functions #####
         # interp ratio
         calc_interp_ratio <- function(trimmed_df, period = NULL){

            if(period == 'annual'){
            no_interp <- trimmed_df %>%
                filter(ms_interp == 0) %>%
                count() %>%
                pull(n)

            interp <- trimmed_df %>%
                filter(ms_interp ==1) %>%
                count() %>%
                pull(n)

            interp_ratio <- interp/(interp+no_interp)

            return(interp_ratio)

            }else if (period == 'month'){
             interp_ratio <- trimmed_df %>%
                    mutate(month = month(datetime)) %>%
                    group_by(month) %>%
                    summarize(n_interp = sum(ms_interp == 1),
                              n_no_interp = sum(ms_interp == 0)) %>%
                    mutate(interp_ratio = n_interp/(n_interp+n_no_interp)) %>%
                    select(month, interp_ratio)

             return(interp_ratio)

            }else{print('Specify period as month or annual.')}

         }

         # status ratio
         calc_status_ratio <- function(trimmed_df, period = NULL){

             if(period == 'annual'){
             no_status <- trimmed_df %>%
                 filter(ms_status == 0) %>%
                 count() %>%
                 pull(n)

             status <- trimmed_df %>%
                 filter(ms_interp == 1) %>%
                 count() %>%
                 pull(n)

             status_ratio <- status/(status+no_status)

             return(status_ratio)}
             else if (period == 'month'){
                 status_ratio <- trimmed_df %>%
                     mutate(month = month(datetime)) %>%
                     group_by(month) %>%
                     summarize(n_stat = sum(ms_status == 1),
                               n_no_stat = sum(ms_status == 0)) %>%
                     mutate(status_ratio = n_stat/(n_stat+n_no_stat)) %>%
                     select(month, status_ratio)

                 return(status_ratio)

             }else{print('Specify period as month or annual.')}
         }

         # missing ratio
         calc_missing_ratio <- function(trimmed_df, period = NULL){

             if(period == 'annual'){
                present <- trimmed_df %>%
                    select(val) %>%
                    na.omit() %>%
                    nrow()

                missing_ratio <- (365-present)/365

                return(missing_ratio)

             }else if(period == 'month'){
                 missing_ratio <- trimmed_df %>%
                     mutate(month = month(datetime)) %>%
                     group_by(month) %>%
                     select(datetime, val) %>%
                     na.omit() %>%
                     summarize(n = n()) %>%
                     mutate(full_days = days_in_month(month),
                            missing_ratio = (full_days-n)/full_days) %>%
                     select(month, missing_ratio)

                 return(missing_ratio)

             }else{
                 print('Specify period as month or annual.')
             }

         }
    ### run functions on data ####
         if(period == 'annual'){

             year_con_df <- raw_con_df %>%
                 mutate(wy = water_year(datetime, origin = 'usgs')) %>%
                 filter(wy == target_year,
                        target_solute == target_solute)

             con_tbl <- tibble(wy = target_year, var = target_solute,
                               ms_interp = calc_interp_ratio(trimmed_df = year_con_df, period = period),
                               ms_status = calc_status_ratio(trimmed_df = year_con_df, period = period),
                               ms_missing = calc_missing_ratio(trimmed_df = year_con_df, period = period),
                               )

             year_q_df <- raw_q_df %>%
                 mutate(wy = water_year(datetime, origin = 'usgs')) %>%
                 filter(wy == target_year)

             q_tbl <- tibble(wy = target_year, var = unique(year_q_df$var)[1],
                               ms_interp = calc_interp_ratio(trimmed_df = year_q_df, period = period),
                               ms_status = calc_status_ratio(trimmed_df = year_q_df, period = period),
                               ms_missing = calc_missing_ratio(trimmed_df = year_q_df, period = period)
                             )

             out_frame <- rbind(con_tbl, q_tbl)
            return(out_frame)
         } else

         if(period == 'month'){

             year_con_df <- raw_con_df %>%
                 mutate(wy = water_year(datetime, origin = 'usgs')) %>%
                 filter(wy == target_year,
                        target_solute == target_solute)

            ms_interp = calc_interp_ratio(trimmed_df = year_con_df, period = period)
            ms_status = calc_status_ratio(trimmed_df = year_con_df, period = period)
            ms_missing = calc_missing_ratio(trimmed_df = year_con_df, period = period)

            con_tbl <- full_join(ms_interp, ms_status, by = 'month') %>%
                full_join(ms_missing, by = 'month') %>%
                mutate(wy = target_year,
                       var = target_solute)

             year_q_df <- raw_q_df %>%
                 mutate(wy = water_year(datetime, origin = 'usgs')) %>%
                 filter(wy == target_year)

             ms_interp = calc_interp_ratio(trimmed_df = year_q_df, period = period)
             ms_status = calc_status_ratio(trimmed_df = year_q_df, period = period)
             ms_missing = calc_missing_ratio(trimmed_df = year_q_df, period = period)

             q_tbl <- full_join(ms_interp, ms_status, by = 'month') %>%
                 full_join(ms_missing, by = 'month') %>%
                 mutate(wy = target_year,
                        var = unique(year_q_df$var)[1])

             out_frame <- rbind(con_tbl, q_tbl)
             return(out_frame)
         }

 }



