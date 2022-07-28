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


 warn_sum <- function(x) {
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
