# Flux Prep and Helpers

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

###### calculate period weighted#########
calculate_pw <- function(chem_df, q_df){
rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df
flux_from_pw <- method1(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
return(flux_from_pw)

###### calculate beale ######
calculate_beale <- function(chem_df, q_df){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
flux_from_beale <- beale.ratio(rl_data, ncomp = 1) %>%
    sum(.)/(1000*area)
return(flux_from_beale)

##### calculate rating #####
calculate_rating <- function(chem_df, q_df){
    rl_data <- prep_raw_for_riverload(chem_df = chem_df, q_df = q_df)
    flux_from_reg <- RiverLoad::rating(rl_data, ncomp = 1) %>%
        sum(.)/(1000*area)
    return(flux_from_reg)


##### calculate WRTDS #####
calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long)
  tryCatch(
    expr = {
      egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long)

              # still looking for reason why wrtds is 1K higher than others
              flux_from_egret <- egret_results$Daily$FluxDay %>%
                warn_sum(.)/(area)

            },
            error = function(e) {
              print('ERROR: WRTDS failed to run')
        })
    return(flux_from_egret)
  }

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

# calculate annual flux from composite
calculate_composite_from_rating_filled_df <- function(rating_filled_df){
        flux_from__comp <- rating_filled_df %>%
            mutate(flux = con_com*q_lps*86400*(1/area)*1e-6) %>%
            group_by(wy) %>%
            summarize(flux = sum(flux)) %>%
            mutate(site_code = site_no)

        return(flux_from__comp)
        }
