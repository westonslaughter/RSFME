library(tidyverse)
library(forecast)
library(feather)
library(xts)
library(imputeTS)
library(here)
library(lfstat)
library(lubridate)
library(ggpubr)
library(patchwork)

set.seed(53045)


source(here('source/flux_methods.R'))

area <- 42.4
site_code = 'w3'

# read in data ####
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))
#slice(1:ts_len)

# subset to 2016 wy
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy) %>%
    mutate(IS_discharge = na.approx(IS_discharge),
           IS_NO3 = na.approx(IS_NO3),
           IS_FDOM = na.approx(IS_FDOM),
           IS_spCond = na.approx(IS_spCond))

# calculate truth ####
chem_df <- dn %>%
    select(date, con = IS_NO3) %>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              con = mean(con)) %>%
    ungroup() %>%
    unique() %>%
    select(date, con) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- dn %>%
    select(date, q_lps = IS_discharge)%>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    mutate(site_code = 'w3', wy = target_wy)

    out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)
truth <- tibble(method = 'truth', estimate = out_val)

# make gradually coarsened chem ####
nth_element <- function(vector, starting_position, n) {
    vector[seq(starting_position, length(vector), n)]
}

coarse_chem <- list()
loopid = 0
for(i in seq(from = 1, to = nrow(dn)/2, by = 4)){
#for(i in (1:186)^2){
loopid <- loopid+1
n = floor(nrow(dn)/i)
coarse_chem[[loopid]] <- tibble(date =  nth_element(dn$datetime, 1, n = n),
                                con = nth_element(dn$IS_NO3, 1, n = n))
names(coarse_chem)[loopid] <- paste0('sample_',n)
}

# apply flux methods to each #####
apply_methods_coarse <- function(chem_df, q_df){
    out <- tibble(method = as.character(), estimate = as.numeric())
    #pw
    out[1,2] <- calculate_pw(chem_df, q_df)
    #beale
    out[2,2] <- calculate_beale(chem_df, q_df)
    #rating
    out[3,2] <- calculate_rating(chem_df, q_df)
    #comp
    out[4,2] <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)

    out$method <- c('pw', 'beale', 'rating', 'composite')
    return(out)
}

out_tbl <- tibble(method = as.character(), estimate = as.numeric(), n = as.integer())
for(i in 2:length(coarse_chem)){

    n <- as.numeric(str_split_fixed(names(coarse_chem[i]), pattern = 'sample_', n = 2)[2])

    chem_df <- coarse_chem[[i]] %>%
        group_by(lubridate::yday(date)) %>%
        summarize(date = date(date),
                  con = mean(con)) %>%
        ungroup() %>%
        unique() %>%
        select(date, con) %>%
        mutate(site_code = 'w3', wy = target_wy)

    out_tbl <- apply_methods_coarse(chem_df, q_df) %>%
        mutate(n = n) %>%
        rbind(., out_tbl)
}


plot_tbl <- out_tbl %>%
    unique() %>%
    mutate(error = ((estimate-truth$estimate[1])/truth$estimate[1])*100,
           error_abs = abs(error),
           method = factor(method, levels = c('pw', 'beale', 'rating', 'composite')),
           percent_coverage = (nrow(dn)/n)/nrow(dn),
           hours = n/4)

plot_tbl %>%
ggplot(., aes(x = hours, y = error_abs))+
    geom_line()+
    geom_point()+
    facet_wrap(vars(method), ncol = 1)+
    scale_y_reverse(limits = c(100,0)) +
    labs(x = 'Frequency (hours)',
         y = '|Accuracy|',
         caption = '15 minute NO3 data from HBEF W3 2016 WY resampled by every nth measurement, compared to truth using every sample and the composite method.
         \n Lines indicate hourly, daily, weekly, biweekly, monthly, and bimonthly intervals.')+
    theme_classic()+
    theme(text = element_text(size = 20))+
    geom_vline(xintercept = 1)+ #hourly
    geom_vline(xintercept = 24)+ #daily
    geom_vline(xintercept = 96)+ #weekly
    geom_vline(xintercept = 192)+ #biweekly
    geom_vline(xintercept = 384)+ #monthly
    geom_vline(xintercept = 768) #bimonthly






