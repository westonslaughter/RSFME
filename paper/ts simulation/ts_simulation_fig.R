library(tidyverse)
library(forecast)
library(feather)
library(xts)
library(imputeTS)
library(here)
library(lfstat)
library(lubridate)
library(ggpubr)

set.seed(53045)


source(here('source/flux_methods.R'))

area <- 42.4
site_code = 'w3'
# loop start #####
loop_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
for(i in 1:10){
#ts_len = 1000
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

mean_q <- mean(dn$IS_discharge)
sum_q <- sum(dn$IS_discharge)

max(dn$datetime)
min(dn$datetime)

max(dn$IS_NO3, na.rm = T)
min(dn$IS_NO3, na.rm = T)
sd(dn$IS_NO3, na.rm = T)
mean(dn$IS_NO3, na.rm = T)
hist(dn$IS_NO3)

summary(lm(log10(dn$IS_NO3)~log10(dn$IS_discharge), data = dn))

summary(lm(log10(dn$IS_spCond)~log10(dn$IS_discharge), data = dn))

summary(lm(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn))

# plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
# plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)
# plot(log10(dn$IS_discharge), log10(dn$IS_spCond), type = 'p', lwd = 2)


# DISCHARGE TS #####
plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)

## fit ARIMA model to series; resample the residuals ####

fit <- auto.arima(xts(dn$IS_discharge, order.by = dn$datetime))
# methods(class='Arima')

lines(dn$datetime, fit$fitted, col = 'blue', lwd = 2)

#ts_len = 1000
resampled_residuals = sample(fit$residuals,
                             #size = ts_len,
                             replace = TRUE)

## create random models
simulated_series = list()
# for(i in seq_len(10)){
#     Sys.sleep(1)
#     resampled_residuals = sample(fit$residuals, size = ts_len, replace = TRUE)
#     simulated_series[[i]] = d$val + resampled_residuals
#     lines(d$datetime, simulated_series[[i]], col = i)
# }

## make TS#####
###unaltered #####
reg <- dn$IS_discharge
reg[dn$IS_discharge < 1] = 1
simulated_series[[1]] = reg + resampled_residuals + 5
simulated_series[[1]][which(simulated_series[[1]] <= 0)] = 1
hold_factor <- (sum_q/sum(simulated_series[[1]]))
simulated_series[[1]] <- simulated_series[[1]]*hold_factor
lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')

###stormflow dominated ####
storm <- dn$IS_discharge
storm[dn$IS_discharge < 1] = 1
simulated_series[[2]] = storm^1.5 + resampled_residuals + 5
simulated_series[[2]][which(simulated_series[[2]] <= 0)] = 1
hold_factor <- (sum_q/sum(simulated_series[[2]]))
simulated_series[[2]] <- simulated_series[[2]]*hold_factor
lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')

###baseflow dominated ####
base <- dn$IS_discharge
base[dn$IS_discharge < 1] = 1
simulated_series[[3]] = storm^0.8 + resampled_residuals + 5
simulated_series[[3]][which(simulated_series[[3]] <= 0)] = 1
hold_factor <- (sum/sum(simulated_series[[3]]))
simulated_series[[3]] <- simulated_series[[3]]*hold_factor
lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')


# plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
# lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')
# lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')
# lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')

# CON TS ####

## make TS#####
# plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)

### chemostatic #####
##### fit ARIMA model to series; resample the residuals ####

fit_n <- auto.arima(xts(dn$IS_NO3, order.by = dn$datetime))

resampled_residuals_n = sample(fit_n$residuals,
                               #size = ts_len,
                               replace = TRUE)

#### fit chemostatic #####
chemo_base <- dn$IS_NO3
simulated_series[[4]] = chemo_base+ resampled_residuals_n
lines(dn$datetime, simulated_series[[4]], col = 'red', type = 'l')

### no pattern ####
# make random sampling function
rtnorm <- function(n, mean, sd, a = 0, b = max(dn$IS_NO3)){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

# apply to make no pattern ts
simulated_series[[5]] <- rtnorm(n = nrow(dn), sd = sd(chemo_base), mean = mean(chemo_base))
lines(dn$datetime, simulated_series[[4]], col = 'blue', type = 'l')

# plot(simulated_series[[5]]~dn$IS_discharge, data = dn)

#### enriching ####
##### fit lm to sp ts ####
fit_fdom <- lm(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)
inter_range <- runif(1000, min = 10^confint(fit_fdom)[1,1], max = 10^confint(fit_fdom)[1,2])
coef_range <- runif(1000, min = 10^confint(fit_fdom)[2,1], max = 10^confint(fit_fdom)[2,2])
error_range <- rnorm(1000, mean = 0, sd = sd(dn$IS_FDOM))
simulated_series[[6]] <- as.numeric()
##### for all #####
for(i in 1:length(simulated_series[[1]])){
    inter <- sample(inter_range, size = 1)
    slope <- sample(coef_range, size = 1)
    error <- sample(error_range, size = 1)
    q <- simulated_series[[1]][i]

    pre_error_val <- (q*slope)+inter

    eps <- pre_error_val*(error/100)

    simulated_series[[6]][i] <- pre_error_val + eps

}
#plot(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)

# check c:q
# plot(log10(simulated_series[[6]])~log10(simulated_series[[1]]))
# summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[1]])))
#
# plot(log10(simulated_series[[6]])~log10(simulated_series[[2]]))
# summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[2]])))
#
# plot(log10(simulated_series[[6]])~log10(simulated_series[[3]]))
# summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[3]])))

# ESTIMATE FLUX#####
# coarsen function
coarsen_data <- function(chem_df){
out <- chem_df %>%
    filter(hour(datetime) %in% c(13:18)) %>%
    filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
    mutate(date = lubridate::date(datetime)) %>%
    distinct(date, .keep_all = T)
return(out)
}

make_q_daily <- function(q_df){
out <- q_df %>%
    group_by(lubridate::yday(datetime)) %>%
    summarize(date = date(datetime),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    select(date, q_lps)
}

calculate_truth <- function(raw_chem_list, q_df, flow_regime = NULL, cq = NULL){
    chem_df <- tibble(datetime = dn$datetime, con = raw_chem_list) %>%
        group_by(lubridate::yday(datetime)) %>%
        summarize(date = date(datetime),
                  con = mean(con)) %>%
        ungroup() %>%
        unique() %>%
        select(date, con) %>%
        mutate(site_code = 'w3', wy = target_wy)

    q_df_add <- q_df %>%
        mutate(site_code = 'w3', wy = target_wy)

    out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df_add, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)
    out <- tibble(method = 'truth', estimate = out_val,
                  flow = flow_regime, cq = cq)

    return(out)
}

apply_methods <- function(chem_df, q_df, flow_regime = NULL, cq = NULL){
    out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
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
    out$flow <- flow_regime
    out$cq <- cq
    return(out)

}
run_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
### chemostatic ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[4]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)
##### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)

##### under storm domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, flow_regime = 'storm', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'storm', cq = 'chemostatic'),
    run_out)

##### under base domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, flow_regime = 'base', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'base', cq = 'chemostatic'),
    run_out)

### no pattern ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[5]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

##### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[5]], q_df, flow_regime = 'unaltered', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'unaltered', cq = 'none'),
    run_out)

##### under storm domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[5]], q_df, flow_regime = 'storm', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'storm', cq = 'none'),
    run_out)

##### under base domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, flow_regime = 'base', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'base', cq = 'none'),
    run_out)

### strong enrich ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[6]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

##### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], q_df, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)

##### under storm domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], q_df, flow_regime = 'storm', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'storm', cq = 'enrich'),
    run_out)

##### under base domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], q_df, flow_regime = 'base', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, flow_regime = 'base', cq = 'enrich'),
    run_out)

### save out ####
loop_out <- rbind(run_out, loop_out)
}

p1 <- ggplot(dn, aes(x = date))+
        geom_line(aes(y = simulated_series[[1]])) +
        ylim(0, 1500)+
        theme_classic()
p1

p2 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[2]])) +
    ylim(0, 1500)+
    theme_classic()
p2

p3 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[3]])) +
    ylim(0, 1500)+
    theme_classic()

p3




