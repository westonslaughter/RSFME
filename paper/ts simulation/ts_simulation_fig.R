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
# loop start #####
loop_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character(), runid = as.integer())
for(i in 1:50){
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


dn %>%
    select(datetime, IS_discharge, IS_NO3, IS_FDOM) %>%
    pivot_longer(cols = -datetime, values_to = 'val', names_to = 'var') %>%
    ggplot(., aes(x = datetime, y = val)) +
    #geom_point()+
    geom_line(lwd = 2)+
    facet_wrap(vars(var), ncol = 1, scales = 'free')+
    theme_classic()+
    theme(text = element_text(size = 20),
          axis.title = element_blank())

# plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
# plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)
# plot(log10(dn$IS_discharge), log10(dn$IS_NO3), type = 'p', lwd = 2)


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
hold_factor <- (sum_q/sum(simulated_series[[3]]))
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
# make random sampling function
rtnorm <- function(n, mean, sd, a = 0, b = max(dn$IS_NO3)){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}
#### apply to make chemo #####
simulated_series[[4]] <- rtnorm(n = nrow(dn), sd = (sd(dn$IS_NO3)/4), mean = mean(dn$IS_NO3))
lines(dn$datetime, simulated_series[[4]], col = 'blue', type = 'l')

### no pattern ####
# apply to make no pattern ts
simulated_series[[5]] <- rtnorm(n = nrow(dn), sd = sd(dn$IS_NO3), mean = mean(dn$IS_NO3))
lines(dn$datetime, simulated_series[[4]], col = 'blue', type = 'l')

# plot(simulated_series[[5]]~dn$IS_discharge, data = dn)

#### enriching ####
##### fit lm to sp ts ####
fit_fdom <- lm(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)
inter_range <- runif(1000, min = confint(fit_fdom)[1,1], max = confint(fit_fdom)[1,2])
coef_range <- runif(1000, min = confint(fit_fdom)[2,1], max = confint(fit_fdom)[2,2])
error_range <- rnorm(1000, mean = 0, sd = sd(dn$IS_FDOM)/2)
simulated_series[[6]] <- as.numeric()
##### for all #####
for(j in 1:length(simulated_series[[1]])){
    inter <- sample(inter_range, size = 1)
    slope <- sample(coef_range, size = 1)
    error <- sample(error_range, size = 1)
    q <- simulated_series[[1]][j]

    pre_error_val <- 10^((log10(q)*slope)+inter)

    eps <- pre_error_val+(error*(pre_error_val)/mean(dn$IS_FDOM))

    simulated_series[[6]][j] <- pre_error_val + eps

}
#plot(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)

# check c:q
#plot(log10(simulated_series[[6]])~log10(simulated_series[[1]]))
#summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[1]])))
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
loop_out <- run_out %>%
        mutate(runid = i) %>%
        rbind(., loop_out)
}

# Figure creation #####
### make header plots #####
# unaltered q plot
p1 <- ggplot(dn, aes(x = date))+
        geom_line(aes(y = simulated_series[[1]])) +
        theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    labs(title = 'Unaltered Flow')+
    scale_y_log10(limits = c(1,1e3))

p1

# storm q plot
p2 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[2]])) +
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          text = element_text(size = 20))+
    labs(title = 'Stormflow Dominated',
         y = 'Q (lps)')+
    scale_y_log10(limits = c(1,1e3))
p2

# base q plot
p3 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[3]])) +
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 20))+
    labs(title = 'Baseflow Dominated')+
    scale_y_log10(limits = c(1,1e3))


p3

# chemo cq
p4 <- tibble(q = simulated_series[[1]], con = simulated_series[[4]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10() +
    labs(title = 'Chemostatic',
         y = 'C')+
    theme(axis.title.x=element_blank(),
          text = element_text(size = 20))

p4

# no pattern cq
p5 <- tibble(q = simulated_series[[1]], con = simulated_series[[5]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10()+
    labs(title = 'No Pattern',
         x = 'Q')+
    theme(axis.title.y=element_blank(),
          text = element_text(size = 20))
p5

# enrich cq
p6 <- tibble(q = simulated_series[[1]], con = simulated_series[[6]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10()+
    labs(title = 'Enriching',
         x = 'Q')+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          text = element_text(size = 20))
p6

### make row 1 plots ####
ymin = 0
ymax = .85

ymin_en = 300
ymax_en = 600
# unaltered flow w/ chemo data
p7_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'chemostatic') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p7_data$method <- factor(p7_data$method, levels = c("pw", "beale", "rating", 'composite'))

p7 <- ggplot(p7_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'unaltered' & loop_out$cq == 'chemostatic'])+
        geom_boxplot()+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin, ymax)+

p7

# unaltered flow w/ no pattern data
p8_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'none') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p8_data$method <- factor(p8_data$method, levels = c("pw", "beale", "rating", 'composite'))

p8 <- ggplot(p8_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'unaltered' & loop_out$cq == 'none'])+
    geom_boxplot()+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin, ymax)
p8

# unaltered flow w/ enrich data
p9_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'enrich') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p9_data$method <- factor(p9_data$method, levels = c("pw", "beale", "rating", 'composite'))

p9 <- ggplot(p9_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'unaltered' & loop_out$cq == 'enrich'])+
    geom_boxplot()+
   # ylim(ymin, ymax)+
    geom_hline(aes(yintercept=0))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.title.y = element_blank(),
          text = element_text(size = 20))+
    ylim(ymin_en, ymax_en)
p9

### make row 2 plots ####
# storm flow w/ chemo data
p10_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'chemostatic') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p10_data$method <- factor(p10_data$method, levels = c("pw", "beale", "rating", 'composite'))

p10 <- ggplot(p10_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'storm' & loop_out$cq == 'chemostatic'])+
    geom_boxplot()+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin,ymax)+
    labs(y = 'Load (kg/ha/yr)')
p10

# storm flow w/ no pattern data
p11_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'none') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p11_data$method <- factor(p11_data$method, levels = c("pw", "beale", "rating", 'composite'))

p11 <- ggplot(p11_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'storm' & loop_out$cq == 'none'])+
    geom_boxplot()+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin, ymax)

p11

# stormflow w/ enrich data
p12_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'enrich') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p12_data$method <- factor(p12_data$method, levels = c("pw", "beale", "rating", 'composite'))

p12 <- ggplot(p12_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'storm' & loop_out$cq == 'enrich'])+
    geom_boxplot()+
    geom_hline(aes(yintercept=0))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin_en, ymax_en)

p12

### make row 3 plots ####
# base flow w/ chemo data
p13_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'chemostatic') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p13_data$method <- factor(p13_data$method, levels = c("pw", "beale", "rating", 'composite'))

p13 <- ggplot(p13_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'base' & loop_out$cq == 'chemostatic'])+
    geom_boxplot()+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin, ymax)+
    geom_hline(aes(yintercept=0))
p13

# base flow w/ no pattern data
p14_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'none') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p14_data$method <- factor(p14_data$method, levels = c("pw", "beale", "rating", 'composite'))

p14 <- ggplot(p14_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'base' & loop_out$cq == 'none'])+
    geom_boxplot()+
    theme_classic() +
    ylim(ymin, ymax)+
    theme(axis.title.y=element_blank(),
          text = element_text(size = 20)) +
    labs(x = 'Method')
p14

# base flow w/ enrich data
p15_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'enrich') %>%
    pivot_wider(names_from = method, values_from = estimate, id_cols = runid, values_fn = mean) %>%
    # mutate(pw = ((pw-truth)/truth)*100,
    #        beale = ((beale - truth)/truth)*100,
    #        rating = ((rating-truth)/truth)*100,
    #        composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = everything() ,names_to = 'method', values_to = 'error')

p15_data$method <- factor(p15_data$method, levels = c("pw", "beale", "rating", 'composite'))

p15 <- ggplot(p15_data, aes(x = method, y = error)) +
    geom_hline(yintercept = loop_out$estimate[loop_out$method == 'truth' & loop_out$flow == 'base' & loop_out$cq == 'enrich'])+
    geom_boxplot()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    ylim(ymin_en, ymax_en)

p15

## make mega fig ####
(plot_spacer() | p4 | p5 | p6)/
(p1 | p7 | p8 | p9)/
(p2 | p10 | p11 | p12)/
(p3 | p13 | p14 | p15)

