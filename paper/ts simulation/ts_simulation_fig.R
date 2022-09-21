library(tidyverse)
library(forecast)
library(feather)
library(xts)
library(imputeTS)
library(here)
library(lfstat)

set.seed(53045)

#ts_len = 1000
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))
    #slice(1:ts_len)

# subset to 2015 wy
dn <- d %>%
    filter(wy == 2016) %>%
    mutate(IS_discharge = na.approx(IS_discharge),
           IS_NO3 = na.approx(IS_NO3),
           IS_FDOM = na.approx(IS_FDOM),
           IS_spCond = na.approx(IS_spCond))

mean_q <- mean(dn$IS_discharge)

max(dn$datetime)
min(dn$datetime)

max(dn$IS_NO3, na.rm = T)
min(dn$IS_NO3, na.rm = T)
sd(dn$IS_NO3, na.rm = T)
mean(dn$IS_NO3, na.rm = T)
hist(dn$IS_NO3)

summary(lm(log(dn$IS_NO3)~log(dn$IS_discharge), data = dn))

summary(lm(log(dn$IS_spCond)~log(dn$IS_discharge), data = dn))

summary(lm(log(dn$IS_FDOM)~log(dn$IS_discharge), data = dn))

plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)
plot(dn$IS_discharge, dn$IS_FDOM, type = 'p', lwd = 2)


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
hold_factor <- (mean_q/mean(simulated_series[[1]]))
simulated_series[[1]] <- simulated_series[[1]]*hold_factor
lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')

###stormflow dominated ####
storm <- dn$IS_discharge
storm[dn$IS_discharge < 1] = 1
simulated_series[[2]] = storm^1.5 + resampled_residuals + 5
simulated_series[[2]][which(simulated_series[[2]] <= 0)] = 1
hold_factor <- (mean_q/mean(simulated_series[[2]]))
simulated_series[[2]] <- simulated_series[[2]]*hold_factor
lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')

###baseflow dominated ####
base <- dn$IS_discharge
base[dn$IS_discharge < 1] = 1
simulated_series[[3]] = storm^0.8 + resampled_residuals + 5
simulated_series[[3]][which(simulated_series[[3]] <= 0)] = 1
hold_factor <- (mean_q/mean(simulated_series[[3]]))
simulated_series[[3]] <- simulated_series[[3]]*hold_factor
lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')


plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')
lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')
lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')

# CON TS ####

## make TS#####
plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)

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

plot(simulated_series[[5]]~dn$IS_discharge, data = dn)

#### diluting ####
##### fit ARIMA model to series; resample the residuals ####
fit_sc <- auto.arima(xts(dn$IS_spCond, order.by = dn$datetime))

resampled_residuals_sc = sample(fit_n$residuals,
                               #size = ts_len,
                               replace = TRUE)
## for unaltered #####
summary(lm(log(dn$IS_FDOM)~log(dn$IS_discharge), data = dn))

enrich_base <- dn$IS_spCond
simulated_series[[6]] <- enrich_base + resampled_residuals_sc

# check c:q
plot(log(simulated_series[[6]])~log(simulated_series[[1]]))
summary(lm(log(simulated_series[[6]])~log(simulated_series[[1]])))

plot(log(simulated_series[[6]])~log(simulated_series[[2]]))
summary(lm(log(simulated_series[[6]])~log(simulated_series[[2]])))

plot(log(simulated_series[[6]])~log(simulated_series[[3]]))
summary(lm(log(simulated_series[[6]])~log(simulated_series[[3]])))

# estimate flux
