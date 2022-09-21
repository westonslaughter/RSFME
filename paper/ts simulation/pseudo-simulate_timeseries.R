library(tidyverse)
library(forecast)
library(feather)
library(xts)
library(imputeTS)
library(here)

ts_len = 1000
d = read_feather(here('data/ms/hbef/discharge/w1.feather')) %>%
    slice(1:ts_len)

plot(d$datetime, d$val, type = 'l', lwd = 2)

#simulation approach 1: fit ARIMA model to series; resample the residuals ####

fit <- auto.arima(xts(d$val, order.by = d$datetime))
# methods(class='Arima')

lines(d$datetime, fit$fitted, col = 'blue', lwd = 2)

simulated_series = list()
for(i in seq_len(10)){
    Sys.sleep(1)
    resampled_residuals = sample(fit$residuals, size = ts_len, replace = TRUE)
    simulated_series[[i]] = d$val + resampled_residuals
    lines(d$datetime, simulated_series[[i]], col = i)
}

#corollary to approach 1: mess with stuff ####

resampled_residuals = sample(fit$residuals, size = ts_len, replace = TRUE)

#scale
simulated_series[[11]] = d$val * 5 + resampled_residuals
plot(d$datetime, simulated_series[[11]], col = 'red', type = 'l')
simulated_series[[12]] = d$val * 0.1 + resampled_residuals
plot(d$datetime, simulated_series[[12]], col = 'green', type = 'l')

#trend
simulated_series[[13]] = d$val + seq(1, 10, length.out = ts_len) + resampled_residuals
plot(d$datetime, simulated_series[[13]], col = 'blue', type = 'l')
simulated_series[[14]] = d$val + sort(resampled_residuals)
plot(d$datetime, simulated_series[[14]], col = 'brown', type = 'l')

#regime
series_subset_length = 100
d2 = select(d, datetime, val) %>% slice(1:series_subset_length)
d2$datetime = d$datetime[seq(1, ts_len, length.out = series_subset_length)]
simulated_series[[15]] = d2 %>%
    complete(datetime = seq(min(datetime), max(datetime), 'day')) %>%
    mutate(val = imputeTS::na_seadec(val)) %>%
    pull(val)
lines(d$datetime, simulated_series[[15]], col = 'orange', type = 'l')


