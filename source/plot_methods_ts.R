library(tidyverse)
library(feather)
library(glue)
library(ggthemes)
library(viridis)


read_add_site <- function(flux_f){

    site_code <- str_split_fixed(flux_f, '/', n = Inf)[,7]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    fluxes <- read_feather(flux_f) %>%
        mutate(site_code = !!site_code) %>%
        mutate(wy = as.factor(wy),
               flux_annual_kg_ha = as.numeric(flux_annual_kg_ha),
               method = as.character(method),
               site_code = as.character(site_code))

    return(fluxes)
}

#### Nitrate Nitrite as N ####
## True fluxes
real_fluxes_f <- list.files('data/fluxes/annual/true/nitrate_nitrite_mgl', full.names = TRUE)

real_fluxes <- tibble()
for(i in 1:length(real_fluxes_f)){
    site_code <- str_split_fixed(real_fluxes_f[i], '/', n = Inf)[,6]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    real_fluxes_this <- read_feather(real_fluxes_f[i]) %>%
        mutate(site_code = !!site_code)

    real_fluxes <- rbind(real_fluxes, real_fluxes_this)
}


sites <- unique(real_fluxes$site_code)
site_data <- read_csv('data/general/site_data.csv') %>%
    filter(site_code %in% !!sites) %>%
    sf::st_as_sf(., coords = c('long', 'lat'), crs = 4326)


## Daily
daily_fluxes_f <- list.files('data/fluxes/annual/daily/nitrate_nitrite_mgl',
                              full.names = TRUE,
                              recursive = TRUE)

daily_fluxes_f <- grep('[.]feather', daily_fluxes_f, value = TRUE)

daily_fluxes <- map_dfr(daily_fluxes_f, read_add_site) %>%
    mutate(thinning = 'daily')

## Weekely
weekly_fluxes_f <- list.files('data/fluxes/annual/weekly/nitrate_nitrite_mgl',
                              full.names = TRUE,
                              recursive = TRUE)

weekly_fluxes_f <- grep('[.]feather', weekly_fluxes_f, value = TRUE)

weekly_fluxes <- map_dfr(weekly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'weekly')


## Biweekly
biweekly_fluxes_f <- list.files('data/fluxes/annual/biweekly/nitrate_nitrite_mgl',
                              full.names = TRUE,
                              recursive = TRUE)

biweekly_fluxes_f <- grep('[.]feather', biweekly_fluxes_f, value = TRUE)

biweekly_fluxes <- map_dfr(biweekly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'biweekly')

## Monthly
monthly_fluxes_f <- list.files('data/fluxes/annual/monthly/nitrate_nitrite_mgl',
                                full.names = TRUE,
                                recursive = TRUE)

monthly_fluxes_f <- grep('[.]feather', monthly_fluxes_f, value = TRUE)

monthly_fluxes <- map_dfr(monthly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'monthly')

## Combine data
real_fluxes <- real_fluxes %>%
    rename(real = flux_annual_kg_ha) %>%
    select(-method)

comp <- rbind(daily_fluxes, weekly_fluxes, biweekly_fluxes, monthly_fluxes) %>%
    merge(real_fluxes, by = c('wy', 'site_code')) %>%
    mutate(dif = real-flux_annual_kg_ha) %>%
    mutate(percent_dif = (dif/((real+flux_annual_kg_ha)/2))*100)


comp$thinning <- factor(comp$thinning , levels=c('daily', 'weekly', 'biweekly', 'monthly'))

# plotting globals
pd <- position_dodge(0.4)

# Time Series of Annual Flux Estimates % Differences from Real Value by Thinning Method
n_all <- ggplot(comp, aes(x = wy, y = percent_dif)) +
    stat_summary(aes(group = thinning, color = thinning),
               fun = mean, na.rm = TRUE,
               geom = "line",
               size = .75, alpha = .2,
               position = pd
               ) +
    stat_summary(aes(group = thinning, color = thinning),
                 fun = mean,
                 fun.min = function(x) min(x),
                 fun.max = function(x) max(x),
                 na.rm = TRUE,
                 geom = "pointrange",
                 position = pd
               ) +
    geom_hline(yintercept = 0, linetype = 'solid', color = 'gray', alpha = .6) +
    facet_wrap(~method, ncol  = 1) +
    scale_color_viridis(discrete=TRUE, option="turbo") +
    theme_few() +
    theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank()) +
    labs(x = 'Water Year', y = 'Percent Dif. from True Nitrate Flux') +
    ggtitle("Annual Flux Estimate % Difference from True Nitrate Flux by Calculation Method \nFrom Daily Data 2011-2021")

ggsave(plot = n_all, filename = 'plots/ts/ts_n_all_methods.png', height = 8)
