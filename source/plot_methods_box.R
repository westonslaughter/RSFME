library(tidyverse)
library(feather)
library(glue)
library(ggthemes)


read_add_site <- function(flux_f){
    
    site_code <- str_split_fixed(flux_f, '/', n = Inf)[,7]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    fluxes <- read_feather(flux_f) 
    if('flux_annual_kg_ha' %in% colnames(fluxes)){
       fluxes <- fluxes %>%
            rename(flux = flux_annual_kg_ha)
        }
    fluxes <- fluxes %>%
        mutate(site_code = !!site_code) %>%
        mutate(wy = as.factor(wy),
               flux = as.numeric(flux),
               method = as.character(method),
               site_code = as.character(site_code))
    
    return(fluxes)
}

read_add_site2 <- function(flux_f){
    
    site_code <- str_split_fixed(flux_f, '/', n = Inf)[,8]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    fluxes <- read_feather(flux_f) 
    if('flux_annual_kg_ha' %in% colnames(fluxes)){
        fluxes <- fluxes %>%
            rename(flux = flux_annual_kg_ha)
    }
    fluxes <- fluxes %>%
        mutate(site_code = !!site_code) %>%
        mutate(wy = as.factor(wy),
               flux = as.numeric(flux),
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

## Weekly 
weekly_fluxes_f <- list.files('data/fluxes/annual/weekly/nitrate_nitrite_mgl/', 
                              full.names = TRUE,
                              recursive = TRUE)

weekly_fluxes_f <- grep('[.]feather', weekly_fluxes_f, value = TRUE)

weekly_fluxes <- map_dfr(weekly_fluxes_f, read_add_site2) %>%
    mutate(thinning = 'weekly')


## Biweekly
biweekly_fluxes_f <- list.files('data/fluxes/annual/biweekly/nitrate_nitrite_mgl/', 
                              full.names = TRUE,
                              recursive = TRUE)

biweekly_fluxes_f <- grep('[.]feather', biweekly_fluxes_f, value = TRUE)

biweekly_fluxes <- map_dfr(biweekly_fluxes_f, read_add_site2) %>%
    mutate(thinning = 'biweekly')

## Monthly
monthly_fluxes_f <- list.files('data/fluxes/annual/monthly/nitrate_nitrite_mgl/', 
                                full.names = TRUE,
                                recursive = TRUE)

monthly_fluxes_f <- grep('[.]feather', monthly_fluxes_f, value = TRUE)

monthly_fluxes <- map_dfr(monthly_fluxes_f, read_add_site2) %>%
    mutate(thinning = 'monthly')

## Combine data
real_fluxes <- real_fluxes %>%
    rename(real = flux) %>%
    select(-method) %>%
    filter(real != 0)

comp <- full_join(rbind(daily_fluxes, weekly_fluxes, biweekly_fluxes,monthly_fluxes), real_fluxes, by = c('wy', 'site_code')) %>%
    mutate(dif = (real-flux)) %>%
    mutate(percent_dif = (dif/((real+flux)/2))*100)


test <- rbind(daily_fluxes, weekly_fluxes, biweekly_fluxes,monthly_fluxes) %>%
    filter(site_code == "01463500")

test2 <-  real_fluxes %>%
    filter(site_code == "01463500") 

test3 <- comp %>%
    filter(site_code == "01463500")
    
ggplot(test3, aes(x = wy, y = flux, color = method))+geom_point()+


comp$thinning <- factor(comp$thinning , levels=c('daily', 'weekly', 'biweekly', 'monthly'))

n_all <- ggplot(comp, aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    facet_wrap(~method, ncol  = 1) +
    theme_few() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_all, filename = 'plots/n_all_methods.png', height = 6)

# Bear
n_bear <- comp %>%
    filter(method == 'bear') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_bear, filename = 'plots/n_bear_methods.png', height = 6)

# Fernow
n_fernow <- comp %>%
    filter(method == 'fernow') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_fernow, filename = 'plots/n_fernow_methods.png', height = 6)

# hbef
n_hbef <- comp %>%
    filter(method == 'hbef') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_hbef, filename = 'plots/n_hbef_methods.png', height = 6)

# Santee
n_santee <- comp %>%
    filter(method == 'santee') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_santee, filename = 'plots/n_santee_methods.png', height = 6)

# wrtds
n_wrtds <- comp %>%
    filter(method == 'wrtds') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True Nitrate Flux')

ggsave(plot = n_wrtds, filename = 'plots/n_wrtds_methods.png', height = 6)

#### Spefic COnductivity ####
## True fluxes 
real_fluxes_f <- list.files('data/fluxes/annual/true/spcond_uscm/', full.names = TRUE)

real_fluxes <- tibble()
for(i in 1:length(real_fluxes_f)){
    site_code <- str_split_fixed(real_fluxes_f[i], '/', n = Inf)[,7]
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
daily_fluxes_f <- list.files('data/fluxes/annual/daily/spcond_uscm/', 
                             full.names = TRUE,
                             recursive = TRUE)

daily_fluxes_f <- grep('[.]feather', daily_fluxes_f, value = TRUE)

daily_fluxes <- map_dfr(daily_fluxes_f, read_add_site) %>%
    mutate(thinning = 'daily')

## Weekely 
weekly_fluxes_f <- list.files('data/fluxes/annual/weekly/spcond_uscm/', 
                              full.names = TRUE,
                              recursive = TRUE)

weekly_fluxes_f <- grep('[.]feather', weekly_fluxes_f, value = TRUE)

weekly_fluxes <- map_dfr(weekly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'weekly')


## Biweekly
biweekly_fluxes_f <- list.files('data/fluxes/annual/biweekly/spcond_uscm/', 
                                full.names = TRUE,
                                recursive = TRUE)

biweekly_fluxes_f <- grep('[.]feather', biweekly_fluxes_f, value = TRUE)

biweekly_fluxes <- map_dfr(biweekly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'biweekly')

## Monthly
monthly_fluxes_f <- list.files('data/fluxes/annual/monthly/spcond_uscm/', 
                               full.names = TRUE,
                               recursive = TRUE)

monthly_fluxes_f <- grep('[.]feather', monthly_fluxes_f, value = TRUE)

monthly_fluxes <- map_dfr(monthly_fluxes_f, read_add_site) %>%
    mutate(thinning = 'monthly')

## Combine data
real_fluxes <- real_fluxes %>%
    rename(real = flux) %>%
    select(-method)

comp <- full_join(rbind(daily_fluxes, weekly_fluxes, biweekly_fluxes,monthly_fluxes), real_fluxes) %>%
    mutate(dif = real-flux_annual_kg_ha) %>%
    mutate(percent_dif = (dif/((real+flux_annual_kg_ha)/2))*100) %>%
    filter(!is.na(method))


comp$thinning <- factor(comp$thinning , levels=c('daily', 'weekly', 'biweekly', 'monthly'))

spcond_all <- ggplot(comp, aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    facet_wrap(~method, ncol  = 1) +
    theme_few() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_all, filename = 'plots/spcond_all_methods.png', height = 6)

# Bear
spcond_bear <- comp %>%
    filter(method == 'bear') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_bear, filename = 'plots/spcond_bear_methods.png', height = 6)

# Fernow
spcond_fernow <- comp %>%
    filter(method == 'fernow') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_fernow, filename = 'plots/spcond_fernow_methods.png', height = 6)

# hbef
spcond_hbef <- comp %>%
    filter(method == 'hbef') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_hbef, filename = 'plots/spcond_hbef_methods.png', height = 6)

# Santee
spcond_santee <- comp %>%
    filter(method == 'santee') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_santee, filename = 'plots/spcond_santee_methods.png', height = 6)

# wrtds
spcond_wrtds <- comp %>%
    filter(method == 'wrtds') %>%
    ggplot(aes(x = thinning, y = percent_dif)) + 
    geom_boxplot() +
    labs(x = 'Thinning Interval', y = 'Percent Dif. from True spCond Flux')

ggsave(plot = spcond_wrtds, filename = 'plots/spcond_wrtds_methods.png', height = 6)



