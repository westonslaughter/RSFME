
# pick site
site_no <-  '03346500'
area <- 604247

# grab 'truth'
true <- read_feather(paste0('data/fluxes/daily/true/nitrate_nitrite_mgl/', site_no,'.feather')) %>%
    mutate(method = 'true',
           flux = flux/!!area) %>%
    select(date = datetime, flux, method, wy)
    
# aggregate all data from biweekly thinned data
beale <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/beale/', site_no,'.feather')) %>%
    rename(date = datetime) %>%
    mutate(method = 'beale',
           flux = flux/1000) %>%
    select(date, flux, method, wy)

# bear should be fixed now
bear <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/bear/', site_no,'.feather')) %>%
    mutate(flux = flux_kg_ha*60*15) %>%
    select(date, flux, method,wy)

egret <- readRDS(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/egret/', site_no,'.rds')) %>%
    .$Daily %>%
    mutate(flux = FluxDay/(!!area*100), method = 'wrtds',
           wy = water_year(Date, origin = 'usgs')) %>%
    select(date = Date, flux, method, wy)
    
fernow <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/fernow/', site_no,'.feather')) %>%
    mutate(flux = flux_weekly_kg_ha/100) %>%
    select(date, flux, method, wy)

hbef <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/hbef/', site_no,'.feather')) %>%
    mutate(flux = flux_daily_kg_ha/100) %>%
    select(date, flux, method, wy)

rating <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/rating/', site_no,'.feather')) %>%
    mutate(flux = flux_kg_ha/100) %>%
    select(date, flux, method, wy)

santee <- read_feather(paste0('data/fluxes/daily/biweekly/nitrate_nitrite_mgl/santee/', site_no,'.feather')) %>%
    mutate(flux = flux_kg_ha/100) %>%
    select(date, flux, method, wy)

# combine all data into a single frame and filter to one water year
dat <- rbind(true, beale, bear, egret, fernow, hbef, rating, santee) %>%
    filter(wy == 2016)

ggplot(dat, aes(x = date, y = flux, color = method))+
    geom_point()+
    scale_y_log10()


### main loop tests
daily_egret_flux <- egret_flux$Daily%>%
    mutate(flux = FluxDay/(!!area), method = 'wrtds',
           wy = water_year(Date, origin = 'usgs')) %>%
    select(date = Date, flux, method, wy)
test <- rbind(true_flux_daily, daily_beale_flux, bear_flux_daily, hbef_flux_daily, daily_santee_flux, daily_rating_flux, daily_egret_flux)
ggplot(test, aes(x = date, y = flux, color = method)) + geom_point() #+ scale_y_log10()
