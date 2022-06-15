library(sf)
library(dataRetrieval)
library(tidyverse)
library(feather)
library(glue)
library(macrosheds)
library(ggplot2)
library(ggthemes)

source('source/helper_functions.R')

## thinning_intervals <- c('daily', 'weekly', 'biweekly', 'monthly', 'quarterly')
thinning_intervals <- c('daily', 'weekly', 'biweekly', 'monthly')
vars <- c('nitrate_nitrite_mgl')

# Read in site data
# (Sites were identified with the identify_usgs_gauges.R script)
site_var_data <- read.csv('data/general/site_var_data.csv', colClasses = 'character') %>%
    filter(parm_cd %in% c('99133')) %>%
    distinct(site_code, parm_cd, .keep_all = T)

n_sites <- site_var_data %>%
    filter(parm_cd == '99133')

# Only want sites with continuous N
site_var_data <- site_var_data %>%
    filter(site_code %in% !!n_sites$site_code)

write.csv(site_var_data, "nitrate_usgs_sites.csv")

# list USGS sites w nitrate
sites <- unique(site_var_data$site_code)

# read in USGS site data
# watershed info
ws <- read.csv("data/general/site_data.csv", colClasses = "character") %>%
    filter(site_code %in% !!sites) %>%
    sf::st_as_sf(., coords = c('long', 'lat'), crs = 4326)

# ecoregion
## eco <- read.csv("data/general/site_eco.csv", colClasses = "character") %>%
##   rename(ecoregion = NA_L2NAME) %>%
##   select(site_code, ecoregion)

## site_codes <- c()
## for(row in 1:nrow(eco)) {
##   row_code <- eco[i,]$site_code
##   if(nchar(row_code) != 8) {
##     code <- paste0("0", row_code)
##     site_codes <- c(site_codes, code)
##   } else {
##     site_codes <- c(site_codes, row_code)
##   }
## }

## eco$site_code <- site_codes

# make a feature collection of sites and their ecoregions
epa_eco_ii <- read_sf("data/spatial/NA_CEC_Eco_Level2.shp") %>%
  select(NA_L2NAME, geometry) %>%
  st_transform(crs = 4326) %>%
  st_make_valid(epa_eco_ii)

site_eco <- st_intersection(ws, epa_eco_ii)

# flashiness (RBI)
rbi <- data.frame()

for(i in 1:nrow(site_eco)) {

    site <- site_eco$site_code[i]
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site}.feather')))

    if(inherits(q_data, 'try-error')) {
      print("no Q data for site")
      output <- c(site, NA)
      rbi <- rbind(rbi, output)
      next
    }

    site_RBI <- ContDataQC::RBIcalc(q_data$val)
    output <- c(site, site_RBI)
    rbi <- rbind(rbi, output)
}


colnames(rbi) <- c("site_code", "RBI")

# merged
usgs <- site_eco %>%
  left_join(rbi, by = c('site_code'))

usgs <- usgs %>%
  arrange(desc(RBI))

usgs$dataset <- "USGS"

usgs.df <- usgs %>%
      mutate(long = unlist(map(usgs$geometry,1)),
             lat = unlist(map(usgs$geometry,2))) %>%
  st_drop_geometry()

write.csv(usgs.df,
         "streamlined/data/site/usgs_nitrate_sites.csv")



ms_sites <- ms_download_site_data()
ms.df <- ms_sites %>%
  mutate(dataset = "macrosheds") %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  select(site_code, ws_area_ha, dataset, latitude, longitude)

ms_ws <- ms.df %>%
    sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)


ms_eco <- st_intersection(ms_ws, epa_eco_ii)

# get MS nitrate sites
ms_cat <- ms_catalog()

# Watershed Attribute Density Plots: MacroSheds and USGS
# MacroSheds
# merge ws area df
usgs.ha <- usgs %>%
  as.data.frame() %>%
  select(site_code, ws_area_ha, dataset) %>%
  mutate(ws_area_ha = as.numeric(ws_area_ha))

ws.area <- usgs.ha %>%
  rbind(ms.ha)

# density plots

# watershed area USGS vs MS
p <- ggplot(ws.area, aes(x=log(ws_area_ha), color = dataset)) +
  geom_density() +
  theme_few()


# flashiness
