library(sf)
library(dataRetrieval)
library(dplyr)

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

# Site Density PLots
# bring in hand-fluxed data
## flux_f <- daily_fluxes_f[1]
read_add_site <- function(flux_f){

    site_code <- str_split_fixed(flux_f, '/', n = Inf)[,7]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    fluxes <- read_feather(flux_f) %>%
      mutate(site_code = !!site_code) %>%
        mutate(wy = as.factor(wy),
               flux_annual_kg_ha = as.numeric(flux),
               ## method = as.character(method),
               site_code = as.character(site_code))

    return(fluxes)
}

real_fluxes_f <- list.files('data/tmp/fluxes/annual/true/nitrate_nitrite_mgl', full.names = TRUE)

real_fluxes <- tibble()
for(i in 1:length(real_fluxes_f)){
    site_code <- str_split_fixed(real_fluxes_f[i], '/', n = Inf)[,7]
    site_code <- str_split_fixed(site_code, '[.]', n = Inf)[,1]
    real_fluxes_this <- read_feather(real_fluxes_f[i]) %>%
        mutate(site_code = !!site_code)

    real_fluxes <- rbind(real_fluxes, real_fluxes_this)
}


sites <- unique(real_fluxes$site_code)

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

epa_eco <- read_sf("data/spatial/NA_CEC_Eco_Level2.shp")
epa_eco_ii <- epa_eco %>%
  select(NA_L2NAME, geometry) %>%
  st_transform(crs = 4326)

# if needed
## sf_use_s2(FALSE)
epa_eco_ii <- st_make_valid(epa_eco_ii)
ws_sf <- st_make_valid(ws)

# make a feature collection of sites and their ecoregions
site_eco <- st_intersection(ws, epa_eco_ii)

# flashiness
## rbi <- read.csv("data/general/site_RBI_data.csv")
site_RBI_df <- data.frame()

for(i in 1:nrow(df)) {

    site <- df$site_code[i]
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site}.feather')))

    if(inherits(q_data, 'try-error')) next

    ## parm_cd <- df$parm_cd[i]
    ## var <- variable_data %>%
    ##     filter(usgs_parm_cd == !!parm_cd) %>%
    ##   pull(var)

    site_RBI <- ContDataQC::RBIcalc(q_data$val)
    output <- c(site, site_RBI)
    site_RBI_df <- rbind(site_RBI_df, output)
}

rbi <- site_RBI_df

colnames(rbi) <- c("site_code", "RBI")

## site_codes <- c()
## for(row in 1:nrow(rbi)) {
##   row_code <- rbi[i,]$site_code
##   if(nchar(row_code) != 8) {
##     code <- paste0("0", row_code)
##     site_codes <- c(site_codes, code)
##   } else {
##     site_codes <- c(site_codes, row_code)
##   }
## }

## rbi$site_code <- site_code

# merged
usgs <- site_eco %>%
  ## left_join(eco, by = c('site_code')) %>%
  left_join(rbi, by = c('site_code'))
usgs$dataset <- "USGS"

usgs.ha <- usgs %>%
  select(site_code, ws_area_ha, dataset) %>%
  mutate(ws_area_ha = as.numeric(ws_area_ha)) %>%
  as.data.frame()

# MacroSheds
library(macrosheds)

ms_sites <- ms_download_site_data()
ms.ha <- ms_sites %>%
  mutate(dataset = "macrosheds") %>%
  select(site_code, ws_area_ha, dataset)

# merge ws area df
ws.area <- usgs.ha %>%
  rbind(ms.ha)

# density plots
p <- ggplot(ws.area, aes(x=log(ws_area_ha), color = dataset)) +
  geom_density() +
  theme_few()

# watershed area


# flashiness
