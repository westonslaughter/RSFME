
library(ggplot2)
library(dplyr)
library(tidyr)
library(feather)
library(macrosheds)
library(RColorBrewer)
library(stringr)

source('source/helper_functions.R')
source('source/usgs_helpers.R')

# NOTE: seems flux calc sript is ocmppiling sites cumulatively, fix later

# HBEF Flux Estimates by Various Methods from MacroSheds RSFME Project
hbef_flux <- read_feather('data/ms/hbef/stream_flux/w9.feather') %>%
  mutate(var = ms_drop_var_prefix(var)) %>%
  select(-ms_recommended) %>%
  distinct(wy, site_code, method, var, .keep_all = TRUE) %>%
  pivot_wider(names_from = var, values_from = val,
              id_cols = c('wy', 'site_code', 'site_code', 'method'))

# Compile published flux from HBEF
hbef_pubs_dir <- file.path(getwd(), 'data/raw/hbef_published_flux/')
link_w1 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.3.17&entityid=520d38828fe2356314e51008a5059dd4'
link_w2 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.4.17&entityid=a6aeef15070be913ee2f06f431b9b7a7'
link_w3 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.5.17&entityid=82d579e0262732d4bc996890c0f4dbd3'
link_w4 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.6.17&entityid=54b3ae4a45a2bb6c7006c2ab45cf63b9'
link_w5 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.7.17&entityid=c08ebaccab4fee5fb60f4eee77f06cb3'
link_w6 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.8.17&entityid=3312389e77cc5fd06bc8a7c9019de0ed'
link_w7 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.9.18&entityid=11eb156e027c3af4e19ae48e335f35b2'
link_w8 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.10.18&entityid=f93eb6d324536491dd042c5496289dec'
link_w9 <- 'https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-hbr.11.17&entityid=bab6ac4dd3349bfd5cba711ecfd3d74f'
hbef_links <- c('w1'=link_w1,'w2'=link_w2,'w3'=link_w3,'w4'=link_w4,'w5'=link_w5,'w6'=link_w6,'w7'=link_w7,'w8'=link_w8,'w9'=link_w9)

for(ws in names(hbef_links)) {
  ws_filename <- paste0(ws, '.csv')
  ws_fp = file.path(hbef_pubs_dir, ws_filename)

  download.file(url = hbef_links[ws], ws_fp)

  ws_pub <- read.csv(ws_fp)
  ws_flux <- hbef_flux %>%
    filter(site_code == site)

  # solutes in the macrosheds df
  ws_ms_solutes <- colnames(ws_flux)[4:length(colnames(ws_flux))]
  ws_ms_solutes <- str_extract(ws_ms_solutes, "[^_]+")
  colnames(ws_flux)[4:length(colnames(ws_flux))] <- ws_ms_solutes

  # solutes in the published df
  ws_pub_solutes <- colnames(ws_pub)[4:length(colnames(ws_pub))]
  ws_pub_solutes <- str_extract(ws_pub_solutes, "[^_]+")
  colnames(ws_pub)[4:length(colnames(ws_pub))] <- ws_pub_solutes


  # solutes shared?
  ws_solutes <- ws_pub_solutes[ws_pub_solutes %in% ws_ms_solutes]
  ms_solutes <- ws_ms_solutes[ws_ms_solutes %in% ws_solutes]
  pub_solutes <- ws_pub_solutes[ws_pub_solutes %in% ws_solutes]

  # filter dataframes
  ws_ms <- ws_flux %>%
    select(wy, site_code, method, any_of(ms_solutes))

  # for each solute in ws record
  for(solute in unique(ws_flux$))

    ws_solute_flux <- ws_flux %>%
      select(wy, site_code, method, !!solute) %>%
      filter(!is.na(!!solute))

# create annual flux comparison plots, including published fluxes
#   add water flux as area behind
# create pairwise comparison regressions between methods
# create all-time cumulative yield bar plot
}
