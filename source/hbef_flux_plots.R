
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
ms_data <- list.files('data/nice/stream_flux/hbef', full.names = TRUE)
hbef_raw <- do.call(rbind,lapply(ms_data,read_feather))

hbef_flux <- hbef_raw %>%
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


# create detect outlier function (if x is 1.5x the interquartile range higher than 3rd Quart or lower than 1st quart)
detect_outlier <- function(x) {
    # calculate first quantile
    Quantile1 <- quantile(x, probs=.25)
    # calculate third quantile
    Quantile3 <- quantile(x, probs=.75)
    # calculate inter quartile range
    IQR = Quantile3-Quantile1
    # return true or false
    x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}

# plotting function
flux_compare_plot <- function(data, watershed, solute) {
  # look at flux time series
  fluxpal <- rev(brewer.pal(n=7, name='Dark2'))
  data <- data %>%
    rename(
      solute = !!solute
    )

  solute_plot <- ggplot() +
    geom_point(data %>% filter(method == 'published'), mapping = aes(x=wy, y=solute, color = method, size = method, shape = method, fill = method)) +
    geom_line(data %>% filter(method != 'published'), mapping = aes(x=wy, y=solute, group = method, color = method)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
        ## legend.position="none",
        text = element_text(size = 26),
        plot.title = element_text(size = 24, face = "bold")) +
    scale_shape_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'published'),
                     values = c(16, 16, 16, 16, 16, 16, 24)) +
    scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                     values = fluxpal) +
    scale_fill_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                     values = fluxpal) +
    scale_size_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                    values = c(2, 2, 2, 2, 2, 2, 6)) +
    ggtitle(glue("Hubbard Brook, {watershed}, 1960's-2020,\n {solute} Flux Estimates by Various Methods ",
                 watershed = watershed, solute = solute
                 )) +
    ylab('Annual Cumulative Solute Flux (kg/ha) \n ') +
    xlab('\n Water Year (Oct-Sep)') +
    scale_x_discrete(breaks=seq(1960, 2020, 10))

  return(solute_plot)
}

library(plotly)
network = 'lter'
domain = 'hbef'
for(ws in names(hbef_links)) {
  plot_fp <- glue('flux_plots/{network}/{domain}/{site}/flux_methods_timeseries/',
                  network = network,
                  domain = domain,
                  site = watershed)
  dir.create(plot_fp, recursive = TRUE)

  ws_filename <- paste0(ws, '.csv')
  ws_fp = file.path(hbef_pubs_dir, ws_filename)

  download.file(url = hbef_links[ws], ws_fp)

  ws_pub <- read.csv(ws_fp)
  ws_flux <- hbef_flux %>%
    filter(site_code == ws)

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

  ws_pub <- ws_pub %>%
    mutate(site_code = ws,
           wy = water_year(paste0(Year_Month, '-01')),
           method = 'published'
           )
    select(wy, site_code, method, any_of(ms_solutes))

  fluxes <- colnames(ws_flux)[4:ncol(ws_flux)]

  # for each solute in ws record
  for(solute in fluxes) {
    ws_pub <-

    ws_solute_flux <- ws_flux %>%
      select(wy, site_code, method, !!solute) %>%
      filter(!is.na(!!solute))

    site_solute_plot <- flux_compare_plot(ws_solute_flux, ws, solute)
    ggsave(paste0(plot_fp, 'solute.png'), site_solute_plot, bg = 'white')
  }


# create annual flux comparison plots, including published fluxes
#   add water flux as area behind
# create pairwise comparison regressions between methods
# create all-time cumulative yield bar plot
}
