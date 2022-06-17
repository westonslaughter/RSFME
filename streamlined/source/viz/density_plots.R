# script plotting USGS and MacroSHeds density by attributes
library(macrosheds)
library(ggplot2)
library(ggthemes)
library(wesanderson)

source("streamlined/source/usgs_helpers.R")

# USGS
usgs <- read_csv("streamlined/data/site/usgs_nitrate_sites.csv")
usgs_sites <- c(1,3,8,9,12,13,14,17,18)
usgs_sites <- usgs[usgs_sites,]$site_code

usgs.df <- get_site_ecoregions(site_data = "streamlined/data/site/usgs_nitrate_sites.csv",
                    eco_fp = "data/spatial/eco/NA_CEC_Eco_Level2.shp",
                    good_sites = usgs_sites)

usgs.df <- usgs.df %>%
  select(site_code, ws_area_ha, RBI, NA_L2NAME, geometry) %>%
  rename(rbi = RBI, ecoregion = NA_L2NAME)

st_write(usgs.df, "streamlined/data/site/usgs_site_info.csv", layer_options = "GEOMETRY=AS_XY")

# MacroSheds
# load in a df of MS sites with
ms <- ms_download_site_data() %>%
  select(site_code, ws_area_ha, latitude, longitude) %>%
  rename(lat = latitude,
         long = longitude)

# find all ms sites with nitrate
ms_vars <- ms_catalog()
ms_sites <- ms_vars[grep("^Nitrate", ms_vars$variable_name),]

# calculate the amount of days in ten years
ten_yrs_days <- 365 * 10

ms_sites <- ms_sites %>%
    mutate(period = as.numeric(
        difftime(last_record_utc,
                 first_record_utc,
                 units = "days"
                 ))) %>%
    # to get records > 60 years
  filter(period > ten_yrs_days)

# get Q for RBI calc
ms_sites_ls <- unique(ms_sites$site_code)
my_q <- ms_load_product(
    my_ms_dir,
    prodname = "discharge",
    site_codes = ms_sites_ls,
    sort_result = TRUE,
    warn = TRUE
)

# The current dataset has duplicates of observations, the new version of the dataset
# does not have that issue but is not on figshare yet. Sorry about that!
my_q <- my_q %>%
    distinct(datetime, site_code, .keep_all = T)

rbi_q <- my_q %>%
  group_by(site_code) %>%
  summarise(rbi = ContDataQC::RBIcalc(val))

ms.df <- ms %>%
  inner_join(rbi_q)

ms.df <- get_site_ecoregions(site_data = ms.df,
                    eco_fp = "data/spatial/eco/NA_CEC_Eco_Level2.shp")

ms.df <- ms.df %>%
  select(site_code, ws_area_ha, rbi, NA_L2NAME, geometry) %>%
  rename(ecoregion = NA_L2NAME)

st_write(ms.df, "streamlined/data/site/ms_site_info.csv", layer_options = "GEOMETRY=AS_XY")

ms.f <- ms.df %>% st_drop_geometry() %>% mutate(dataset = "macrosheds")
usgs.f <- usgs.df %>% st_drop_geometry() %>% mutate(dataset = "USGS")

sites <- rbind(ms.f, usgs.f)
sites$ws_area_ha <- as.numeric(sites$ws_area_ha)
sites$rbi <- as.numeric(sites$rbi)

# Plots
names(wes_palettes)

# n
ms.n <- nrow(sites[sites$dataset == 'macrosheds',])
usgs.n <- nrow(sites[sites$dataset == 'USGS',])

# Watershed Area
gg_wa <- ggplot(sites, aes(x=log(ws_area_ha), color = dataset)) +
  geom_density(aes(fill= dataset), alpha=0.7) +
  ggtitle("\nDistribution of Watershed Areas in MacroSheds and USGS",
          subtitle = "across sites used in flux analysis\n\n"
          ) +
  scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
  scale_color_manual(values = wes_palette("Royal1", n = 2)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("\nLog of Watershed Area (ha)\n") +
  ylab("\nDensity\n") +
  theme_few() +
  theme(text = element_text(size=28),
        axis.line = element_line(colour = "black"),
        plot.subtitle = element_text(face = "italic"),
        panel.border = element_blank(),
        legend.key.size = unit(1, 'cm'), #change legend key size
        ## legend.key.height = unit(1, 'cm'), #change legend key height
        ## legend.key.width = unit(1, 'cm'), #change legend key width
        ## legend.position = 'bottom',
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=24),
        legend.spacing.y = unit(1.0, 'cm')
        ) +
  annotate("text", x=4, y=0.08, size = 11,
           label= paste0("n = ", ms.n)) +
  annotate("text", x=10, y=0.08, size = 11,
           label= paste0("n = ", usgs.n))

# Flashiness (RBI)
gg_rbi <- ggplot(sites, aes(x=rbi, color = dataset)) +
  geom_density(aes(fill= dataset), alpha=0.7) +
  ggtitle("\nDistribution of Hydrologic Flashiness (RBI) in MacroSheds and USGS",
          subtitle = "across sites used in flux analysis\n\n"
          ) +
  scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
  scale_color_manual(values = wes_palette("Royal1", n = 2)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("\nRichards-Baker (Flashiness) Index (RBI)"^1*"\n")) +
  ylab("\nDensity\n") +
  theme_few() +
  theme(text = element_text(size=28),
        axis.line = element_line(colour = "black"),
        plot.subtitle = element_text(face = "italic"),
        panel.border = element_blank(),
        legend.key.size = unit(1, 'cm'), #change legend key size
        ## legend.key.height = unit(1, 'cm'), #change legend key height
        ## legend.key.width = unit(1, 'cm'), #change legend key width
        ## legend.position = 'bottom',
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=24),
        legend.spacing.y = unit(1.0, 'cm')
        ) +
  annotate("text", x=0.5, y=1.2, size = 11,
           label= paste0("n = ", ms.n)) +
  annotate("text", x=0.13, y=1.2, size = 11,
           label= paste0("n = ", usgs.n)) +
  labs(caption = expression(""^1*"Baker et al. 2004"))

# EcoRegions
ggplot(sites, aes(x = `rbi`, y = `ecoregion`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "rbi", option = "C") +
  labs(title = 'Temperatures in Lincoln NE in 2016') +
  theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12)
    )



ggplot(sites, aes(x = `ws_area_ha`, y = `ecoregion`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "ws_area_ha", option = "C") +
  labs(title = 'Watershed Area by EcoRegion') +
  theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    )
