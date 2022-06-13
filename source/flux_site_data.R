library(sf)

# read in sites
site_data <- read.csv("data/general/site_data.csv")
site_data_sf = st_as_sf(site_data, coords = c("long", "lat"),
                 crs = 4326)


# get site ecogeographic region
## epa_eco_ii <- download.file("https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip",
##                             "data/spatial/epa_eco_ii.zip")
epa_eco <- read_sf("data/spatial/NA_CEC_Eco_Level2.shp")
epa_eco_ii <- epa_eco %>%
  select(NA_L2NAME, geometry) %>%
  st_transform(crs = 4326)

# if needed
## sf_use_s2(FALSE)
## epa_eco_ii <- st_make_valid(epa_eco_ii)
## site_data_sf <- st_make_valid(site_data_sf)

# make a feature collection of sites and their ecoregions
site_eco <- st_intersection(site_data_sf, epa_eco_ii)

# get site flashiness index
