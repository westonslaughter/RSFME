library(sf)

# read in sites
site_data <- read.csv("data/general/site_data.csv")
site_data_sf = st_as_sf(site_data, coords = c("long", "lat"),
                 crs = 4326)


# get site ecogeographic region
## (use CURL options with care)
options(download.file.method="curl", download.file.extra="-k -L")
eco_fp <- paste0(getwd(), "/data/spatial/")
eco_f <- paste0(eco_fp, "epa_eco.zip")

download.file("https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip",
                            eco_f)
unzip(eco_f, exdir = eco_fp)

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

# save
st_write(site_eco, "data/general/site_eco.csv", layer_options = "GEOMETRY=AS_XY")

## get_epa_ecoregion <- function(site_data, ecogeo_data, level) {


## }
