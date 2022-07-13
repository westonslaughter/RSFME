library(macrosheds)

my_ms_dir <- here("streamlined/data/ms")

options(timeout = 60) #default 60 might not be enough if your connection is slow

ms_sites <- ms_download_site_data()
ms_networks <- unique(ms_sites$network)

ms_download_core_data(
    my_ms_dir,
    networks = ms_networks,
    quiet = FALSE
)
