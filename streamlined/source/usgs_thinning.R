#### Thinning Functions for USGS Data ####
# Load in packages


# set thinning interval options and variable
thinning_intervals <- c('daily', 'weekly', 'biweekly', 'monthly')
var <- 'nitrate_nitrite_mgl'

# read in df of USGS sites w continous Nitrate
usgs <- read.csv("streamlined/data/site/usgs_nitrate_sites.csv",
                 colClasses = "character")

# loop thru sites and run thinning functions
for(i in 1:nrow(site_var_data)) {
    # set USGS site
    site <- usgs$site_code[i]

    # get discharge
    q_data <- try(read_feather(glue('data/raw/q_cfs/{site}.feather')))
    if(inherits(q_data, 'try-error')) next

    # get chemistry
    chem_data <- try(read_feather(glue('data/raw/{var}/{site}.feather')))
    if(inherits(chem_data, 'try-error')) next

    for(p in 1:length(thinning_intervals)){

        if(thinning_intervals[p] == 'daily'){
            usgs_thin_daily()
        }

        if(thinning_intervals[p] == 'weekly'){
            usgs_thin_weekly()
        }

        if(thinning_intervals[p] == 'biweekly'){
            usgs_thin_biweekly()
        }

        if(thinning_intervals[p] == 'monthly'){
            usgs_thin_monthly()
        }
    }
}
