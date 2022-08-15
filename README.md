# flux_estimation_methods_evaluation

TODO:
- make single USGS and HBEF annual flux calcs loop, which runs flux calcs by site-year, and also runs calcs for each site with all available years
- run this loop and share / first look at results
- ~~look at HBEF Watershed 3, run full annual flux estimates for all site years, each site-year seperately, and ocmpare to 'true flux' from sensor data for NO3_N~~
- **look at HBEF Watershed 3, run full annual flux estimates for all site years, each site-year seperately, and ocmpare to 'true flux' from sensor data for Ca AND HBEF offficial published estimate for Ca flux**
- WRTDS results (currently) make sense as-is, need to rigorously document UNITS and ensure that it all maths out to kg/day
- look into model paramter 'tuning' does it have built in for settting a Max Q or Max Flux? If not, post processing?
- look at fdaily time series of years where extreme outlier flux overestimates are occuring- what kind of days/flows are these overestimates occuring at?  
- diagnostic plots such as showing daily flux with the days of sampling indicated, and printed statistic about the 'balancedness' of the sampels across the year. 
- compare plot Q:C by WY vs Q:C over all site data
- average annual hydrograph against hydrograph of any particular year
- improve USGS flux comparison plots, currently impossible to tell outliers or diff in good ones even
  


## purpose
The purpose of this repo is to develop and maintain tools to evaluate how to best calculate riverine solute fluxes. 

## structure

### main script
The master script is written to coordinate between all subscripts and data sources.

### data 
 Most data in this analysis will be pulled from USGS servers during each analysis.

All data that is not pulled from in-script web requests will reside in the data folder. The data folder is organized following the same procedure as Macrosheds (https://github.com/MacroSHEDS/macrosheds), nested by domain and site. 

### source
Each method of flux estimation is replicated in it's own subscript. These scripts and any associated helper functions are stored in the source folder.

### out
The output from selected methods and visualizations specificed in the master script are output to the out folder.
