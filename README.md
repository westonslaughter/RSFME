# flux_estimation_methods_evaluation

this repo is undergoing restructuring

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
