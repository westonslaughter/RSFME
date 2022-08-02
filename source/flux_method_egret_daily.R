# calculate annual flux for USGS and MacroSheds sites, by site and site year
# using multiple methods: beale, period weighted, rating, composite, WRTDS

# laod libraries
library(here)
library(dplyr)
library(feather)
library(zoo)
library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(glue)
library(feather)
library(zoo)
library(lfstat)
library(RiverLoad)
library(ContDataQC)
library(EGRET)
library(macrosheds)
library(stringr)

source('source/usgs_helpers.R')
source('source/helper_functions.R')
source('source/flux_methods.R')
source('source/egret_overwrites.R')
# load source files


# MAIN LOOP
# calculate annual flux for USGS and MacroSheds sites
datasets <- c('USGS', 'MacroSheds')

# USGS
# using multiple methods: beale, period weighted, rating, composite, WRTDS
# by site-year
# by site
