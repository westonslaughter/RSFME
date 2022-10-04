library(ggplot2)
library(dplyr)
library(feather)
library(grid)
library(gridExtra)
library(ggrepel)
library(lfstat)

source('./source/flux_methods.R')
source('./source/helper_functions.R')



# all HBEF
hbef_ws <- c('w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9')
ws.data <- read_feather('data/ms/hbef/stream_flux/w9.feather')

ws.data$wy <- as.Date(ws.data$wy, format = "%Y")
breaks.vec <- seq(min(ws.data$wy), max(ws.data$wy), by = "5 years")

# make list of plots
hbef_flux_plots <- list()
index = 1

# ppulatewith each WS data
for(ws in hbef_ws) {
    # find major outliers
    ws.this.outliers <- which(ws.data$site_code == ws & ws.data$val > 200)

    # report number of WRTDS measurements removed
    ws.n.otl <- length(ws.this.outliers[ws.this.outliers])
    ws.n.text <- paste(ws, 'outliers:', ws.n.otl)
    print(ws.n.text)

    # remove these outliers
    ws.data <- ws.data[-ws.this.outliers,]

    # plot remainder
    hbef_flux_plots[[ws]] <- ggplot(ws.data[ws.data$site_code == ws,], aes(x=wy, y = val)) +
      geom_point(aes(color = method), size = 3.5) +
      ## geom_line(aes(group = method, color = method)) +
      theme_minimal() +
      ylim(0, 150) +
      scale_x_date(breaks = breaks.vec, limits = c(breaks.vec[1], breaks.vec[length(breaks.vec)]), date_labels = "%Y") +
      theme(panel.grid.major = element_blank(),
            legend.position="none",
            text = element_text(size = 24),
            plot.title = element_text(size = 24, face = "bold")) +
      ggtitle(ws) +
      annotate(geom="text",
               x=breaks.vec[length(breaks.vec)-2],
               y=120,
               label=ws.n.text,
               size = 6,
               color="black")

  if(!index %in% c(1, 4, 7)){
    hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank())
  }

  if(index %in% c(1, 7)) {
      hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
        theme(axis.title.y = element_blank())
  } else if(index == 4) {
      hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
        ylab('Nitrate Nitrogen (Kg/yr)')
  }

  if(!index %in% c(7, 8, 9) ){
    hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
  }

    index = index + 1
}

plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in HBEF Watersheds 1-9")
hbef_annual_flux_methods <- grid.arrange(
             hbef_flux_plots[['w1']],
             hbef_flux_plots[['w2']],
             hbef_flux_plots[['w3']] + theme(legend.position = 'top'),
             hbef_flux_plots[['w4']],
             hbef_flux_plots[['w5']],
             hbef_flux_plots[['w6']],
             hbef_flux_plots[['w7']],
             hbef_flux_plots[['w8']],
             hbef_flux_plots[['w9']],
             nrow = 3, ncol = 3,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))

# all USGS
ws.data <- read_feather('out/usgs_annual_flux.feather')
ws.data <- ws.data[ws.data$site_code != "site_no",]

ws.data$wy <- as.Date(ws.data$wy, format = "%Y")
breaks.vec <- seq(min(ws.data$wy), max(ws.data$wy), by = "1 years")

usgs_ws <- unique(ws.data$site_code)

# make list of plots
usgs_flux_plots <- list()
index = 1

# loop thru dif thinnings
thinnings <- c('daily', 'weekly', 'biweekly', 'monthly')
for(thinning in thinnings) {
  ## ws.data <- ws.data.all[ws.data.all$thin == thinning,]


  # ppulatewith each WS data
  for(ws in usgs_ws) {
    # find major outliers
    ## ws.this.outliers <- which(ws.data$site_code == ws & ws.data$flux > 200)
    ## # report number of WRTDS measurements removed
    ## ws.n.otl <- length(ws.this.outliers[ws.this.outliers])
    ## ws.n.text <- paste(ws, 'outliers:', ws.n.otl)
    ## print(ws.n.text)
    # remove these outliers
    ## ws.data <- ws.data[-ws.this.outliers,]
    ##

    ## thinning = 'daily'
    ## ws = '07374000'

    # plot remainder
    usgs_flux_plots[[thinning]][[ws]] <- ggplot(ws.data[ws.data$site_code == ws & ws.data$thin == thinning,], aes(x=wy, y = flux)) +
      geom_point(aes(color = method), size = 3.5, alpha = 0.7) +
      theme_minimal() +
      ylim(0, 200) +
      scale_x_date(breaks = breaks.vec, limits = c(breaks.vec[1], max(ws.data$wy)), date_labels = "%Y") +
      theme(panel.grid.major = element_blank(),
            legend.position="none",
            text = element_text(size = 24),
            plot.title = element_text(size = 24, face = "bold")) +
      ggtitle(ws)
      ## +
      ## annotate(geom="text",
      ##          x=breaks.vec[length(breaks.vec)],
      ##          y=120,
      ##          label=ws.n.text,
      ##          size = 6,
      ##          color="black")
     ## usgs_flux_plots[[thinning]][[ws]]

  ## if(!index %in% c(1, 4, 7)){
  ##   usgs_flux_plots[[thinning]][[ws]] <- usgs_flux_plots[[thinning]][[ws]] +
  ##     theme(axis.title.y = element_blank(),
  ##           axis.text.y = element_blank())
  ## }

  ## if(index %in% c(1, 7)) {
  ##     usgs_flux_plots[[thinning]][[ws]] <- usgs_flux_plots[[thinning]][[ws]] +
  ##       theme(axis.title.y = element_blank())
  ## } else if(index == 4) {
  ##     usgs_flux_plots[[thinning]][[ws]] <- usgs_flux_plots[[thinning]][[ws]] +
  ##       ylab('Nitrate Nitrogen (Kg/yr)')
  ## }

  ## if(!index %in% c(7, 8, 9) ){
  ##   usgs_flux_plots[[thinning]][[ws]] <- usgs_flux_plots[[thinning]][[ws]] +
  ##     theme(axis.title.x = element_blank(),
  ##           axis.text.x = element_blank())
  ## }

    index = index + 1
  }
}

plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Monthly Thin")
usgs_annual_flux_methods_monthly <- grid.arrange(
             usgs_flux_plots[['monthly']][[usgs_ws[1]]],
             usgs_flux_plots[['monthly']][[usgs_ws[2]]],
             usgs_flux_plots[['monthly']][[usgs_ws[3]]] + theme(legend.position = 'top'),
             usgs_flux_plots[['monthly']][[usgs_ws[4]]],
             usgs_flux_plots[['monthly']][[usgs_ws[5]]],
             usgs_flux_plots[['monthly']][[usgs_ws[6]]],
             usgs_flux_plots[['monthly']][[usgs_ws[7]]],
             usgs_flux_plots[['monthly']][[usgs_ws[8]]],
             usgs_flux_plots[['monthly']][[usgs_ws[9]]],
             usgs_flux_plots[['monthly']][[usgs_ws[10]]],
             usgs_flux_plots[['monthly']][[usgs_ws[11]]],
             usgs_flux_plots[['monthly']][[usgs_ws[12]]],
             usgs_flux_plots[['monthly']][[usgs_ws[13]]],
             usgs_flux_plots[['monthly']][[usgs_ws[14]]],
             usgs_flux_plots[['monthly']][[usgs_ws[15]]],
             usgs_flux_plots[['monthly']][[usgs_ws[16]]],
             usgs_flux_plots[['monthly']][[usgs_ws[17]]],
             usgs_flux_plots[['monthly']][[usgs_ws[18]]],
             usgs_flux_plots[['monthly']][[usgs_ws[19]]],
             usgs_flux_plots[['monthly']][[usgs_ws[20]]],
             usgs_flux_plots[['monthly']][[usgs_ws[21]]],
             usgs_flux_plots[['monthly']][[usgs_ws[22]]],
             usgs_flux_plots[['monthly']][[usgs_ws[23]]],
             usgs_flux_plots[['monthly']][[usgs_ws[24]]],
             usgs_flux_plots[['monthly']][[usgs_ws[25]]],
             usgs_flux_plots[['monthly']][[usgs_ws[26]]],
             usgs_flux_plots[['monthly']][[usgs_ws[27]]],
             usgs_flux_plots[['monthly']][[usgs_ws[28]]],
             usgs_flux_plots[['monthly']][[usgs_ws[29]]],
             usgs_flux_plots[['monthly']][[usgs_ws[30]]],
             usgs_flux_plots[['monthly']][[usgs_ws[31]]],
             usgs_flux_plots[['monthly']][[usgs_ws[32]]],
             usgs_flux_plots[['monthly']][[usgs_ws[33]]],
             usgs_flux_plots[['monthly']][[usgs_ws[34]]],
             usgs_flux_plots[['monthly']][[usgs_ws[35]]],
             usgs_flux_plots[['monthly']][[usgs_ws[36]]],
             usgs_flux_plots[['monthly']][[usgs_ws[37]]],
             usgs_flux_plots[['monthly']][[usgs_ws[38]]],
             usgs_flux_plots[['monthly']][[usgs_ws[39]]],
             usgs_flux_plots[['monthly']][[usgs_ws[40]]],
             usgs_flux_plots[['monthly']][[usgs_ws[41]]],
             usgs_flux_plots[['monthly']][[usgs_ws[42]]],
             usgs_flux_plots[['monthly']][[usgs_ws[43]]],
             usgs_flux_plots[['monthly']][[usgs_ws[44]]],
             usgs_flux_plots[['monthly']][[usgs_ws[45]]],
             usgs_flux_plots[['monthly']][[usgs_ws[46]]],
             usgs_flux_plots[['monthly']][[usgs_ws[47]]],
             usgs_flux_plots[['monthly']][[usgs_ws[48]]],
             usgs_flux_plots[['monthly']][[usgs_ws[49]]],
             usgs_flux_plots[['monthly']][[usgs_ws[50]]],
             usgs_flux_plots[['monthly']][[usgs_ws[51]]],
             usgs_flux_plots[['monthly']][[usgs_ws[52]]],
             usgs_flux_plots[['monthly']][[usgs_ws[53]]],
             usgs_flux_plots[['monthly']][[usgs_ws[54]]],
             usgs_flux_plots[['monthly']][[usgs_ws[55]]],
             usgs_flux_plots[['monthly']][[usgs_ws[56]]],
             usgs_flux_plots[['monthly']][[usgs_ws[57]]],
             usgs_flux_plots[['monthly']][[usgs_ws[58]]],
             usgs_flux_plots[['monthly']][[usgs_ws[59]]],
             usgs_flux_plots[['monthly']][[usgs_ws[60]]],
             usgs_flux_plots[['monthly']][[usgs_ws[61]]],
             nrow = 8, ncol = 8,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))

plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Weekly Thin")
usgs_annual_flux_methods_monthly <- grid.arrange(
             usgs_flux_plots[['weekly']][[usgs_ws[1]]],
             usgs_flux_plots[['weekly']][[usgs_ws[2]]],
             usgs_flux_plots[['weekly']][[usgs_ws[3]]] + theme(legend.position = 'top'),
             usgs_flux_plots[['weekly']][[usgs_ws[4]]],
             usgs_flux_plots[['weekly']][[usgs_ws[5]]],
             usgs_flux_plots[['weekly']][[usgs_ws[6]]],
             usgs_flux_plots[['weekly']][[usgs_ws[7]]],
             usgs_flux_plots[['weekly']][[usgs_ws[8]]],
             usgs_flux_plots[['weekly']][[usgs_ws[9]]],
             usgs_flux_plots[['weekly']][[usgs_ws[10]]],
             usgs_flux_plots[['weekly']][[usgs_ws[11]]],
             nrow = 8, ncol = 8,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))

plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Biweekly Thin")
usgs_annual_flux_methods_monthly <- grid.arrange(
             usgs_flux_plots[['biweekly']][[usgs_ws[1]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[2]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[3]]] + theme(legend.position = 'top'),
             usgs_flux_plots[['biweekly']][[usgs_ws[4]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[5]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[6]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[7]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[8]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[9]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[10]]],
             usgs_flux_plots[['biweekly']][[usgs_ws[11]]],
             nrow = 8, ncol = 8,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))


plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Daily Thin")
usgs_annual_flux_methods_monthly <- grid.arrange(
             usgs_flux_plots[['daily']][[usgs_ws[1]]],
             usgs_flux_plots[['daily']][[usgs_ws[2]]],
             usgs_flux_plots[['daily']][[usgs_ws[3]]] + theme(legend.position = 'top'),
             usgs_flux_plots[['daily']][[usgs_ws[4]]],
             usgs_flux_plots[['daily']][[usgs_ws[5]]],
             usgs_flux_plots[['daily']][[usgs_ws[6]]],
             usgs_flux_plots[['daily']][[usgs_ws[7]]],
             usgs_flux_plots[['daily']][[usgs_ws[8]]],
             usgs_flux_plots[['daily']][[usgs_ws[9]]],
             usgs_flux_plots[['daily']][[usgs_ws[10]]],
             usgs_flux_plots[['daily']][[usgs_ws[11]]],
             nrow = 8, ncol = 8,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))

# HBEF, Watershed 3 Deep Dive
ws <- 'w3'
ws.data <- read_feather('data/ms/hbef/stream_flux/w9.feather')
ws.data$wy <- as.Date(ws.data$wy, format = "%Y")
breaks.vec <- seq(min(ws.data$wy), max(ws.data$wy), by = "5 years")

ws.q <- read_feather('data/ms/hbef/discharge/w3.feather')

# load in 'raw' sample chem from w3
w3_samples <- read_feather('data/ms/hbef/true//w3_chem_samples.feather')                                        #

# load in daily WRTDS
w3_2012 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2012.feather')
w3_2013 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2013.feather')
w3_2014 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2014.feather')
w3_2015 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2015.feather')
w3_2016 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2016.feather')
w3_2017 <- read_feather('data/ms/hbef/true/w3_dailyWRTDS_2017.feather')

w3_wrtds_daily <- rbind(w3_2012,w3_2013,w3_2014,w3_2015,w3_2016,w3_2017)

w3_wrtds_daily <- w3_wrtds_daily %>%
  mutate(site_code = 'w3',
         ## var = 'GN_NO3_N',
         ## ms_recommended = 0,
         FluxDay = FluxDay,
         method = 'wrtds'
         ) %>%
  select(Date, FluxDay, site_code, method) %>%
  rename(date = Date, val = FluxDay)



# load in 'true flux' for watershed 3
w3_true <- read_feather("data/ms/hbef/true//w3_sensor_wdisch.feather")
w3_true_flux_daily <- w3_true %>%
  group_by(date) %>%
  # convert to Nitrate_Nitrogen
  mutate(IS_NO3 = (IS_NO3 * .2259)) %>%
  summarise(val = sum(IS_NO3, na.rm = TRUE)/1000) %>%
  mutate(site_code = 'w3',
         ## var = 'GN_NO3_N',
         ## ms_recommended = 0,
         method = 'true'
         )

# w3 Daily
w3_daily <- rbind(w3_true_flux_daily, w3_wrtds_daily)

w3_true_flux <- w3_true %>%
  group_by(date) %>%
  # convert to Nitrate_Nitrogen
  mutate(IS_NO3 = IS_NO3 * .2259) %>%
  summarise(val = sum(IS_NO3, na.rm = TRUE)) %>%
  mutate(wy = water_year(date)) %>%
  group_by(wy) %>%
  summarise(val = sum(val, na.rm = TRUE)/1000) %>%
  mutate(site_code = 'w3',
         var = 'GN_NO3_N',
         method = 'true',
         ms_recommended = 0)

w3_true_flux$wy <- as.Date(w3_true_flux$wy, format = "%Y")

# bing true flux onto df
ws.data <- rbind(ws.data, w3_true_flux)

# find major outliers
ws.this.outliers <- which(ws.data$site_code == ws & ws.data$val > 25)
ws.outlier.wys <- ws.data[ws.data$site_code == ws &ws.data$val > 25,]$wy

# report number of WRTDS measurements removed
ws.n.otl <- length(ws.this.outliers[ws.this.outliers])
ws.n.text <- paste(ws, 'outliers:', ws.n.otl)
print(ws.n.text)

# remove these outliers
## ws.data <- ws.data[-ws.this.outliers,]

# bring in monthly flux from HBEF website
hbef_flux_official <- read.csv('data/raw/hbef_published_flux/ws3_stream_monthly_flux_gHa.csv') %>%
  mutate(date = paste0(Year, '-', Month, '-', '01'),
         wy = water_year(date)) %>%
  select(wy, NO3_N_flux)


# colors
library(RColorBrewer)
fluxpal <- brewer.pal(n=7, name='Dark2')
## fluxalpha <- paste0(fluxpal, "90")
## fluxalpha[length(fluxalpha)] = fluxpal[length(fluxpal)]


# plot remainder
# optional, only 'avg 'true' and 'wrtds'
fluxalpha <- fluxpal
fluxalpha[c(2:4, 6)] <- paste0(fluxpal[c(2:4, 6)], "00")
fluxalpha[1] <- '#A7226E'
fluxalpha[5] <- '#2f9599'
fluxalpha[7] <- '#ec2049'

w3_plot <- ggplot(ws.data[ws.data$site_code == ws,], aes(x=wy, y = val)) +
  geom_point(aes(color = method, shape = method, size = method)) +
  ## geom_line(aes(group = method, color = method)) +
  theme_minimal() +
  ylim(0, 15) +
  scale_x_date(breaks = breaks.vec, limits = c(breaks.vec[1], breaks.vec[length(breaks.vec)]), date_labels = "%Y") +
  theme(panel.grid.major = element_blank(),
        ## legend.position="none",
        text = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold")) +
  ggtitle(ws) +
  scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true'),
                     values = fluxalpha) +
  scale_shape_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true'),
                     values = c(20, 20, 20, 20, 20, 20, 4)) +
  scale_size_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true'),
                     values = c(4, 4, 4, 4, 4, 4, 6)) +
  annotate(geom="text",
               x=breaks.vec[length(breaks.vec)-2],
               y=120,
               label=ws.n.text,
               size = 6,
               color="black")+                                                                 # Draw vlines to plot
  geom_vline(xintercept = ws.outlier.wys,
             col = "red", lwd = 0.1)
w3_plot

# Daily
w3_ls <- list()

# get Q as well
ws.q <- ws.q %>%
  mutate(date = lubridate::date(datetime)) %>%
  rename(discharge = val)

# find major outliers
ws.this.outliers <- which(w3_daily$site_code == ws & w3_daily$val > 10)
outlier.dates <- w3_daily[w3_daily$val > 10,]$date

# report number of WRTDS measurements removed
ws.n.otl <- length(ws.this.outliers[ws.this.outliers])
ws.n.text <- paste(ws, 'outliers:', ws.n.otl)
print(ws.n.text)

# identify these outliers
w3_daily <- w3_daily %>%
  mutate(outlier = ifelse(val > 10, 'red', 'black')) %>%
  left_join(ws.q, by = c('date'))

w3_daily <- w3_daily %>%
  mutate(discharge_z = (discharge - mean(discharge)/sd(discharge)))

# plot daily flux values of true and WRTDS flux
w3_flux_daily_plot <- ggplot(w3_daily, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  # add red tick marks at days with flux vals > 10
  ylim(0, 10) +
  # should add flow w second axis as well?
  theme_minimal() +                                                                 # Draw vlines to plot
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.1)

w3_flux_daily_plot

# plot daily flux values of true and WRTDS flux
w3_q_daily_plot <- ggplot(w3_daily, aes(x=date, y = discharge)) +
  geom_point() +
  # add red tick marks at days with flux vals > 10
  ## ylim(-1, 1) +
  # should add flow w second axis as well?
  theme_minimal() +                                                                 # Draw vlines to plot
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.1)

w3_q_daily_plot

w3_ls[['flux']] = w3_flux_daily_plot
w3_ls[['Q']] = w3_q_daily_plot

## plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Daily Thin")
w3_flux_q <- grid.arrange(
             w3_ls[['flux']] + theme(legend.position = 'top'),
             w3_ls[['Q']],
             nrow = 2, ncol = 1,
             ## top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))

# just 2013
w3_13 <- w3_daily %>%
  filter(water_year(date) == "2013")

w3_13_samples <- w3_samples %>%
  filter(wy %in% '2013', !is.na(con))

sampledates <- unique(w3_13_samples$datetime)

w3_13_plot <- ggplot(w3_13, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  theme_minimal() +                                                                 # Draw vlines to plot
  ## ylim(0, 5) +
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.4) +
  geom_vline(xintercept = sampledates,
             col = "blue", lwd = 0.05, linetype = 5) +
  theme(
    panel.grid.major.y  = element_blank()
  ) + ggtitle("WRTDS vs True Flux, Watershed 3, 2013")

w3_13_plot

# just 2014
w3_14 <- w3_daily %>%
  filter(water_year(date) == "2014")

w3_14_samples <- w3_samples %>%
  filter(wy %in% '2014', !is.na(con))

sampledates <- unique(w3_14_samples$datetime)

w3_14_plot <- ggplot(w3_14, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  theme_minimal() +                                                                 # Draw vlines to plot
  ## ylim(0, 5) +
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.4) +
  geom_vline(xintercept = sampledates,
             col = "blue", lwd = 0.05, linetype = 5) +
  theme(
    panel.grid.major.y  = element_blank()
  ) + ggtitle("WRTDS vs True Flux, Watershed 3, 2014")

w3_14_plot


# just 2015
w3_15 <- w3_daily %>%
  filter(water_year(date) == "2015")

w3_15_samples <- w3_samples %>%
  filter(wy %in% '2015', !is.na(con))

sampledates <- unique(w3_15_samples$datetime)

w3_15_plot <- ggplot(w3_15, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  theme_minimal() +                                                                 # Draw vlines to plot
  ## ylim(0, 5) +
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.4) +
  geom_vline(xintercept = sampledates,
             col = "blue", lwd = 0.05, linetype = 5) +
  theme(
    panel.grid.major.y  = element_blank()
  ) + ggtitle("WRTDS vs True Flux, Watershed 3, 2015")

w3_15_plot



# just 2016
w3_16 <- w3_daily %>%
  filter(water_year(date) == "2016")

w3_16_samples <- w3_samples %>%
  filter(wy %in% '2016', !is.na(con))

sampledates <- unique(w3_16_samples$datetime)

w3_16_plot <- ggplot(w3_16, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  theme_minimal() +                                                                 # Draw vlines to plot
  ## ylim(0, 5) +
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.4) +
  geom_vline(xintercept = sampledates,
             col = "blue", lwd = 0.05, linetype = 5) +
  theme(
    panel.grid.major.y  = element_blank()
  ) + ggtitle("WRTDS vs True Flux, Watershed 3, 2016")

w3_16_plot




# just 2017
w3_17 <- w3_daily %>%
  filter(water_year(date) == "2017")

w3_17_samples <- w3_samples %>%
  filter(wy %in% '2017', !is.na(con))

sampledates <- unique(w3_17_samples$datetime)

w3_17_plot <- ggplot(w3_17, aes(x=date, y = val)) +
  geom_point(aes(color = method)) +
  theme_minimal() +                                                                 # Draw vlines to plot
  ## ylim(0, 5) +
  geom_vline(xintercept = outlier.dates,
             col = "red", lwd = 0.4) +
  geom_vline(xintercept = sampledates,
             col = "blue", lwd = 0.05, linetype = 5) +
  theme(
    panel.grid.major.y  = element_blank()
  ) + ggtitle("WRTDS vs True Flux, Watershed 3, 2017")

w3_17_plot


w3_dplt <- list()

              w3_dplt[['2013']] = w3_13_plot
              w3_dplt[['2014']] = w3_14_plot
              w3_dplt[['2015']] = w3_15_plot
              w3_dplt[['2016']] = w3_16_plot
              w3_dplt[['2017']] = w3_17_plot

## plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in USGS Watersheds, Daily Thin")
w3_flux_q <- grid.arrange(
             w3_dplt[['2013']] + theme(legend.position = 'top'),
             w3_dplt[['2014']],
             w3_dplt[['2015']],
             w3_dplt[['2016']],
             w3_dplt[['2017']],
             nrow = 5, ncol = 1,
             ## top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))



# getting w3 flux for all site data at once

target_year <- as.numeric(as.character(good_years))

raw_data_target_year <- raw_data_full %>%
            mutate(wy = as.numeric(as.character(wy))) %>%
            filter(wy %in% target_year)

        q_target_year <- raw_data_target_year %>%
            select(site_code, datetime, q_lps, wy)%>%
            na.omit()

        con_target_year <- raw_data_target_year %>%
            select(site_code, datetime, con, wy) %>%
            na.omit()

        ### calculate annual flux ######
        chem_df <- con_target_year
        q_df <- q_target_year

##### calculate WRTDS #####
calculate_wrtds <- function(chem_df, q_df, ws_size, lat, long, datecol = 'date') {
  tryCatch(
    expr = {
        egret_results <- adapt_ms_egret(chem_df, q_df, ws_size, lat, long, datecol = datecol)

        # still looking for reason why wrtds is 1K higher than others
        flux_from_egret <- egret_results$Daily %>%
          mutate(wy = water_year(Date)) %>%
          group_by(wy) %>%
          summarise(flux = warn_sum(FluxDay)/(area))
        },
    error = function(e) {
            print('ERROR: WRTDS failed to run')
            return(NA)
        })
    return(flux_from_egret)
}

write_feather(egret_results$Daily, 'data/ms/hbef/true/w3_wrtds_allyear_daily.feather')
write_feather(flux_from_egret, 'data/ms/hbef/true/w3_wrtds_allyear_annual.feather')

w3_allyr_flux <- flux_from_egret %>%
  rename(val = flux) %>%
  mutate(site_code = 'w3',
         var = 'GN_NO3_N',
         method = 'all_year',
         ms_recommended = 0)

# merge into ws.data
w3_allyr_flux$wy <- as.Date(w3_allyr_flux$wy, format = "%Y")
ws.data <- rbind(ws.data, w3_allyr_flux)


library(RColorBrewer)
fluxpal <- brewer.pal(n=8, name='Dark2')
## fluxalpha <- paste0(fluxpal, "90")
## fluxalpha[length(fluxalpha)] = fluxpal[length(fluxpal)]


# plot remainder
# optional, only 'avg 'true' and 'wrtds'
fluxalpha <- fluxpal
fluxalpha[c(1, 2, 6, 3, 4, 7)] <- paste0(fluxpal[c(1, 2, 6, 3, 4, 7)], "00")
## fluxalpha[2] <- '#A7226E'
fluxalpha[5] <- '#2f9599'
fluxalpha[7] <- '#ec2049'


w3_plot <- ggplot(ws.data[ws.data$site_code == ws & ws.data$val < 10,], aes(x=wy, y = val)) +
  geom_point(aes(color = method, shape = method, size = method)) +
  geom_line(data = ws.data[ws.data$site_code == ws & ws.data$val < 10 & ws.data$method != 'wrtds',],
            aes(color = method)) +
  theme_minimal() +
  ylim(0, 10) +
  scale_x_date(breaks = breaks.vec, limits = c(breaks.vec[1], breaks.vec[length(breaks.vec)]), date_labels = "%Y") +
  theme(panel.grid.major = element_blank(),
        ## legend.position="none",
        text = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold")) +
  ggtitle(ws) +
  scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true', 'all_year'),
                     values = fluxalpha) +
  scale_shape_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true', 'all_year'),
                     values = c(20, 20, 20, 20, 20, 20, 18, 20)) +
  scale_size_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'true', 'all_year'),
                     values = c(5, 5, 5, 5, 7, 5, 7, 5)) +
  annotate(geom="text",
               x=breaks.vec[length(breaks.vec)-2],
               y=120,
               label=ws.n.text,
               size = 6,
               color="black") +                                                                 # Draw vlines to plot
  geom_vline(xintercept = ws.outlier.wys,
             col = "darkgray", lwd = 0.25, linetype = 5)
w3_plot

# all HBEF
hbef_ws <- c('w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9')
ws.data <- read_feather('data/ms/hbef/stream_flux/w9.feather')

ws.data$wy <- as.Date(ws.data$wy, format = "%Y")
breaks.vec <- seq(min(ws.data$wy), max(ws.data$wy), by = "5 years")

# make list of plots
hbef_flux_plots <- list()
index = 1

# ppulatewith each WS data
for(ws in hbef_ws) {
    # find major outliers
    ws.this.outliers <- which(ws.data$site_code == ws & ws.data$val > 200)

    # report number of WRTDS measurements removed
    ws.n.otl <- length(ws.this.outliers[ws.this.outliers])
    ws.n.text <- paste(ws, 'outliers:', ws.n.otl)
    print(ws.n.text)

    # remove these outliers
    ws.data <- ws.data[-ws.this.outliers,]

    # plot remainder
    hbef_flux_plots[[ws]] <- ggplot(ws.data[ws.data$site_code == ws,], aes(x=wy, y = val)) +
      geom_point(aes(color = method), size = 3.5) +
      ## geom_line(aes(group = method, color = method)) +
      theme_minimal() +
      ylim(0, 150) +
      scale_x_date(breaks = breaks.vec, limits = c(breaks.vec[1], breaks.vec[length(breaks.vec)]), date_labels = "%Y") +
      theme(panel.grid.major = element_blank(),
            legend.position="none",
            text = element_text(size = 24),
            plot.title = element_text(size = 24, face = "bold")) +
      ggtitle(ws) +
      annotate(geom="text",
               x=breaks.vec[length(breaks.vec)-2],
               y=120,
               label=ws.n.text,
               size = 6,
               color="black")

  if(!index %in% c(1, 4, 7)){
    hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank())
  }

  if(index %in% c(1, 7)) {
      hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
        theme(axis.title.y = element_blank())
  } else if(index == 4) {
      hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
        ylab('Nitrate Nitrogen (Kg/yr)')
  }

  if(!index %in% c(7, 8, 9) ){
    hbef_flux_plots[[ws]] <- hbef_flux_plots[[ws]] +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
  }

    index = index + 1
}

plt_title <- paste("Annual Flux of Nitrate Nitrogen (kg) in HBEF Watersheds 1-9")
hbef_annual_flux_methods <- grid.arrange(
             hbef_flux_plots[['w1']],
             hbef_flux_plots[['w2']],
             hbef_flux_plots[['w3']] + theme(legend.position = 'top'),
             hbef_flux_plots[['w4']],
             hbef_flux_plots[['w5']],
             hbef_flux_plots[['w6']],
             hbef_flux_plots[['w7']],
             hbef_flux_plots[['w8']],
             hbef_flux_plots[['w9']],
             nrow = 3, ncol = 3,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))
