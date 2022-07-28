library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)

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

    # plot remainder
    usgs_flux_plots[[thinning]][[ws]] <- ggplot(ws.data[ws.data$site_code == ws & ws.data$thin == thinning,], aes(x=wy, y = flux)) +
      geom_point(aes(color = method), size = 3.5) +
      theme_minimal() +
      ylim(0, 160) +
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
     usgs_flux_plots[[thinning]][[ws]]

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
             nrow = 4, ncol = 4,
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
             nrow = 4, ncol = 4,
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
             nrow = 4, ncol = 4,
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
             nrow = 4, ncol = 4,
             top = textGrob(plt_title, gp=gpar(fontsize=26, font=3), vjust = -2),
             vp = viewport(width=0.9, height=0.9))
