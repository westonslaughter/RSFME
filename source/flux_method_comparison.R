library(ggplot2)
library(grid)
library(gridExtra)

hbef_ws <- c('w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9')


# all HBEF
ws.data <- w9

# remove incredible outliers
ws.outliers <- ws.data$val > 200
ws.data <- ws.data[!ws.outliers,]

ws.data$wy <- as.Date(ws.data$wy, format = "%Y")
breaks.vec <- seq(min(ws.data$wy), max(ws.data$wy), by = "5 years")

# make list of plots
hbef_flux_plots <- list()
index = 1

# ppulatewith each WS data
for(ws in hbef_ws) {
    # title

    # report % of WRTDS measurements removed
    print(length(ws.outliers[ws.outliers == TRUE]))


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
      ggtitle(ws)

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
  vp=viewport(width=0.9, height=0.9))
