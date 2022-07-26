library(ggplot2)
library(gridExtra)

hbef_ws <- c('w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9')

hbef_flux_dfs <- list()
hbef_flux_plots <- list()

for(ws in hbef_ws) {
  fp <- paste0('data/ms/hbef/stream_flux/', ws, '.feather')
  assign(ws, read_feather(fp))

  hbef_flux_dfs[[ws]] <- get(ws)
}

for(ws in hbef_ws) {

    # remove incredible outliers
    ws.outliers <- hbef_flux_dfs[[ws]]$val > 200
    hbef_flux_dfs[[ws]] <- hbef_flux_dfs[[ws]][!ws.outliers,]

    # report % of WRTDS measurements removed
    print(length(ws.outliers[ws.outliers == TRUE]))

    # plot remainder
     hbef_flux_plots[[ws]]<- ggplot(hbef_flux_dfs[[ws]], aes(x=wy, y = val)) +
      geom_point(aes(color = method), size = 3.5) +
      ## geom_line(aes(group = method, color = method)) +
      theme_minimal() +
      theme(panel.grid.major = element_blank())
}

hbef_annual_flux_methods <- grid.arrange(hbef_flux_plots[['w1']],
             hbef_flux_plots[['w2']],
             hbef_flux_plots[['w3']],
             hbef_flux_plots[['w4']],
             hbef_flux_plots[['w5']],
             hbef_flux_plots[['w6']],
             hbef_flux_plots[['w7']],
             hbef_flux_plots[['w8']],
             hbef_flux_plots[['w9']],
             nrow = 3, ncol = 3)
