library(ggplot2)
library(dplyr)

# HBEF Flux Method Comparison
w3_flux <- read_feather('data/ms/hbef/stream_flux/w3.feather') %>%
  mutate(var = ms_drop_var_prefix(var)) %>%
  select(-ms_recommended) %>%
  pivot_wider(names_from = var, values_from = val, id_cols = c('wy', 'site_code', 'method'))

# create multiple linear model
lm_fit <- lm(Ca ~ spCond, data=w3_flux)
summary(lm_fit)

ggplotRegression <- function (fit) {
   ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
     geom_point() +
     stat_smooth(method = "lm", col = "red") +
     labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                        "Intercept =",signif(fit$coef[[1]],5 ),
                        " Slope =",signif(fit$coef[[2]], 5),
                        " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(lm_fit)

library(RColorBrewer)
fluxpal <- brewer.pal(n=6, name='Dark2')

# look at flux time series
ggplot(w3_flux, aes(x = wy, y= spCond)) +
  geom_point(aes(color = method, size = method)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        ## legend.position="none",
        text = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold")) +
  scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite'),
                     values = fluxpal) +
  scale_size_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite'),
                     values = c(5, 5, 5, 5, 5, 5, 5, 5))
