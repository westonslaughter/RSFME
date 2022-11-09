library(ggplot2)
library(dplyr)
library(tidyr)
library(feather)
library(macrosheds)
library(RColorBrewer)

# HBEF Flux Method Comparison
w3_flux <- read_feather('data/ms/hbef/stream_flux/w3.feather') %>%
  mutate(var = ms_drop_var_prefix(var)) %>%
  select(-ms_recommended) %>%
  distinct(wy, site_code, method, var, .keep_all = TRUE) %>%
  pivot_wider(names_from = var, values_from = val, id_cols = c('wy', 'site_code', 'method')) %>%
  filter(!is.na(Ca), !is.na(spCond))

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
                        " P =",signif(summary(fit)$coef[2,4], 5))) +
     theme(text = element_text(size = 26))
}

ggplotRegression(lm_fit)

# fluc methods sf
w3_flux_methods <- w3_flux %>%
  select(-Ca)

# bring in 'true' flux SpCond
w3_flux_true <- read_feather("data/ms/hbef/true/w3_sensor_wdisch.feather") %>%
  mutate(wy = water_year(date)) %>%
  group_by(wy) %>%
  summarise(spCond = sum(IS_spCond, na.rm = TRUE)/1000) %>%
  mutate(site_code = 'w3',
         method = 'true') %>%
  select(-spCond, spCond)

# bring in published flux SpCond
w3_flux_pub <- read.csv('data/raw/hbef_published_flux/ws3_stream_monthly_flux_gHa.csv') %>%
  mutate(date = paste0(Year, '-', Month, '-', '01'),
         wy = water_year(date)) %>%
  select(wy, SpecCond_volwt) %>%
  group_by(wy) %>%
  summarize(spCond = sum(SpecCond_volwt)) %>%
  filter(spCond > 0) %>%
  mutate(site_code = 'w3',
         method = 'published')%>%
  select(-spCond, spCond)

w3_spcond <- rbind(rbind(w3_flux_methods, w3_flux_pub), w3_flux_true)

# look at flux time series
fluxpal <- brewer.pal(n=8, name='Dark2')
ggplot(w3_spcond, aes(x = wy, y= spCond)) +
  geom_point(aes(color = method, size = method)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        ## legend.position="none",
        text = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold")) +
  scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'true', 'published', 'composite'),
                     values = fluxpal) +
  scale_size_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'true', 'published', 'composite'),
                     values = c(5, 5, 5, 5, 5, 5, 5, 5))
