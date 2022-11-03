library(ggplot2)
library(dplyr)
library(tidyr)
library(feather)
library(macrosheds)
library(RColorBrewer)

source('source/helper_functions.R')
source('source/usgs_helpers.R')

# HBEF Flux Method Comparison
w6_flux <- read_feather('data/ms/hbef/stream_flux/w6.feather') %>%
  mutate(var = ms_drop_var_prefix(var)) %>%
  select(-ms_recommended) %>%
  distinct(wy, site_code, method, var, .keep_all = TRUE) %>%
  pivot_wider(names_from = var, values_from = val, id_cols = c('wy', 'site_code', 'method')) %>%
  filter(!is.na(Ca))

# create multiple linear model

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

## lm_fit <- lm(Ca ~ spCond, data=w6_flux)
## summary(lm_fit)
## ggplotRegression(lm_fit)

# fluc methods sf
w6_flux_methods <- w6_flux %>%
  group_by(wy, method, site_code) %>%
  filter(site_code == 'w6') %>%
  summarize(Ca = sum(Ca))

# bring in published flux SpCond
w6_flux_pub <- read.csv('data/raw/hbef_published_flux/ws6_stream_monthly_flux_gHa.csv') %>%
  mutate(date = paste0(Year, '-', Month, '-', '01'),
         wy = water_year(date)) %>%
  group_by(wy) %>%
  # NOTE: this unit conversion right?
  summarize(Ca = sum(Ca_flux)) %>%
  ## filter( > 0) %>%
  mutate(site_code = 'w6',
         method = 'published',
         # and convert from grams to kg
         Ca = (Ca)/1000
         )

w6_spcond <- bind_rows(w6_flux_methods, w6_flux_pub)

# look at flux time series
fluxpal <- rev(brewer.pal(n=7, name='Dark2'))
# trim outliers
# create detect outlier function (if x is 1.5x the interquartile range higher than 3rd Quart or lower than 1st quart)
detect_outlier <- function(x) {
    # calculate first quantile
    Quantile1 <- quantile(x, probs=.25)
    # calculate third quantile
    Quantile3 <- quantile(x, probs=.75)
    # calculate inter quartile range
    IQR = Quantile3-Quantile1
    # return true or false
    x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}
w6 <- w6_spcond[!detect_outlier(w6_spcond$Ca),]
w6_nopub <- w6[w6$method != 'published',]
w6_pub <- w6[w6$method == 'published',]

transparent <- fluxpal
transparent[length(transparent)] <- '#FFFFFF'

ggplot() +
  geom_point(w6 %>% filter(method == 'published'), mapping = aes(x=wy, y=Ca, color = method, size = method, shape = method, fill = method)) +
  geom_line(w6 %>% filter(method != 'published'), mapping = aes(x=wy, y=Ca, group = method, color = method)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        ## legend.position="none",
        text = element_text(size = 26),
        plot.title = element_text(size = 24, face = "bold")) +
  scale_shape_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite', 'published'),
                     values = c(16, 16, 16, 16, 16, 16, 24)) +
  scale_color_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                     values = fluxpal) +
  scale_fill_manual(breaks = c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                     values = fluxpal) +
  scale_size_manual(breaks =  c('average', 'pw', 'beale', 'rating', 'wrtds', 'composite',  'published'),
                    values = c(2, 2, 2, 2, 2, 2, 6)) +
  ggtitle("Hubbard Brook, Watershed 6, 1960's-2020,\nCalcium Flux Estimates by Various Methods ") +
  ylab('Annual Cumulative Calcium Flux (kg/ha)\n ') +
  xlab('\n Water Year (Oct-Sep)') +
  scale_x_discrete(breaks=seq(1960, 2020, 10))
