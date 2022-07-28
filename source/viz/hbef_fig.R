# MacroSheds
library(devtools)
#install_github("https://github.com/MacroSHEDS/macrosheds.git")
library(macrosheds)
library(here)
library(RiverLoad)

source('streamlined/source/usgs_helpers.R')
source('source/helper_functions.R')

#declare ms core dir
my_ms_dir <- here('streamlined/data/ms')

# load in a df of MS sites with
ms <- ms_download_site_data() %>%
    select(site_code, ws_area_ha, latitude, longitude) %>%
    rename(lat = latitude,
           long = longitude)

# get Q for RBI calc
ms_sites_ls <- unique(ms_sites$site_code)
raw_data_q <- ms_load_product(
    my_ms_dir,
    prodname = "discharge",
    site_codes = ms_sites_ls,
    sort_result = TRUE,
    warn = F
)

#acutal_ms_site_ls <- 
actual_files <- list.files('streamlined/data/ms', full.names = T, recursive = T)
actual_files_ls <- tibble(file = actual_files) %>%
    filter(grepl('stream_chemistry', file),
           !grepl('documentation', file)) %>%
    mutate(site_file = str_split_fixed(file, '/', n = Inf)[,6],
           site_code =  str_split_fixed(site_file, '[.]', n = Inf)[,1]) %>%
    filter(site_code %in% unique(raw_data_q$site_code)) %>%
    pull(site_code)

##### for ws2 #####
# for each MS sites
s <- actual_files_ls[53] #58
ms_dir <- my_ms_dir

# load in site area
area <- ms[ms$site_code == s,]$ws_area_ha

# load in site Q
raw_data_q <- ms_load_product(
    ms_dir,
    prodname = "discharge",
    site_codes = s,
    sort_result = TRUE,
    warn = TRUE
)

# load in stream chemistry
raw_data_chem <- ms_load_product(
    ms_dir,
    prodname = "stream_chemistry",
    site_codes = s,
    sort_result = TRUE,
    warn = TRUE
)

# join data

ms_q <- raw_data_q %>%
    #filter(ms_status == 0) %>%
    pivot_wider(id_cols = c(datetime, site_code),
                names_from = var,
                values_from = val)

ms_chem <- raw_data_chem %>%
    filter(ms_status == 0,
           grepl('NO3', var)) %>%
    pivot_wider(id_cols = c(datetime, site_code),
                names_from = var,
                values_from = val)

ms_comb <- full_join(ms_q, ms_chem, by = c('site_code', 'datetime')) %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

cq_by_wy <- ms_comb %>%
    mutate(q_log = log10(IS_discharge),
           c_log = log10(GN_NO3_N)) %>%
    filter(is.finite(c_log),
           is.finite(q_log)) %>%
    group_by(wy) %>%
    summarize(n = length(datetime),
              cq_r2 = summary(lm(c_log~q_log, singular.ok = T))$r.squared) %>%
    mutate(wy = as.integer(as.character(wy)))

ggplot(cq_by_wy, aes(x = wy, y = cq_r2)) +
    geom_point()+
    geom_line()
    
out_frame <- tibble(wy = as.integer(),
                    flux = as.numeric(),
                    site_code = as.character(),
                    method = as.character())
### calc different flux methods
for(y in unique(ms_comb$wy)) {
    ms_df_wy <- ms_comb %>%
        filter(wy == y) %>%
        rename(date = datetime,
               nitrate_n_mgl = GN_NO3_N,
               q_lps = IS_discharge)
    
    chem_df <- ms_df_wy %>%
        ungroup()%>%
        select(date, con = nitrate_n_mgl, site_code, wy) %>%
        na.omit()
    
    q_df <- ms_df_wy %>%
        ungroup() %>%
        select(date, q_lps, site_code, wy)
    
    chem_check <- nrow(ms_df_wy[!is.na(ms_df_wy$nitrate_n_mgl),]) > 4
    q_check <- nrow(q_df) > 311
    align_check <- ms_df_wy %>%
        select(nitrate_n_mgl, q_lps) %>%
        na.omit() %>%
        nrow(.) > 2
    com_check <- chem_check*q_check*align_check
    
    if(com_check == 0 ) {
        print(paste(s, y, 'insufficient data, skipped'))
        next
    } else {
        
        pw_wy <- calculate_pw(chem_df = chem_df, q_df = q_df)
        
        beale_wy <- calculate_beale(chem_df = chem_df, q_df = q_df)
        
        rating_wy <- calculate_rating(chem_df = chem_df, q_df = q_df)
        
        composite_wy <- calculate_composite_from_rating_filled_df(
            generate_residual_corrected_con(chem_df = chem_df, q_df = q_df)
        )
        
        append_frame <- tibble(wy = y,
                            flux = c(pw_wy[[1]], beale_wy[[1]], 
                                     rating_wy [[1]], composite_wy$flux[[1]],
                                     NA), 
                            site_code = s,
                            method = c('pw', 'beale', 'rating', 'composite', 'mean'))
        append_frame$flux[5] <- append_frame %>%
            filter(method != 'mean') %>%
            summarize(flux = mean(flux))%>%
            .$flux
        
        out_frame <- rbind(out_frame, append_frame)
        
    }
    
}

w1_dif_by_wy <- out_frame %>%
    mutate(wy = as.integer(as.character(wy))) %>%
    pivot_wider(names_from = 'method', values_from = 'flux') %>%
    na.omit() %>%
    mutate(pw = pw-mean,
           beale = beale-mean,
           rating = rating-mean,
           composite = composite-mean) %>%
    select(-mean) %>%
    pivot_longer(names_to = 'method', values_to = 'load',
                 cols = -c(wy, site_code)) %>%
    left_join(., cq_by_wy, by = 'wy')

w1_load_by_wy <- out_frame %>%
    mutate(wy = as.integer(as.character(wy))) %>%
    filter(wy != 2021,
           method != 'mean') %>%
    pivot_wider(names_from = method,
                values_from = flux) %>%
    left_join(., cq_by_wy, by = 'wy') %>%
    mutate(best_value = ifelse(cq_r2 > 0.3, composite, pw),
           best_label = ifelse(cq_r2 > 0.3, 'composite', 'pw')) %>%
    pivot_longer(names_to = 'method', values_to = 'load',
                 cols = c(pw, beale, rating, composite, best_value))



ggplot(filter(w1_load_by_wy, method %in% c('pw', 'composite')), aes(x = as.integer(as.character(wy)),
                          y = load, fill = method))+
    geom_col()+
    theme_few()
    facet_wrap(~method,
               ncol = 1)
