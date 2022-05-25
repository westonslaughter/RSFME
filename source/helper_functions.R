# simple linear interpolation between con samples at frequency of q
# in:
# out: 
interpolate_con_to_q <- function(chem_df, q_df){
 library(imputeTS)
 # Join chem and q
 join_df <- full_join(chem_df, q_df, by = 'date') %>%
   arrange(date) %>%
   mutate(interp = 0)
 join_df$interp[is.na(join_df$con)] <- 1
 
 #interpolate chem to match q
 interp_df <- na_interpolation(join_df, option = "linear", maxgap = Inf)
 return(interp_df)
}
