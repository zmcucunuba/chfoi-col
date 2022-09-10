
# Obtaining R2
rm(list=ls())
library(tidyverse)
library(ggpmisc)

res_folder <- "res/COL-2022-09-08 00:51:17/"
post_folder <- paste0(res_folder, "posterior/")
datasets <- readRDS(paste0(res_folder, "summary-COL.RDS"))$survey


df <- data.frame()

for (i in datasets) {
  
  dd <- readRDS(paste0(post_folder, i, ".RDS"))
  df_i <- data.frame(y_obs = dd$mod_1$prev_expanded$pobs[!is.na(dd$mod_1$prev_expanded$pobs)],
                     y_pred = dd$mod_1$prev_expanded$medianv[!is.na(dd$mod_1$prev_expanded$pobs)],
                     survey = dd$dat$survey[1])
  
  df <- rbind(df, df_i)
  rm(dd)
  
}



png(filename = paste0(res_folder, "r2.png"), 
           width = 250 * 8, height = 250 * 8)

ggplot(data = df, aes(x = y_obs, y = y_pred)) +
  stat_poly_line(method = "lm") +
  stat_poly_eq() +
  facet_wrap(~ survey, scales = "free") +
  geom_point()

dev.off()






# 
# RSquare <- function(y_actual, y_predict){
#   cor(y_actual,y_predict)^2
# }
