

library(tidyverse)
library(readxl)
library(epitrix)
rm(list=ls())

results_folder <- "col_2022_09_09_21_28_15"

codes <- read_excel("data/FOI_metadata_codes.xlsx") %>% arrange(DICTUM_Code)
codes_paper <- codes$DICTUM_Code[codes$used_burden_paper == 1] 


prev_plot_data <- data.frame()


# i <- "COL-035-18"
for (i in codes_paper) {
  paper_code_text_i <- paste0 ("ID:", codes$code_burden_paper[codes$DICTUM_Code==i])
  res_i <- readRDS(paste0("res/",results_folder,"/posterior/", i, ".RDS"))
  prev_plot_data_i <- res_i$PPC0$plots$plot_prev$data
  prev_plot_data_i$paper_code_text <- paste(paper_code_text_i)
  
  prev_plot_data_i$p_obs_l   <- prev_plot_data_i$pobs
  prev_plot_data_i$p_obs_u   <- prev_plot_data_i$pobslo
  prev_plot_data_i$p_obs_m   <- prev_plot_data_i$pobsup
  prev_plot_data_i$n_sample  <- prev_plot_data_i$sample_size

  
  if (is.data.frame(res_i$PPC1$prev_obs_binned_5)) {
    
    prev_plot_data_i$p_obs_l   <- NA
    prev_plot_data_i$p_obs_u   <- NA
    prev_plot_data_i$p_obs_m   <- NA
    prev_plot_data_i$n_sample  <- NA
    
    prev_plot_data_i$p_obs_m[prev_plot_data_i$age %in% res_i$PPC1$prev_obs_binned_5y$age] <-
      res_i$PPC1$prev_obs_binned_5y$pobs[res_i$PPC1$prev_obs_binned_5y$age %in% prev_plot_data_i$age]
    prev_plot_data_i$p_obs_l[prev_plot_data_i$age %in% res_i$PPC1$prev_obs_binned_5y$age] <-
      res_i$PPC1$prev_obs_binned_5y$pobslo[res_i$PPC1$prev_obs_binned_5y$age %in% prev_plot_data_i$age]
    prev_plot_data_i$p_obs_u[prev_plot_data_i$age %in% res_i$PPC1$prev_obs_binned_5y$age] <-
      res_i$PPC1$prev_obs_binned_5y$pobsup[res_i$PPC1$prev_obs_binned_5y$age %in% prev_plot_data_i$age]
    prev_plot_data_i$n_sample[prev_plot_data_i$age %in% res_i$PPC1$prev_obs_binned_5y$age] <-
      res_i$PPC1$prev_obs_binned_5y$sample_size[res_i$PPC1$prev_obs_binned_5y$age %in% prev_plot_data_i$age]

  } 
  
  prev_plot_data <- rbind(prev_plot_data, prev_plot_data_i)
  print(paste("finished ", i))
  rm(paper_code_text_i, res_i, prev_plot_data_i,  i)
  
}

prev_plot <-
  ggplot(prev_plot_data) +
  geom_ribbon(aes( x= age, ymin = plower, ymax = pupper), fill = '#c994c7') +
  geom_line(aes( x= age, y = medianv), colour = '#7a0177') +
  geom_errorbar(aes(age, ymin = p_obs_l, ymax = p_obs_u), width = 0.1) +
  geom_point(aes(age, p_obs_m, size = n_sample), fill = '#7a0177', colour = 'black', shape = 21) +
  theme_linedraw() +
  # theme_bw(strip.background = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid")) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
  theme(legend.position = 'none') +
  ylab ('Sero-positivity') + xlab("Age") + 
  facet_wrap(~survey)


path_plots <- paste0("res/", results_folder, "/")

png(paste0(path_plots, "summary_prev_plot.png"), width = 200 * 3, height = 250 * 4)
prev_plot
dev.off()


