rm(list=ls())
library(rstan)
library(tidyverse)
library(reshape2)
library(bayesplot)
library(loo)
library(pracma)
library(cowplot)
library(grid)
library(gridExtra)
library(Hmisc)

source('fun/Utility.R')
source('fun/Fitting.R')
source('fun/RunModels.R')
source('fun/Checking.R')


# ---- Models 
# For these to run well, I need to re-compile the RDSs files 
# source('stanmodels-col/compile_stan_files_col_julia.R')

# ---- Models 
MConstant       <- readRDS('stanmodels-col/MConstant.RDS')
MNormal        <- readRDS('stanmodels-col/MContNormal.RDS')



#------ Data
dat0 <- readRDS('data/clean_data_total_models_COL.RDS')
only_one_age_class <- dat0 %>% group_by(survey) %>% summarise(age_classes = n()) %>%
  filter(age_classes <2)
only_one_age_class <- only_one_age_class$survey
not_working <- NA
not_working_or_1_age_class <- c(only_one_age_class, not_working)

dat0 <- filter(dat0, !survey %in% not_working_or_1_age_class)

dat0$counts <- round(dat0$counts, 0)
dat0$total <- round(dat0$total, 0)
dat0 <- dat0 %>% filter(country == "COL")
(datasets <- as.character(unique(dat0$survey)))
(ld <- length(datasets))

my_dir <- 'COL-2021'
dir_results(my_dir) # from file XXX 

saveRDS(dat0, paste0('res/', my_dir, "/data-", my_dir,  '.RDS'))

summary <-
  dat0 %>% group_by(survey, tsur, country, test, setting, loc_type, 
                    ADM1, ADM2, ADM3, lat_dec, long_dec, source_type,
                    year_init, year_end, n_ages, sample_size) %>%  summarise()

saveRDS(summary, paste0('res/', my_dir, "/summary-", my_dir,  '.RDS'))




for (i in datasets) {
  
  print(paste('======================= reading', i))
  
  RunSaveModels(my_dir    = my_dir,
                suv       = i,
                dat0      = dat0,
                n_iters   = 2200,
                # n_warmup = 300,
                MConstant = MConstant,
                MTimeVarying   = MNormal,
                
  )
  
  
}

dd <- readRDS("res/COL-2021/posterior/COL-001-01.RDS")
dd$model_comp

