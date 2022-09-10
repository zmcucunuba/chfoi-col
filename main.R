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
library(dplyr)
library(vscDebugger)
library(epitrix)
library(gsubfn) # obtain limits from applying cut function

source('fun/fUtility.R')
source('fun/fFitting.R')
source('fun/fRun.R')
source('fun/fCheck.R')
source('fun/fPlot.R')


# ---- Models 
# For these to run well, I need to re-compile the RDSs files 
# source('stanmodels-col/compile_stan_files_col_julia.R')

# ---- Models 
MConstant           <- readRDS('stanmodels-col/MConstant.RDS')
MContinuousNormal   <- readRDS('stanmodels-col/MContNormal.RDS')


dat0 <- readRDS("data/data-COL-2021.RDS")
(datasets <- as.character(unique(dat0$survey)))
(ld <- length(datasets))


# Automated name of the folder where results will be stored
my_dir <- epitrix::clean_labels(paste0('COL-', Sys.time()))
dir_results(my_dir) # from file XXX 
saveRDS(dat0, paste0('res/', my_dir, "/data-COL.RDS"))
summary_rep <-
  dat0 %>% dplyr::group_by(survey, tsur, country, test, setting, loc_type, 
                    ADM1, ADM2, ADM3, lat_dec, long_dec, source_type,
                    year_init, year_end, n_ages, sample_size) %>%  dplyr::summarise()

saveRDS(summary_rep, paste0('res/', my_dir, "/summary-COL.RDS"))



for (i in datasets) {
  # i = "COL-035-18" #Test survey
  print(paste('======================= reading', i))
  
  RunSaveModels(my_dir    = my_dir,
                suv       = i,
                dat0      = dat0,
                n_iters   = 3000,
                MConstant = MConstant,
                MContinuous   = MContinuousNormal
                )
  
  
}


options(error = browser) #Cuando uso este, VSC me deja ver la líne y script del error. THANK GOD!!

# xx <- readRDS("res/COL-2022-09-09 19:24:08/posterior/COL-001-01.RDS")

