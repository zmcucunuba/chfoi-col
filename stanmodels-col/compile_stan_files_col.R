# compile_stan_files 

rm(list=ls())
library(rstan)
library(tidyverse)


MConstant     <- stan_model('stanmodels-col/constant_foi_Bi.stan')
saveRDS(MConstant, 'stanmodels-col/MConstant.RDS')

MContNormal    <- stan_model('stanmodels-col/continuous_foi_normal_Bi.stan')
saveRDS(MContNormal, 'stanmodels-col/MContNormal.RDS')

