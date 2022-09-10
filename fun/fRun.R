RunSaveModels <- function(my_dir,
                          suv,
                          dat0,
                          n_iters,
                          n_warmup,
                          MConstant,
                          MContinuous)
{
  
  
  t0 <- Sys.time()
  my_dir0 <- paste0('res/', my_dir)
  
  
  
  dat <- filter(dat0, survey == suv) %>% arrange(age_mean_f) %>%
    mutate(birth_year = tsur - age_mean_f)
  
  
  
  mod_0   <- fFitModel(model= MConstant, dat,
                       m_name = 'Constant', n_iters = n_iters) #, n_warmup = n_warmup);   
  print(paste0(suv, ' finished ------ ConstantFOI'))
  
  mod_1   <- fFitModel(model= MContinuous, dat,
                       m_name ='ContinuousNormal', n_iters = n_iters) #, n_warmup = n_warmup);   
  print(paste0(suv,     ' finished ------ MNormal'))
  
  
  # ---------------------  
  age_max <- max(dat$age_mean_f)
  
  
  name_plot  <- paste0(my_dir0, '/plots/', suv, '.png')
  name_posterior  <- paste0(my_dir0, '/posterior/',suv, '.RDS')
  
  
  foi_mod <- rstan::extract(mod_1$fit, 'foi', inc_warmup = FALSE)[[1]]
  max_lambda <-  (as.numeric(quantile(foi_mod, 0.95))) * 1.3
  lambda_sim <- NA

  
  #  ---- Plotting
  PPC0    <- fCombinedPlots(res=mod_0, dat, lambda_sim, max_lambda)
  PPC1    <- fCombinedPlots(res=mod_1, dat, lambda_sim, max_lambda)
  
  
  mod_0$prev_expanded <- PPC0$prev_expanded
  mod_1$prev_expanded <- PPC1$prev_expanded
  
  # dif_m <- loo::compare (mod_0$loo_fit, mod_1$loo_fit) # loo::compare deprecated
  dif_m <- loo_compare (mod_0$loo_fit, mod_1$loo_fit)
  
  
  PPC0$summary_mod$difference <- 0; PPC0$summary_mod$diff_se <- 1; 
  PPC1$summary_mod$difference <- dif_m[1];   PPC1$summary_mod$diff_se <- dif_m[2];
  
  # browser()
  
  model_comparison <- rbind(PPC0$summary_mod, PPC1$summary_mod)
  
  
  res_comp <- compare_and_save_best_model( survey = suv, 
                                           model_comparison = model_comparison, 
                                           name_posterior, 
                                           mod_0, mod_1)
  
  model_comp    <- res_comp$model_comp
  mod_comp_plot <- get_model_comparison_plot(res_comp)
  
  
  
  parrange_0 <- vertical_plot_arrange_per_model(PPC0)
  parrange_1 <- vertical_plot_arrange_per_model(PPC1)
  

  png(name_plot, width = 200 * 7, height = 250 * 8)
  
  grid.arrange(parrange_0, 
               parrange_1,
               # mod_comp_plot, # SAD!! tuve que quitarlo porque me ocasionaba un error raro en el plot
               nrow = 1)
  dev.off()

  print(paste('end', suv))
  
  
  t5 <- Sys.time()
  time_taken <- t5-t0
  print(time_taken)
  
  
  
  res_survey <- list(dat = dat,
                     time_taken = time_taken,
                     mod_0 = mod_0, 
                     mod_1 = mod_1,
                     model_comp = model_comp,
                     PPC0 = PPC0,
                     PPC1 = PPC1)
  
  
  
  saveRDS(res_survey, name_posterior)
  
  
}


#Hello




