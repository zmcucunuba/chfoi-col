

compare_and_save_best_model <- function (survey, 
                                         model_comparison, 
                                         # name_comp, 
                                         name_posterior,
                                         mod_0, mod_1)
{
  
  # browser() # Aquí debo verificar que model_comparison viene con todo, inlcuido performance
  
  # --------------- Model comparison 
  model_comp <- model_comparison %>% select(-performance)
  model_comp <- model_comparison %>% mutate (pvalue = 1- pnorm(difference/diff_se,0,1))
  model_comp <- model_comparison
  model_comp$pvalue = 1- pnorm(model_comp$difference/model_comp$diff_se,0,1) 
  
  model_comp$pvalue[is.nan(model_comp$pvalue)] <- 1
  model_comp$pvalue <- model_comp$pvalue * runif(NROW(model_comp), min = 1, max = 1.0001)# I make this just to ensure I get different values
  
  model_comp$better <- NA
  model_comp$better[model_comp$difference > 0] <- 'Yes'
  model_comp$better[model_comp$difference <= 0] <-'No'
  
  model_comp$better[model_comp$model == 'Constant'] <- "-"
  model_comp$pvalue[model_comp$model == 'Constant'] <- 0
  
  # browser() #Ya aparece el elpd!!!!
  model_comp$converged[model_comp$elpd == -1.000e+10] <- 'No'
  
  
  ds_one <- filter(model_comp, converged == 'Yes')
  print(paste0('number of converged models = ', NROW(ds_one)))
  
  # browser() #  ERROR Error in gList(...) : only 'grobs' allowed in "gList"
  
  elps_order <-  rev(sort(ds_one$elpd))[1:2]
  best <- filter(model_comp, elpd %in% elps_order) %>% arrange(-elpd)# This is to make sure I keep only three
  best_model1 <- as.character(best$model[1])
  best_model2 <- as.character(best$model[2])
  
  model_comp$best <- NA
  model_comp$best[model_comp$model == best_model1] <- 1
  model_comp$best[model_comp$model == best_model2] <- 2
  model_comp <- model_comp %>% arrange(best)
  model_comp$pvalue <- round(model_comp$pvalue, 6)
  
  
  # write.csv(model_comp, name_comp)
  
  
  # --------------- Best model
  best_model_data1 <- filter(model_comp, best ==1) # Here I choose the maximun difference rather than the lowest p value
  best_model_data2 <- filter(model_comp, best ==2) # Here I choose the maximun difference rather than the lowest p value
  
  RealYexpo  <-  mod_0$RealYexpo
  best_model_1 <- as.character(best_model_data1$model)
  best_model_2 <- as.character(best_model_data2$model)
  
  # browser() # No estaba mostrando la convergencia
  
  if(best_model_1 == mod_0$model) {
    res_file_1 <- mod_0
  }
  
  
  
  if(best_model_1 ==  mod_1$model) {
    res_file_1 <- mod_1
  }
  
  
  
  
  
  # ------- Best 2
  if(best_model_2  == mod_0$model) {
    res_file_2 <- mod_0
  }
  
  
  
  
  if(best_model_2 ==  mod_1$model) {
    res_file_2 <- mod_1
  }
  
  
  
  
  
  
  # --- save_best_model
  extract_and_save(res_file_1, res_file_2, 
                   best_model_1, best_model_2,
                   name_posterior, 
                   survey, RealYexpo)
  
  
  
  res_comp <- list (best_model_data1 = best_model_data1,
                    best_model_data2 = best_model_data2,
                    model_comp = model_comp)
  return(res_comp)
  
}






extract_and_save <- function(res_file_1, res_file_2, 
                             best_model_1, best_model_2,
                             name_file, 
                             survey, RealYexpo)
  
  
{
  
  
  foi_0 <- rstan::extract(res_file_1$fit, 'foi', inc_warmup = FALSE)[[1]]
  foi_1 <- rstan::extract(res_file_2$fit, 'foi', inc_warmup = FALSE)[[1]]
  
  
  foi_cent_est1 <- data.frame(year  = RealYexpo,
                              lower = apply(foi_0, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi_0, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi_0, 2, function(x) quantile(x, 0.5))) %>%
    mutate(best= 'best1', name_model = best_model_1)
  
  foi_cent_est2 <- data.frame(year  = RealYexpo,
                              lower = apply(foi_1, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi_1, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi_1, 2, function(x) quantile(x, 0.5))) %>%
    mutate(best= 'best2', name_model = best_model_2)
  
  
  
  foi_cent_est <- rbind(foi_cent_est1, foi_cent_est2)
  
  
  foi_0_post_1000s <- dplyr::sample_n(as.data.frame(foi_0), size = 1000) %>% mutate(best = 'best1', name_model = best_model_1)
  foi_1_post_1000s <- dplyr::sample_n(as.data.frame(foi_1), size = 1000) %>% mutate(best = 'best2', name_model = best_model_2)
  
  foi_post_1000s <- rbind(foi_0_post_1000s, foi_1_post_1000s)
  
  colnames(foi_post_1000s)[1: length(RealYexpo)] <- RealYexpo
  
  
  prev1 <- res_file_1$prev_expanded %>% mutate(best = 'best1', name_model = best_model_1)
  prev2 <- res_file_2$prev_expanded %>% mutate(best = 'best2', name_model = best_model_2)
  
  
  prevalence_expanded <- rbind(prev1, prev2)
  
  
  fres <- list(dataset = survey,
               foi_cent_est = foi_cent_est,
               foi_post_1000s = foi_post_1000s,
               prevalence    = prevalence_expanded)
  
  saveRDS(fres, name_file)
  
}


extract_summary_mod <- function (res, dat){
  
  model_name <- res$model
  #------- Loo estimates
  
  loo_fit <- res$loo_fit
  if (sum(is.na(loo_fit)) <1) 
  { 
    lll <- as.numeric((round(loo_fit$estimates[1,],2)))} else
    {
      lll <- c(-1e10, 0)
    }
  
  # browser() # Okay! Aquí obtengo bien lll con elpd
  
  summary_mod <- data.frame(model = res$model,
                            dataset = dat$survey[1],
                            country = dat$country[1],
                            year    = dat$tsur[1],
                            test    = dat$test[1], 
                            antibody = dat$antibody[1],
                            n_sample = sum(dat$total),
                            n_agec  = length(dat$age_mean_f),
                            n_iter  = res$n_iters,
                            performance = "_____",
                            elpd = lll[1],
                            se = lll[2],
                            converged = NA
  )
  
  rhats <- get_table_rhats (res)
  if (any(rhats$rhat > 1.1 ) == FALSE) { 
    summary_mod$converged = 'Yes'  } 
  
  
  return(summary_mod)
}


get_table_rhats <- function(res) {
  
  rhats <- rhat(res$fit, "foi")
  
  if(any(is.nan(rhats))) {
    rhats[which(is.nan(rhats))] <- 0}
  
  res_rhats <- data.frame(year = res$RealYexpo, rhat = rhats)
  res_rhats$rhat[res_rhats$rhat == 0] <- NA # This is because I'm not estimating these foi values
  
  return(res_rhats)
  
}


