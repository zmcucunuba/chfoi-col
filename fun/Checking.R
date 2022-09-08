

compare_and_save_best_model <- function (survey, 
                                         model_comparison, 
                                         name_comp, name_posterior,
                                         res1, res3)
{
  

  # --------------- Model comparison 
  model_comp <- model_comparison %>% select(-performance)
  model_comp <- model_comp %>% mutate (pvalue = 1- pnorm(diff/diff_se,0,1)) 
  model_comp$pvalue[is.nan(model_comp$pvalue)] <- 1
  model_comp$pvalue <- model_comp$pvalue * runif(NROW(model_comp), min = 1, max = 1.0001)# I make this just to ensure I get different values
  
  model_comp$better <- NA
  model_comp$better[model_comp$diff > 0] <- 'Yes'
  model_comp$better[model_comp$diff <= 0] <-'No'
  
  model_comp$better[model_comp$model == 'Constant'] <- "-"
  model_comp$pvalue[model_comp$model == 'Constant'] <- 0
  model_comp <- data.frame(model_comp)
  
  model_comp$converged[model_comp$elpd == -1.000e+10] <- 'No'
  model_comp$elpd <- model_comp$elpd 
  
  ds_one <- filter(model_comp, converged == 'Yes')
  print(paste0('number of converged models = ', NROW(ds_one)))
  
  elps_order <-  rev(sort(ds_one$elpd))[1:2]
  best <- filter(model_comp, elpd %in% elps_order) %>% arrange(-elpd)# This is to make sure I keep only three
  best_model1 <- as.character(best$model[1])
  best_model2 <- as.character(best$model[2])
    
  model_comp$best <- NA
  model_comp$best[model_comp$model == best_model1] <- 1
  model_comp$best[model_comp$model == best_model2] <- 2
  model_comp <- model_comp %>% arrange(best)
  model_comp$pvalue <- round(model_comp$pvalue, 6)
  

  write.csv(model_comp, name_comp)


  # --------------- Best model
  best_model_data1 <- filter(model_comp, best ==1) # Here I choose the maximun difference rather than the lowest p value
  best_model_data2 <- filter(model_comp, best ==2) # Here I choose the maximun difference rather than the lowest p value

  RealYexpo  <-  res1$RealYexpo
  best_model_1 <- as.character(best_model_data1$model)
  best_model_2 <- as.character(best_model_data2$model)
  

  if(best_model_1 == res1$model) {
    res_file_1 <- res1
  }
  

  
  if(best_model_1 ==  res3$model) {
    res_file_1 <- res3
  }
  
  
  
  

  # ------- Best 2
  if(best_model_2  == res1$model) {
    res_file_2 <- res1
  }
  
  
  
  
  if(best_model_2 ==  res3$model) {
    res_file_2 <- res3
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
  

  foi1 <- rstan::extract(res_file_1$fit, 'foi', inc_warmup = FALSE)[[1]]
  foi2 <- rstan::extract(res_file_2$fit, 'foi', inc_warmup = FALSE)[[1]]

 
  foi_cent_est1 <- data.frame(year  = RealYexpo,
                              lower = apply(foi1, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi1, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi1, 2, function(x) quantile(x, 0.5))) %>%
    mutate(best= 'best1', name_model = best_model_1)
  
  foi_cent_est2 <- data.frame(year  = RealYexpo,
                              lower = apply(foi2, 2, function(x) quantile(x, 0.1)),
                              upper = apply(foi2, 2, function(x) quantile(x, 0.9)),
                              median = apply(foi2, 2, function(x) quantile(x, 0.5))) %>%
    mutate(best= 'best2', name_model = best_model_2)
  
  
  
  foi_cent_est <- rbind(foi_cent_est1, foi_cent_est2)
  

  foi1_post_1000s <- dplyr::sample_n(as.data.frame(foi1), size = 1000) %>% mutate(best = 'best1', name_model = best_model_1)
  foi2_post_1000s <- dplyr::sample_n(as.data.frame(foi2), size = 1000) %>% mutate(best = 'best2', name_model = best_model_2)

  foi_post_1000s <- rbind(foi1_post_1000s, foi2_post_1000s)
 
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

