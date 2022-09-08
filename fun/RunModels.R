RunSaveModels <- function(my_dir,
                          suv,
                          dat0,
                          n_iters,
                          n_warmup,
                          MConstant,
                          MTimeVarying)
{

  # my_dir    = my_dir
  # suv       = "COL-035-59"
  # dat0      = dat0
  # n_iters   = 1500
  # n_warmup = 400
  # MConstant = MConstant
  # MTimeVarying = MNormal
  
  # browser()
  
  t0 <- Sys.time()
  my_dir0 <- paste0('res/', my_dir)

  
  
  dat <- filter(dat0, survey == suv) %>% arrange(age_mean_f) %>%
    mutate(birth_year = tsur - age_mean_f)
  
  
  # fFitModel <- function(model, dat, m_name, 
  #                       n_iters = 3000, 
  #                       n_warmup = 1000,
  #                       n_thin = 2, 
  #                       delta = 0.90, 
  #                       mtreed = 10, 
  #                       Decades = 0){
  
  mod_0   <- fFitModel(model= MConstant, dat,
                      m_name = 'Constant', n_iters = n_iters) #, n_warmup = n_warmup);   
  print(paste0(suv, ' finished ------ ConstantFOI'))

  mod_1   <- fFitModel(model= MTimeVarying, dat,
                          m_name ='MNormal', n_iters = n_iters) #, n_warmup = n_warmup);   
  print(paste0(suv,     ' finished ------ MNormal'))
  
  
  # ---------------------  
  age_max <- max(dat$age_mean_f)
 

  name_plot  <- paste0(my_dir0, '/plots/', suv, '.png')
  name_posterior  <- paste0(my_dir0, '/posterior/',suv, '.RDS')
  # name_comp  <- paste0(my_dir0, '/comp/', suv, '.csv')


  foi_mod <- rstan::extract(mod_1$fit, 'foi', inc_warmup = FALSE)[[1]]
  max_lambda <-  (as.numeric(quantile(foi_mod, 0.95))) * 1.3
  lambda_sim <- NA
  

  
  #  ---- Plotting
  PPC1    <- fPCheck(mod_0, dat, lambda_sim, max_lambda)
  PPC3    <- fPCheck(mod_1, dat, lambda_sim, max_lambda)
 

  mod_0$prev_expanded <- PPC1$prev_expanded
  mod_1$prev_expanded <- PPC3$prev_expanded
 
  
  
  dif_m3 <- loo::compare (mod_0$loo_fit, mod_1$loo_fit)


  
  PPC1$lll$diff <- 0; PPC1$lll$diff_se <- 1; 
  PPC3$lll$diff <- dif_m3[1];   PPC3$lll$diff_se <- dif_m3[2]; 

  
  # Compare and save summary posterior of the best model
  model_comparison <- data.frame(rbind(PPC1$lll,
                                       PPC3$lll
  ))
  
  
  res_comp <- compare_and_save_best_model( survey = suv, 
                                           model_comparison = model_comparison, 
                                           # name_comp, 
                                           name_posterior, 
                                           mod_0, mod_1)
  
  model_comp    <- res_comp$model_comp
  mod_comp_plot <- get_model_comparison_plot(res_comp)
  
  dev.off()
  png(name_plot, width = 250 * 8, height = 250 * 8)
  grid.arrange(PPC1$plots, 
               PPC3$plots,
               mod_comp_plot,
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
                     model_comp = model_comp)
  
  
  
  saveRDS(res_survey, name_posterior)
  
  
}


#-----------------------


plot_info_table <- function(info){
  
  dato <- data.frame(
    y = NROW(info):1,
    text = paste0(rownames(info), ': ',info[,1])
  )
  p <- ggplot(dato, aes(x=1, y=y)) + 
    scale_y_continuous(limits=c(0, NROW(info) +1), breaks=NULL) +
    scale_x_continuous(breaks=NULL) + 
    theme_void() +
    geom_text(aes(label=text), size = 10, fontface = 'bold') 
  
  return(p)
}




#-------------------------

fPCheck <- function(res, dat, lambda_sim = NA, max_lambda) {
  
  

  
  fit <- res$fit
  
  if (is.character (res$fit) == FALSE)  {
    if  (class(fit@sim$samples)  != "NULL" ) {
      
      foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]
      model_name <- res$model
      #------- Loo estimates
      
      loo_fit <- res$loo_fit
      if (sum(is.na(loo_fit)) <1) 
      { 
        lll <- as.numeric((round(loo_fit$estimates[1,],2)))} else
        {
          lll <- c(-1e10, 0)
        }
      

      lll <- data.frame(model = res$model,
                        dataset = dat$survey[1],
                        c = dat$country[1],
                        tsur    = dat$tsur[1],
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
      
      
      scaleFUN <- function(x) sprintf("%.0f", x)
      X_sim <- rstan::extract(fit, 'P_sim')[[1]]
      lower <- apply(X_sim, 2, function(x) quantile(x, 0.1))
      upper <- apply(X_sim, 2, function(x) quantile(x, 0.9))
      medianv  <- apply(X_sim, 2, function(x) quantile(x, 0.5))
      pred_post <- data.frame(age = dat$age_mean_f,
                              plower = lower,
                              pupper = upper,
                              medianv = medianv,
                              pobs = dat$prev_obs,
                              sample_size = dat$total,
                              pobslo = dat$prev_obs_lower,
                              pobsup  = dat$prev_obs_upper)
      
      prev_expanded <- get_prev_expanded(foi, dat)
      prev_expanded$survey <- dat$survey[1]
      prev_plot2 <- ggplot(prev_expanded) +
        geom_ribbon(aes( x= age, ymin = plower, ymax = pupper), fill = 'pink') +
        geom_point(aes(age, pobs, size = sample_size), fill = 'white', colour = 'black') +
        geom_errorbar(aes(age, ymin = pobslo, ymax = pobsup), width = 0.1) +
        geom_line(aes(x = age, y = medianv), colour = 'darkred') +
        theme_bw(25) +
        coord_cartesian(xlim = c(0, 60), ylim = c(0,1)) +
        theme(legend.position = 'none') +
        ylab ('Sero-positivity') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        scale_x_continuous(labels=scaleFUN)
      
      
      
      
      
      foi_dat <- res$foi_cent_est
      
      #-------- This bit is to get the length of the foi data  
      
      if (!is.na(lambda_sim)) {
        lambda_mod_length <- NROW(foi_dat) 
        lambda_sim_length <- length(lambda_sim)
        
        if (lambda_mod_length < lambda_sim_length) {
          remove_x_values <- lambda_sim_length - lambda_mod_length
          lambda_sim <- lambda_sim[-c(1:remove_x_values)]
        }
        
        foi_dat$simulated <- lambda_sim
        
      }
      
      
      #--------  
      foi_dat$medianv[1] <- NA
      foi_dat$lower[1] <- NA
      foi_dat$upper[1] <- NA
      
      lambda_plot <- ggplot(foi_dat) +
        geom_ribbon(aes( x= year, ymin = lower, ymax = upper), fill = '#2ca25f', alpha = 0.5) +
        geom_line(aes(x = year, y = medianv), colour = 'darkgreen', size = 1.5) +
        theme_bw(25) +
        coord_cartesian(ylim = c(0, max_lambda)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        scale_x_continuous(labels=scaleFUN) +
        ylab ('FOI') 
      
      
      if (!is.na(lambda_sim)) {
        lambda_plot <= lambda_plot + geom_line(aes(x = year, y = simulated), colour = 'red', size = 1.5) 
      }
      
      
      
      
      rhats <- rhat(res$fit, "foi")
      if(any(is.nan(rhats))) {
        rhats[which(is.nan(rhats))] <- 0
        
      }
      
      res_rhats <- data.frame(year = res$RealYexpo, rhat = rhats)
      res_rhats$rhat[res_rhats$rhat == 0] <- NA # This is because I'm not estimating these foi values
      
      phats <- ggplot(res_rhats, aes(year, rhat)) +
        geom_line(colour = 'purple') +
        geom_point() +
        coord_cartesian(ylim = c(0.7, 2)) +
        geom_hline(yintercept = 1.1, colour = 'blue', size = 2) +
        theme_bw(25) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        scale_x_continuous(labels=scaleFUN) +
        ylab ('R^ hat') 
      
      if (any(rhats > 1.1 ) == TRUE) { 
        lll$converged = 'No' 
      } 
      
      if (any(rhats > 1.1 ) == FALSE) { 
        lll$converged = 'Yes' 
      } 
      
      
      
      t_lll <- t(lll)
      data_plot <- plot_info_table(t_lll) 
      
      
      plots <- grid.arrange(data_plot, 
                            prev_plot2, 
                            lambda_plot, 
                            phats, nrow = 4,
                            heights = c(1.5, 1, 1, 1))
      res_p <- list (plots = plots,
                     lll     = lll,
                     prev_expanded = prev_expanded)
      
      
      
    } } else
      
    {
      print ('model did not run')
      print_warning <- 'errors'
      df <- data.frame() 
      
      g0 <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x = 4, y = 5, label = print_warning) +
        theme_bw(25) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
        ylab('') + xlab('') 
      g1 <- g0
      g0 <- g0 +   labs(subtitle = res$model)+
        theme(plot.title = element_text(size=10))
      
      plots <- grid.arrange(g0, g1, g1, g1, g1, nrow = 5)
      res_p <- list (plots = plots,
                     loo_fit = c('Not available'),
                     lll     = c('Not available'),
                     prev_expanded = c('Not available')
      )
      
    }
  
  
  
  
  
  return(res_p)
  
  
}



#========================= PLOT comparison
get_model_comparison_plot <- function(res_comp) {
  model_comp  <- res_comp$model_comp
  best_model  <- as.character(res_comp$best_model_data1$model)
  best_modelP <- as.numeric(res_comp$best_model_data1$pvalue)
  
  emptyp <- ggplot(data = data.frame()) +
    geom_point() +
    xlim(0,1) + ylim (0,1) + theme_void() +
    annotate('text',  x = .5, y = .6, label = best_model, size = 17) +
    annotate('text',  x = .5, y = .55, label = '(best model)', size = 15) 
  
  
  infot <- filter(model_comp, converged == 'Yes') %>% select(model, diff, diff_se, pvalue, best) %>%
    mutate(diff = round(diff, 2), 
           diff_se = round(diff_se, 2),
           pvalue = round(pvalue, 4))
  
  blank <- data.frame(x= 1:10, y = 1:100)
  table_pars <- 
    ggplot(blank, aes(x, y)) +
    geom_blank() + ylab('') + xlab ('') + 
    annotation_custom(tableGrob(d= infot,
                                theme = ttheme_default(base_size = 20))) +
    theme_void(30)
  
  pf <- plot_grid(emptyp, table_pars, nrow = 2 )
  
  return(pf)
  
}

