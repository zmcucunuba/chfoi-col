
# --------- Get exposure matrix
get_exposure_matrix <- function(dat, yexpo)
{
  age_class <- dat$age_mean_f
  ly  <- length(yexpo)
  exposure       <- matrix(0,nrow=length(age_class), ncol=ly)   
  for (k in 1:length(age_class)) exposure[k,(ly-age_class[k]+1):ly] <-  1 
  exposure_output <- exposure
  return(exposure_output)
  
}



get_prev_expanded <- function(foi, dat) {
  
  ndat <- data.frame(age= 1:80)
  dim_foi <- dim(foi)[2] 
  if(dim_foi< 80)  {
    oldest_year <- 80 - dim_foi + 1
    foin <- matrix(NA, nrow = dim(foi)[1], 80)
    foin[, oldest_year: 80] <- foi
    foin[, 1:(oldest_year-1) ] <- rowMeans(foi[,1:5])
  } else {
    foin <- foi}
  
  foi_expanded <- foin
  
  
  age_class <- 1:NCOL(foi_expanded)
  ly  <- NCOL(foi_expanded)
  exposure       <- matrix(0,nrow=length(age_class), ncol=ly)   
  for (k in 1:length(age_class)) exposure[k,(ly-age_class[k]+1):ly] <-  1 
  exposure_expanded <- exposure
  
  
  iterf <- NROW(foi_expanded)
  age_max <- NROW(exposure_expanded)
  PrevPn <- matrix(NA, nrow = iterf, ncol = age_max)
  for (i in 1:iterf){
    PrevPn[i,] <- 1 - exp( - exposure_expanded %*% foi_expanded[i,])
  }
  
  lower <- apply(PrevPn, 2, function(x) quantile(x, 0.1))
  upper <- apply(PrevPn, 2, function(x) quantile(x, 0.9)) 
  medianv  <- apply(PrevPn, 2, function(x) quantile(x, 0.5)) 
  
  
  pred <- data.frame(age = 1:80, 
                     plower = lower, 
                     pupper = upper, 
                     medianv = medianv)
  
  dat_to_merge <- data.frame(age = dat$age_mean_f,
                             pobs = dat$prev_obs,
                             pobslo = dat$prev_obs_lower,
                             pobsup  = dat$prev_obs_upper,
                             sample_size = dat$total)
  
  pred_final <- merge(pred, dat_to_merge, by = 'age', all.x = TRUE)
  
  return(pred_final)
  
}




make_yexpo <- function(dat) {
  yexpo <- (seq_along(min(dat$birth_year):dat$tsur[1]))
}



get_posterior_summary <- function(results_chain)
{
  res <- sapply(results_chain, 
                function(i) c(quantile(i, c(0.5, 0.025, 0.975))))
  row.names(res) <- c('Median', 'Lower', 'Upper')
  return(res)
}


obtain_prevalence_extended <- function(dato, exposure, ly, nbreaks, lambdaYexpo) {
  
  ly            <- length(yexpo)
  
  if(nbreaks == 0 & ly < 100) 
  {
    ly = 100
    lambdaYexpo <- matrix(lambdaYexpo, nrow = length(lambdaYexpo), ncol = ly )
  }
  
  new_dat       <- data.frame(age = 1:ly) 
  exposure_new  <- matrix(0,nrow=length(new_dat$age), ncol=ly)   
  for (k in 1:length(new_dat$age)) exposure_new[k,(ly - new_dat$age[k]+1):ly] <-  1 
  
  olders <- data.frame(age = (ly+1):99) 
  exposure_olders <- matrix(1,nrow=length(olders$age), ncol=ly)   
  
  exposure_total <- rbind(exposure_new, exposure_olders)
  new_age_clases <- c(new_dat$age,   olders$age)
  
  iterf <- nrow(lambdaYexpo)
  PrevPn <- matrix(NA, nrow = iterf, ncol = length(new_age_clases))
  for (i in 1:iterf){
    PrevPn[i,] <- 1 - exp( - exposure_total %*% lambdaYexpo[i,])}
  new_PPP <- matrix(NA, ncol = 3, nrow = length(new_age_clases))
  for (j in seq_along(new_age_clases)){
    new_PPP[j,] <- quantile(PrevPn[,j], c(.025, .5, .975))}
  new_PPP        <- as.data.frame(new_PPP)
  new_PPP$age    <- new_age_clases
  names(new_PPP) <- c("L", "M", "U", "age")
  
  return(new_PPP)
  
}

thin_chain <- function(results_chain, thin = 10)
{
  results_chain <- results_chain[seq(1, nrow(results_chain), thin),]
  return (results_chain)
}


#  Streaming Stan Fit file
cleanObject <- function(object, pars){
  pars <- c(pars,'lp__')
  nn <- paste0('^',pars,'(\\[|$)',collapse="|")
  ids <-  grep(nn,  object@sim$fnames_oi)
  ids.2 <- which(names(object@par_dims) %in% pars)
  for(i in 1:4){
    a <- attributes(object@sim$samples[[i]])
    x <- object@sim$samples[[i]][ids]
    for(j in c('names','inits','mean_pars'))
      a[[j]] <- a[[j]][ids]
    attributes(x) <- a
    object@sim$samples[[i]] <- x
  }
  object@par_dims <- object@par_dims[ids.2]
  object@sim$dims_oi <-   object@sim$dims_oi[ids.2]  
  object@sim$pars_oi<- object@sim$pars_oi[ids.2]
  object@sim$fnames_oi <-  object@sim$fnames_oi[ids]
  object@sim$n_flatnames <- length(object@sim$fnames_oi)
  
  return(object)
  
}


sub_sample <- function(stanfit, n, keep_warmup = TRUE) {
  sim <- stanfit@sim
  samp <- sim$samples
  W <- sim$warmup
  I <- sim$iter
  sel <- c(if (keep_warmup) 1:W, sample((W + 1):I, size = n))
  subsamp <- lapply(samp, function(chain_samp) {
    lapply(chain_samp, function(x) x[sel])
  })
  stanfit@sim$samples <- subsamp
  stanfit
} 


dir_results <- function(name_dir)
  
  
{
  
  my_dir <- paste0('res/', name_dir)
  dir_plots <- paste0(my_dir, '/plots')
  dir_comp   <- paste0(my_dir, '/comp')
  dir_posterior  <- paste0(my_dir, '/posterior')
  
  if (dir.exists(my_dir) == FALSE) {
    dir.create(my_dir)
  } 
  
  if(dir.exists(dir_plots) == FALSE) {
    dir.create(dir_plots)
  } 
  
  
  if(dir.exists(dir_comp) == FALSE) {
    dir.create(dir_comp)
  } 
  
  if(dir.exists(dir_posterior) == FALSE) {
    dir.create(dir_posterior)
  } 
  
  
}
