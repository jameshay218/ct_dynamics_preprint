
likelihood_vl <- function(obs, predicted, sd, min_vl=0, limit_quantification=2){
  ## Change this if dealing with data with a different LOD and LOQ
  #normal_vl <- which(obs > min_vl)
  normal_vl <- which(obs > limit_quantification)
  unquantified_vl <- which(obs <= limit_quantification & obs > min_vl)
  low_vl <- which(obs <= min_vl)
  
  liks <- numeric(length(obs))
  liks[normal_vl] <- dnorm(obs[normal_vl],predicted[normal_vl],sd,TRUE)
  liks[unquantified_vl] <- log(pnorm(limit_quantification, predicted[unquantified_vl], sd, lower.tail=TRUE) - pnorm(min_vl, predicted[unquantified_vl], sd, lower.tail=TRUE))
  liks[low_vl] <- pnorm(min_vl, predicted[low_vl], sd, lower.tail=TRUE, TRUE)
  
  ## If want to fit Gumbel
  #liks[normal_vl] <- dgumbel(obs[normal_vl],predicted[normal_vl],sigma=sd,TRUE)
  #liks[low_vl] <- pgumbel(min_vl,predicted[low_vl],sigma=sd, lower.tail=TRUE, TRUE)
  
  liks
}



create_lik_func_vl <- function(parTab, data, PRIOR_FUNC=NULL, ver="model",model_ver="gamma",test_incu_period=FALSE,...) {
  par_names <- parTab$names
  ts <- data$t
  obs <- data$obs
  
  incu_pars <- which(parTab$names %like% "incu")
  
  
  f <- function(pars){
    if(test_incu_period){
      incu_dats <- tibble(i=seq_along(incu_pars),incu=pars[incu_pars])
      
      ts <- data %>% left_join(incu_dats,by="i") %>%
        mutate(t = t + incu) %>%
        pull(t)
    }
    names(pars) <- par_names
    if(model_ver == "gamma"){
      vl <- viral_load_func(pars, ts, FALSE)
    } else {
      vl <- viral_load_func_hinge(pars, ts, FALSE)
    }
    if(ver=="model"){
      return(vl)
    }
    prior <- 0
    if(!is.null(PRIOR_FUNC)){
      prior <- PRIOR_FUNC(pars)
    }
    lik <- likelihood_vl(obs, vl, pars["obs_sd"], pars["LOD"])*data$scale
    lik <- sum(lik)
    lik + prior
  }
  f
}