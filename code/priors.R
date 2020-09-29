prior_func_beta <- function(pars){
  prior1 <- dnorm(pars["beta"], 0, 0.2,log=TRUE)
  obs_sd_prior <- dnorm(pars["obs_sd"], 7, 1,log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], 7, 1,log=TRUE)
  t_wane_prior <- dnorm(pars["t_wane"], 27, 3,log=TRUE)
  true_0_prior <- dnorm(pars["true_0"], -7.5, 1,log=TRUE)
  prior1 + obs_sd_prior + viral_peak_prior + t_wane_prior + true_0_prior 
}


generate_all_priors <- function(prior_table){
  all_priors <- NULL
  for(i in 1:nrow(prior_table)){
      if(prior_table$prior_ver[i] == "normal"){
      all_priors[[i]] <- tibble(variable=prior_table$par[i], 
                                value=rnorm(10000, mean=prior_table$mean[i],sd=prior_table$sd[i]),
                                ver="Prior")
      } else {
        beta1_mean <- prior_table$mean[i]
        beta1_sd <- prior_table$sd[i]
        beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
        beta_beta <- beta_alpha*(1/beta1_mean - 1)
        all_priors[[i]] <- tibble(variable=prior_table$par[i], 
                                  value=rbeta(10000, beta_alpha,beta_beta),
                                  ver="Prior")
        
      }
  }
  all_priors <- do.call("bind_rows", all_priors)
  ## variable, value, ver
  all_priors
}

#obs_sd_prior <- 0
#beta1_mean <- 0.1
#beta1_sd <- 1
#beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
#beta_beta <- beta_alpha*(1/beta1_mean - 1)
#wane_2_prior <- dbeta(pars["wane_rate2"],beta_alpha,beta_beta,log=TRUE)

prior_func_hinge2 <- function(pars){
  prior1 <- dnorm(pars["beta"], 0, 0.5,log=TRUE)
  obs_sd_prior <- dnorm(pars["obs_sd"], 6, 2,log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], 9, 3,log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],30,5,log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],16,3,log=TRUE)
  level_prior <- dnorm(pars["level_switch"],3.7,0.5,log=TRUE)
  beta1_mean <- 0.1
  beta1_sd <- 0.1
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  #beta_alpha <- 1
  #beta_beta <- 1
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  #beta_prior <- dnorm(pars["prob_detect"], 0.1, 0.1, log=TRUE)
  prior1 + obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior + level_prior + beta_prior
}

prior_func_hinge3 <- function(pars){
  #prior1 <- dnorm(pars["beta"], 0, 0.5,log=TRUE)
  obs_sd_prior <- dnorm(pars["obs_sd"], 6, 2,log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], 9, 3,log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],30,5,log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],16,3,log=TRUE)
  level_prior <- dnorm(pars["level_switch"],3.7,0.5,log=TRUE)
  beta1_mean <- 0.1
  beta1_sd <- 0.1
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  #beta_alpha <- 1
  #beta_beta <- 1
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  #beta_prior <- dnorm(pars["prob_detect"], 0.1, 0.1, log=TRUE)
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior + level_prior + beta_prior
}
prior_func <- function(pars){
  viral_loads <- viral_load_func(pars, ages_observed)
  obs <- prop_detectable(ages_observed,pars,viral_loads)
  sum(dnorm(obs, desired_probs, sd=0.1,log=TRUE)) + dnorm(pars["beta"],0,0.5,log=TRUE)
}

prior_func_wane <- function(pars){
  wane <- pars["t_wane"]
  sum(dnorm(wane, 21, sd=0.1,log=TRUE))
}