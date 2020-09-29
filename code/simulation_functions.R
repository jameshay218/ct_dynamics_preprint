simulate_epidemic_process <- function(n_indivs, growth_rate, times){
  incidence <- exp(growth_rate*times)/n_indivs
  incidence <- incidence - incidence[1]
  incidence[length(incidence)] <- 0
  p_infected <- sum(incidence)
  print(paste0("Cumulative incidence to date: ", signif(p_infected,3)))
  p1 <- ggplot(tibble(time=times[1:(length(times)-1)],incidence=incidence[1:(length(times)-1)])) +
    ylab("Daily probability of infection") +
    xlab("Time") +
    theme_bw() +
    geom_line(aes(x=time,y=incidence))
  return(list(plot=p1, incidence=incidence, overall_prob_infection=p_infected))
}
seir_ode_virus <- function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  inc<-Y[5]
  pcr_exposed <- Y[6]
  pcr_pos <- Y[7]
  N <- sum(Y[1:4])
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  tg <- par[4] ## How long until detectable?
  omega <- par[5] ## How long until undetectable?
  
  dYdt<-vector(length=7)
  dYdt[1]= -beta*I*S/N 
  dYdt[2]= beta*I*S/N - sigma*E
  dYdt[3]= sigma*E - gamma*I
  dYdt[4]= gamma*I
  dYdt[5] = beta*I*S/N
  dYdt[6] = beta*I*S/N - tg*pcr_exposed
  dYdt[7] <- tg*pcr_exposed - omega*pcr_pos
  
  return(list(dYdt))
}

seir_ode <- function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  inc<-Y[5]
  N <- sum(Y[1:4])
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  dYdt<-vector(length=4)
  dYdt[1]= -beta*I*S/N 
  dYdt[2]= beta*I*S/N - sigma*E
  dYdt[3]= sigma*E - gamma*I
  dYdt[4]= gamma*I
  dYdt[5] = beta*I*S/N
  
  return(list(dYdt))
}

simulate_seir_process <- function(n_indivs, pars, times){
  # Set parameter values
  R0 <- pars["R0"]
  sigma<-pars["sigma"];
  gamma<-pars["gamma"];
  beta <- R0*gamma
  N <- n_indivs
  I0 <- pars["I0"]
  recovered0 <- pars["recovered0"]
  
  init <- c(N-I0-recovered0,0,I0,0,recovered0)
  t<-times
  par<-c(beta,sigma,gamma)
  # Solve system using lsoda
  sol<-lsoda(init,t,seir_ode,par)
  # Plot solution
  sol <- as.data.frame(sol)
  colnames(sol) <- c("time","S","E","I","R","cumulative_incidence")
  incidence <- diff(c(0,sol$cumulative_incidence/N))
  prevalence <- (sol$E + sol$I)/N
  sol <- reshape2::melt(sol, id.vars="time")
  sol$value <- sol$value/N
  
  p <- ggplot(sol) + 
    geom_line(aes(x=time,y=value,col=variable)) + 
    ylab("Per capita") + 
    xlab("Date") +
    theme_bw()
  p_inc <- ggplot(data.frame(x=t,y=incidence,y1=prevalence)) + geom_line(aes(x=x,y=y),col="red") +
    geom_line(aes(x=x,y=y1),col="blue") +
    ylab("Per capita incidence (red) and prevalence (blue)") + 
    xlab("Date") +
    theme_bw()
  
  return(list(plot=p, incidence_plot=p_inc, incidence=incidence, 
              seir_outputs=sol,prevalence=prevalence,
              overall_prob_infection=sum(incidence)))  
}

simulate_infection_times <- function(n, p_infected, incidence){
  scaled_incidence <- incidence/sum(incidence)
  are_infected <- numeric(n)
  infection_times <- numeric(n)
  for(i in 1:n){
    infection <- rbinom(1,1, p_infected)
    are_infected[i] <- infection
    if(infection == 1){
      t_inf <- sample(1:length(incidence), 1, prob=scaled_incidence)
      infection_times[i] <- t_inf
    } else {
      infection_times[i] <- -1
    }
  }
  return(infection_times)
}

simulate_symptom_onsets <- function(infection_times, incubation_period_par1=1.62, incubation_period_par2=0.418){
  onset_times <- numeric(length(infection_times))
  
  for(i in seq_along(onset_times)){
    if(infection_times[i] > 0){
      tmp_onset_time <- rlnorm(1, incubation_period_par1, incubation_period_par2)
      onset_times[i] <- infection_times[i] + tmp_onset_time
    } else {
      onset_times[i] <- -1
    }
  }
  return(onset_times)
}


simulate_viral_loads <- function(infection_times, times, kinetics_pars, ver="hinge",additional_detect_process=FALSE){
 viral_loads <- matrix(-100, nrow=length(infection_times), ncol=length(times))
  n <- length(infection_times)
  
  ## Pre-compute negative binomial draws for loss in detectability
  if(additional_detect_process){
    ## How many days do you remain detectable after hitting the switch point?
    days_still_detectable <- rnbinom(n, size=1,prob=kinetics_pars["prob_detect"])
  }
  
  for(i in seq_along(infection_times)){
    if(i %% 1000 == 0) print(i)
    if(infection_times[i] > 0){
      pars <- kinetics_pars
      mod_probs <- rep(1, length(times))
      if(ver == "hinge"){
        pars["incu"] <- pars["incu"] + infection_times[i]
        pars["tshift"] <- pars["tshift"] + infection_times[i]
        pars["desired_mode"] <- pars["desired_mode"] + infection_times[i]
        vl <- viral_load_func_hinge(pars, times, FALSE)
      } else if(ver == "gamma") {
        pars["desired_mode"] <- pars["desired_mode"] + infection_times[i]
        pars["tshift"] <- pars["tshift"] + infection_times[i]
        vl <- viral_load_func(pars, times, FALSE)
      } else {
        pars["tshift"] <- pars["tshift"] + infection_times[i]
        vl <- viral_load_func_asymp(pars, times, FALSE)
        
        ## Additional way that individuals can become undetectable
        if(additional_detect_process){
          ## How long to wait until undetectable?
          mod_probs[which(times >= pars["tshift"] + pars["desired_mode"] + pars["t_switch"] + days_still_detectable[i])] <- 0
        }
      }
      vl[vl < pars["true_0"]] <- pars["true_0"]
      vl[which(times < infection_times[i])] <- -100
      vl <- vl * mod_probs
      vl[which(mod_probs == 0)] <- -100
      viral_loads[i,] <- vl
    }
  }
  return(viral_loads=viral_loads)
}



sim_recover_func <- function(N, ages, truebeta, parTab, filename, nchains,
                             chainwd, savewd, mcmcPars,mcmcPars_multivariate=NULL,
                             PRIOR_FUNC=NULL) {
  pars_fitted <- parTab$values
  names(pars_fitted) <- parTab$names
  
  prob_infection <- exp(-1*truebeta*ages)/sum(exp(-1*truebeta*ages))
  
  ## Draw infection times for each person
  inf_times <- sample(ages, size=N, prob=prob_infection, replace=TRUE)
  
  ## Solve model, add noise and only use observable Cts
  sim_vls <- viral_load_func(pars_fitted, inf_times) + rnorm(N, 0, pars_fitted["obs_sd"])
  sim_vls <- sim_vls[sim_vls <  pars_fitted["intercept"]]
  sim_vls[sim_vls < 0] <- 0
  
  p_growth_hist <- ggplot() + 
    geom_histogram(data=data.frame(sim_vls=sim_vls),aes(x=sim_vls),binwidth=1) +
    export_theme +
    xlab("Ct") +
    ylab("Count") +
    scale_x_continuous(trans="reverse") +
    scale_y_continuous(expand=c(0,0))
  
  to.png(print(p_growth_hist), paste0(savewd,"/",filename,"_hist.png"),height=3,width=5,res=300,units="in")
  
  ## RUN MCMC
  filenames <- paste0(filename,"_",1:nchains)
  chainwd1 <- paste0(chainwd, "/",filename)
  
  
  #res <- foreach(j=seq_along(filenames),.packages = c("lazymcmc","rethinking","extraDistr","ggthemes")) %dopar% {
  #  source("~/Google Drive/nCoV/viral_loads_time/code_20200616/model_funcs.R")
  for(j in seq_along(filenames)){
    if(!file.exists(chainwd1)){
      dir.create(chainwd1)
    }
    setwd(chainwd1)
    ## The MCMC code uses the parameter table. Here, we should specify some random starting
    ## points in the "values" column.
    startTab <- parTab
    for(i in 1:nrow(startTab)){
      if(startTab$fixed[i] == 0) {
        startTab$values[i] <- runif(1, startTab$lower_start[i],startTab$upper_start[i])
      }
    }
    output <- run_MCMC(parTab=startTab, data=sim_vls, ages=ages, mcmcPars=mcmcPars, filename=filenames[j], 
                       CREATE_POSTERIOR_FUNC=create_lik_func, mvrPars=NULL, 
                       PRIOR_FUNC = PRIOR_FUNC, OPT_TUNING=0.2)
    
    if(!is.null(mcmcPars_multivariate)){
      chain <- read.csv(output$file)
      bestPars <- get_best_pars(chain)
      chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
      covMat <- cov(chain)
      mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
      
      startTab <- parTab
      startTab$values <- bestPars
      output <- run_MCMC(parTab=startTab, data=sim_vls, ages=ages, mcmcPars=mcmcPars_multivariate, filename=filenames[j], 
                         CREATE_POSTERIOR_FUNC=create_lik_func, mvrPars=mvrPars, 
                         PRIOR_FUNC = PRIOR_FUNC, OPT_TUNING=0.2)
    }
  }
  
  chains <- load_mcmc_chains(chainwd1,parTab = parTab,burnin=mcmcPars["adaptive_period"],multi=!is.null(mcmcPars_multivariate))
  par(ask=FALSE,mfrow=c(1,1))
  to.pdf(plot(chains[[1]]), paste0(savewd,"/",filename,"_traceplot.pdf"))
  par(ask=FALSE,mfrow=c(1,1))
  print(truebeta)
  return(res)
}


construct_beta <- function(arnaught, t_I, n_t) {
  beta_t_all <- arnaught / t_I
  if(length(arnaught) == 1) {
    function(t) beta_t_all
  } else {
    approxfun(0:n_t, beta_t_all)
  }
}

#' Simulate a discrete-time approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#' 
#' Ed Baskerville
#' 15 April 2020
#' 
#' No age structure.
#' 
#' @param n_t Number of units of time to simulate.
#' @param n_steps_per_t Number of discrete simulation steps to take per
#' unit of simulation time.
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
# # for debugging
# arnaught
# N = N
# E_init = 0          # Initial in E
# I_init = 1e-5 * N     # Initial in I
# t_E = 3               # time from exposed to infected
# t_I = 5               # time from infected to recovery
# n_t = 1000            # total simulation time
# max_R0 = 2.0          # Initial, max R0
# min_R0 = 0.9         # Minimum with intervention
# end_max_time = 50     # time of intervention
# start_min_time = 70  #time when R0 hits lowest value
# n_steps_per_t = 10
# CONTINUOUS = TRUE       # If TRUE, R decreases according to arctan after intervention. If FALSE, R


simulate_seir_stochastic <- function(
  arnaught, 
  t_E, 
  t_I,
  N, 
  S_init, 
  E_init, 
  I_init,
  n_t, 
  n_steps_per_t
) {
  # Precompute a few things
  delta_t <- 1 / n_steps_per_t
  
  # Draws a binomial based on a rate
  draw <- function(n, rate) {
    p <- 1 - exp(-rate * delta_t)
    rbinom(1, n, p)
  }
  
  # Function to compute beta at a particular time
  beta <- construct_beta(arnaught, t_I, n_t)
  
  # Step forward from t to t + delta_t
  step <- function(t, S_prev, E_prev, I_prev) {
    dS <- draw(S_prev, beta(t) * I_prev / N)
    dIR <- draw(I_prev, 1 / t_I)
    
    # SEIR model
    dEI <- draw(E_prev, 1 / t_E)
    list(
      S = S_prev - dS,
      E = E_prev + dS - dEI,
      I = I_prev + dEI - dIR,
      dS = dS,
      dEI = dEI,
      dIR = dIR
    )
  }
  
  # Set up state vectors over time
  S <- numeric(n_t + 1)
  S[1] <- S_init
  
  E <- numeric(n_t + 1)
  I <- numeric(n_t + 1)
  
  E[1] <- E_init
  I[1] <- I_init
  
  # Track transitions over time
  dS <- rep(NA, n_t + 1)
  dEI <- rep(NA, n_t + 1)
  dIR <- rep(NA, n_t + 1)
  
  # Simulate
  for(tt in 1:(n_t-1)) {
    S_prev <- S[tt]
    E_prev <- E[tt]
    I_prev <- I[tt]
    
    dS[tt+1] <- 0
    dEI[tt+1] <- 0
    dIR[tt+1] <- 0
    for(i in 1:n_steps_per_t) {
      state_next <- step(tt + delta_t * (i - 1), S_prev, E_prev, I_prev)
      S_prev <- state_next$S
      E_prev <- state_next$E
      I_prev <- state_next$I
      dS[tt+1] <- dS[tt+1] + state_next$dS
      dEI[tt+1] <- dEI[tt+1] + state_next$dEI
      dIR[tt+1] <- dIR[tt+1] + state_next$dIR
    }
    
    S[tt+1] <- S_prev
    E[tt+1] <- E_prev
    I[tt+1] <- I_prev
  }
  
  # Return each compartment over time
  data.frame(
    time = 0:n_t,
    S = S,
    E = E,
    I = I,
    R = N - S - E - I,
    dS = dS,
    dEI = dEI,
    dIR = dIR
  )
}
