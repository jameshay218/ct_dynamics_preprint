## Versions ode, odin, baskerville
simulate_seir_wrapper <- function(population_n, solve_times, pars, version="ode",tstep=1){
  if(version == "ode"){
    ####################################################
    ## Deterministic model
    ####################################################
    seir_pars <- c("R0"=pars["R0"],"gamma"=pars["gamma"],"sigma"=pars["sigma"],
                   "I0"=pars["I0"]*population_n,"recovered0"=0)
    names(seir_pars) <- c("R0","gamma","sigma","I0","recovered0")
    epidemic_process <- simulate_seir_process(population_n,seir_pars,solve_times)
    res <- epidemic_process$seir_outputs %>% pivot_wider(names_from="variable",values_from="value")
    res$Rt <- (res$S) * pars["R0"]
    res$inc <- c(0,diff(res$cumulative_incidence))
    res$inc <- round(res$inc*population_n)
    res <- res %>% rename(step=time)
    incidence <- epidemic_process$incidence
    overall_prob <- epidemic_process$overall_prob_infection
  } else if (version == "odin"){
    ####################################################
    ## Stochastic model
    ####################################################
    beta1 <- pars["R0"]*pars["gamma"]
    gamma1 <- pars["gamma"]
    sigma1 <- pars["sigma"]
    I0 <- ceiling(pars["I0"]*population_n)
    
    seir <- seir_generator(beta=beta1,sigma=sigma1,gamma=gamma1,S_ini=population_n-I0,I_ini=I0)
    res <- seir$run(solve_times)
    
    ## Make sure we get a simulation with an outbreak
    while(max(res[,"I"]) <= I0) res <- seir$run(solve_times)
    incidence <- res[,"inc"]/population_n
    overall_prob <- max(res[,"R"])/population_n
    res <- as.data.frame(res)
    res$Rt <- (res$S/population_n) * pars["R0"]
  } else {
    I0 <- ceiling(pars["I0"]*population_n)
    res <- simulate_seir_stochastic(pars["R0"],1/pars["gamma"],1/pars["sigma"],
                                    population_n,
                                    population_n-I0,0,I0,length(solve_times),tstep)
    res <- res[1:length(solve_times),]
    incidence <- res[,"dS"]/population_n
    overall_prob <- max(res[,"R"])/population_n
    res <- as.data.frame(res)
    res$step <- res$time
    res$Rt <- (res$S/population_n) * pars["R0"]
  }
  incidence <- incidence * population_n
  res_melted <- res %>% pivot_longer(-step)
  p_compartments <- res_melted %>% filter(name %in% c("S","E","I","R","cumulative_incidence")) %>%
    ggplot() + geom_line(aes(x=step,y=value,col=name)) +
    ylab("Per capita") + 
    xlab("Date") +
    theme_bw() +
    theme(legend.position="top")
  
  p_inc <- ggplot(data.frame(x=solve_times,y=incidence)) + 
    geom_line(aes(x=x,y=y),col="red") +
    ylab("True incidence") + 
    xlab("Date") +
    theme_bw()
  p_rt <- res_melted %>% filter(name == "Rt") %>%
    ggplot() +
    geom_line(aes(x=step,y=value),col="blue") +
    scale_y_continuous(limits=c(0,pars["R0"]+1)) +
    geom_hline(yintercept=1,linetype="dashed") +
    ylab("Rt") +
    xlab("Date") +
    theme_bw()
  inc_plot <- p_compartments / p_inc / p_rt
  
  ## Get growth rates
  GR_daily <- log(incidence[2:length(incidence)]/incidence[1:(length(incidence)-1)])
  GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
  GR_daily_dat <- data.frame(t=solve_times[2:length(solve_times)],GR=GR_daily,ver="daily")
  
  lastdays <- c(seq(10,50,by=10),35)
  
  GR_all <- GR_daily_dat
  for(t in seq_along(lastdays)){
    GR_full <- NULL
    lastday <- lastdays[t]
    for (i in (lastday+1):length(solve_times)) {
      #for (i in 1:length(GR_daily)) {
      end_index <- i-1
      start_index <- max(1, (i-lastday))
      GR_full <- c(GR_full,mean(GR_daily[start_index:end_index], na.rm=TRUE))
    }
    GR_full_dat <- data.frame(t=(lastday+1):length(solve_times), GR=GR_full,ver=as.character(lastday))
    GR_all <- bind_rows(GR_all, GR_full_dat)
  }
  gr_crossover <- GR_all %>% filter(ver == "daily") %>% 
    filter(t < 250 & t > 100) %>%
    mutate(abs_gr = abs(GR)) %>% 
    filter(abs_gr == min(abs_gr, na.rm=TRUE)) %>% pull(t)
  
  p_gr <- ggplot(GR_all %>% filter(ver != "daily")) +
    geom_line(aes(x=t,y=GR,col=ver)) +
    #geom_point(aes(x=t,y=GR,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    theme_bw()
    
  p_gr1 <- ggplot(GR_all %>% filter(ver == "daily")) +
    geom_line(aes(x=t,y=GR),col="black") +
    #geom_point(aes(x=t,y=GR,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    theme_bw()
  
  
  list(seir_outputs=res, incidence=incidence,overall_prob=overall_prob,plot=inc_plot, growth_rates=GR_all,
       growth_rate_p =p_gr1/p_gr)
}


simulate_observations_wrapper <- function(
  incidence, times, symp_frac=0.35, 
  population_n=length(incidence),
  incu_period_par1=1.621,incu_period_par2=0.418,
  conf_delay_par1=5,conf_delay_par2=2){
  
  not_infected <- population_n-sum(incidence)
  inc_dat <- tibble(infection_time=c(NA,times),inc=c(not_infected,incidence))
  inc_dat <- inc_dat %>% uncount(inc)
  
  inc_dat <- inc_dat %>% 
    mutate(i=1:n()) %>%
    mutate(is_infected = ifelse(is.na(infection_time), 0, 1)) %>% 
    ## Is this individual going to be symptomatic?
    mutate(is_symp=ifelse(is_infected, rbinom(n(), 1, symp_frac), 0)) %>%
    ## Symptom onset time
    mutate(incu_period=ifelse(is_infected & is_symp, rlnorm(n(), incu_period_par1, incu_period_par2), NA),
           onset_time=infection_time+round(incu_period)) %>%
    ## Confirmation time
    mutate(confirmation_delay=extraDistr::rdgamma(n(),conf_delay_par1,conf_delay_par2))
  inc_dat
}


## Test an increasing fraction of people per day
logistic_func <- function(t,start_prob,end_prob, growth_rate, switch_point=100){
  start_prob + (end_prob-start_prob)/(1 + exp(-growth_rate*(t-switch_point)))
}

simulate_reporting <- function(individuals,
                               frac_report=1,
                               timevarying_prob=NULL, 
                               solve_times, 
                               symptomatic=FALSE,
                               contact_tracing=FALSE, 
                               contact_tracing_penalty=2){
  
  ## The basic form is that a constant fraction of individuals are reported each day
  n_indivs <- nrow(individuals)
  sampled_individuals <- individuals
  
  if(contact_tracing & symptomatic){
    print("Note, with `contact_tracing` turned on, surveillance will be based on symptomatic individuals regardless of the 
          `symptomatic` flag")
  }
  
  ## Flat reporting rate
  if(is.null(timevarying_prob)){
    ## Base case, just sampled individuals at random
    if(!symptomatic){
      sampled_individuals <- sampled_individuals %>% 
        sample_frac(frac_report) %>%
        group_by(i) %>%
        mutate(sampled_time=sample(solve_times,n()),
               ## We sample but then have to wait for a result
               confirmed_time=sampled_time+confirmation_delay) %>%
        ungroup()
    } else {
    ## Symptomatic based surveillance. Observe individuals based on symptom onset date
      sampled_individuals <- sampled_individuals %>% 
        filter(is_symp==1) %>% 
        sample_frac(frac_report) %>%
        mutate(sampled_time=onset_time+confirmation_delay,
               ## We sample individuals some number of days after symptom onset
               confirmed_time=sampled_time)
    }
    ## Reporting rate varies over time
  } else {
    if(!symptomatic){
      ## How many individuals are we going to observe by the end?
      frac_report_overall <- 1-prod(1-timevarying_prob$prob)
      
      ## We will first get the fraction of individuals we will sample over the whole period
      ## Then, we will change the time-varying reporting probability to the relative number
      ## sampled on each day
      scaled_timevarying_prob <- timevarying_prob$prob/frac_report_overall
      
      sampled_individuals <- sampled_individuals %>%
        sample_frac(frac_report_overall) %>%
        group_by(i) %>%
        mutate(sampled_time = sample(timevarying_prob$t, n(), replace=FALSE, prob=scaled_timevarying_prob),
               confirmed_time = sampled_time + confirmation_delay) %>%
        ungroup()
      } else {
        ## This is quite different - if you have symptom onset on day t, there is a probability that you will be observed
        ## Symptomatic based surveillance. Observe individuals based on symptom onset date
        sampled_individuals <- sampled_individuals %>% 
          filter(is_symp==1) %>% 
          left_join(timevarying_prob %>% rename(onset_time=t) %>% dplyr::select(-ver), by="onset_time") %>%
          group_by(i) %>%
          mutate(is_reported=rbinom(n(), 1, prob)) %>%
          filter(is_reported == 1) %>%
          mutate(sampled_time=onset_time+confirmation_delay,
                 ## We sample individuals some number of days after symptom onset
                 confirmed_time=sampled_time)
      }
  }
  
    
  ## Plot incidence of infections,  onsets, confirmations and number sampled per day
  grouped_dat <- sampled_individuals %>% 
    dplyr::select(i, infection_time, onset_time, sampled_time, confirmed_time) %>%
    pivot_longer(-i) %>% 
    drop_na() %>%
    group_by(name, value) %>% 
    tally() %>%
    rename(var=name,
           t=value,
           n=n) %>%
    complete(var, nesting(t),fill = list(n = 0)) %>%
    mutate(ver="Sampled individuals")
  
  grouped_dat_all <- individuals %>% 
    dplyr::select(i, infection_time, onset_time) %>%
    pivot_longer(-i) %>% 
    drop_na() %>%
    group_by(name, value) %>% 
    tally() %>%
    rename(var=name,
           t=value,
           n=n) %>%
    complete(var, nesting(t),fill = list(n = 0)) %>%
    mutate(ver="All individuals")
  
  grouped_dat_combined <- bind_rows(grouped_dat, grouped_dat_all)
  
  p_all <- grouped_dat_combined %>% ggplot() +
    geom_line(aes(x=t,y=n,col=var)) +
    theme_bw() +
    ylab("Number of individuals") +
    xlab("Time") +
    theme(legend.position="bottom") +
    facet_wrap(~ver,ncol=1, scales="free_y")
  
  ## Get day-by-day growth rate
  grouped_dat_combined <- grouped_dat_combined %>% group_by(var, ver) %>% 
    mutate(gr_daily=log(n/lag(n,1))) %>%
    mutate(gr_daily = ifelse(is.finite(gr_daily), gr_daily, NA)) %>%
    ungroup()
  
  ## Get rolling average growth rate
  growth_rates_all <- expand_grid(grouped_dat_combined, window=seq(10,50,by=10)) %>% 
    arrange(var, ver, window, t) %>%
    group_by(var, ver, window) %>%
    mutate(gr_window=zoo::rollmean(gr_daily, window,align="right",fill=NA)) %>%
    mutate(window = as.factor(window))
  
  p_gr <- ggplot(growth_rates_all) +
    geom_line(aes(x=t,y=gr_window,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    
    theme_bw() +
    theme(legend.position="bottom") +
    facet_grid(window~var)
  
  
  list(sampled_individuals=sampled_individuals,plot=p_all, plot_gr=p_gr)
  
}


simulate_viral_loads_wrapper <- function(indiv_dat,solve_times=0:365,
                             kinetics_pars,viral_load_func_ver="hinge2",obs_ver="gumbel",
                             additional_detect_process=TRUE){
    vl_dat <- indiv_dat %>%
    group_by(i) %>% 
    mutate(last_detectable_day = infection_time + rnbinom(n(), 1, prob=kinetics_pars["prob_detect"]) + 
             kinetics_pars["tshift"] + kinetics_pars["desired_mode"] + kinetics_pars["t_switch"],
           vl=viral_load_func_asymp(kinetics_pars, sampled_time, FALSE, infection_time)) %>%
    mutate(vl = ifelse(sampled_time > last_detectable_day, -100, vl),
           vl = ifelse(is_infected == 0, -100, vl)) %>%
    mutate(ct = kinetics_pars["intercept"] - log2(10)*(vl - kinetics_pars["LOD"]),
           ct_obs = rgumbel(n(), ct, kinetics_pars["obs_sd"])) %>%
    mutate(ct_obs = pmin(ct_obs, kinetics_pars["intercept"]))
    
    return(viral_loads=vl_dat)
}
