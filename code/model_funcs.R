# library(ggplot2)
# library(ggthemes)
# library(ggsci)


get_index_pars <- function(chain, index){
  tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
  pars <- as.numeric(chain[index,2:(ncol(chain)-1)])
  names(pars) <- tmpNames
  return(pars)
}

export_theme <- theme_tufte() + 
  theme(
    text=element_text(family="sans"),
    
    ## Axis text and titles
    axis.text.x = element_text(size=8, family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=10,family="sans",vjust=-1),
    axis.title.y=element_text(size=10,family="sans"),
    
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    
    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    
    ## Legends
    legend.title=element_text(size=8,family="sans",face="italic"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=6,family="sans"),
    strip.background=element_rect(fill="#f0f0f0")
    )

ages_observed <- c(0,5,10,15,20,25,30,35,40)
desired_probs <- c(0,1,0.9,0.72,0.49,0.29,0.14,0.04,0)

viral_load_func_asymp <- function(pars, obs_t, convert_ct=TRUE, infection_time=0){
  tshift <- pars["tshift"] + infection_time ## Days post infection until growth
  desired_mode <- pars["desired_mode"] #+ infection_time## Days post growth to peak
  t_switch <- pars["t_switch"] ## Days post peak to switch point
  
  height <- pars["viral_peak"] ## Peak viral load
  obs_sd <- pars["obs_sd"] ## Variance about mean
  level_switch <- pars["level_switch"]
  
  true_0 <- pars["true_0"] ## y-axis shift
  yintercept <- pars["intercept"] ## Ct intercept
  
  wane_rate <- (height-level_switch)/t_switch
  wane_rate2 <- (level_switch - pars["LOD"])/pars["wane_rate2"]
  growth_rate <- (height-true_0)/desired_mode
  
  y <- numeric(length(obs_t))
  y[obs_t <= tshift] <- true_0
  y[obs_t > tshift & obs_t <= (desired_mode + tshift)] <- growth_rate * (obs_t[obs_t > tshift & obs_t <= (desired_mode + tshift)] - tshift) + true_0
  y[obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)] <- height - wane_rate * (obs_t[obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)] - (desired_mode+tshift))
  y[obs_t > (desired_mode + tshift + t_switch)] <- level_switch - wane_rate2*(obs_t[obs_t > (desired_mode + tshift + t_switch)] - (desired_mode + tshift + t_switch))
  ct <- y
  if(convert_ct){
    ct <- yintercept - log2(10) * (ct-pars["LOD"])
  }
  ct
}


viral_load_func <- function(pars, obs_t, convert_ct=TRUE){
  tshift <- pars["tshift"]
  desired_mode <- pars["desired_mode"]
  gamma_mode <- desired_mode - tshift
  
  gamma_sd <- pars["gamma_sd"]
  height <- pars["viral_peak"]
  obs_sd <- pars["obs_sd"]
  true_0 <- pars["true_0"]
  yintercept <- pars["intercept"]
  
  # Here are the corresponding rate and shape parameter values:
  ra <- ( gamma_mode + sqrt( gamma_mode^2 + 4*gamma_sd^2 ) ) / ( 2 * gamma_sd^2 )
  sh <- 1 + gamma_mode * ra
  scale_factor <- dgamma(gamma_mode, shape=sh,rate=ra)
  gamma_dist <- (height-true_0)*(dgamma(obs_t - tshift, shape=sh,rate=ra)/scale_factor) + true_0
  ct <- gamma_dist
  if(convert_ct){
    ct <- yintercept - log2(10) * (ct-pars["LOD"])
  }
  ct
}

viral_load_func_hinge <- function(pars, obs_t, convert_ct=TRUE){
  tshift <- pars["tshift"]
  desired_mode <- pars["desired_mode"]
  
  height <- pars["viral_peak"]
  obs_sd <- pars["obs_sd"]
  true_0 <- pars["true_0"]
  yintercept <- pars["intercept"]
  t_onset <- pars["incu"]
  tw <- pars["t_wane"]
  
  y <- numeric(length(obs_t))
  wane_rate <- height / (t_onset - desired_mode + tw)
  
  y[obs_t <= tshift] <- true_0
  y[obs_t > tshift & obs_t <= desired_mode] <- height * (obs_t[obs_t > tshift & obs_t <= desired_mode]-tshift) / (desired_mode-tshift)
  y[obs_t > desired_mode] <- height - wane_rate * (obs_t[obs_t > desired_mode] - desired_mode)
  
  ct <- y
  if(convert_ct){
    ct <- yintercept - log2(10) * (ct-pars["LOD"])
  }
  ct
}

## Function to give probability of observing x given age a and the viral kinetics curve
p_a <- function(x, a, pars,viral_loads,obs_model="normal",log_prob=FALSE, additional_detect_process=FALSE) {
  viral_load_sd <- pars["obs_sd"]
  LOD <- pars["intercept"]
  additional_probs <- rep(1, length(a))
  if(additional_detect_process){
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    ts <- a[a > t_switch] - t_switch
    additional_probs[which(a > t_switch)] <-  (1-pars["prob_detect"]*pars["t_unit"])^ts
  }
  ## If Cts are normally or Gumbelly distributed
  if(obs_model=="normal"){
    ## If a vector or single observation
    if(length(x) > 1){
      probs <- numeric(length(x))
      probs[x >= LOD] <- pnorm(LOD, mean=viral_loads, sd=viral_load_sd,lower.tail=FALSE,log.p=log_prob) + (1-additional_probs)*pnorm(LOD, mean=viral_loads, sd=viral_load_sd, lower.tail=TRUE, log.p=log_prob)
      probs[x < LOD] <- dnorm(x[x < LOD], mean=viral_loads, sd=viral_load_sd, log=log_prob)*additional_probs
    } else {
      if(x >= LOD) {
        probs <- pnorm(LOD, mean=viral_loads, sd=viral_load_sd,lower.tail=FALSE,log.p=log_prob) + (1-additional_probs)*pnorm(LOD, mean=viral_loads, sd=viral_load_sd, lower.tail=TRUE, log.p=log_prob)
      } else {
        probs <- dnorm(x, mean=viral_loads, sd=viral_load_sd, log=log_prob)*additional_probs
      }
    }
  } else {
    ## If a vector or single observation
    if(length(x) > 1){
      probs <- numeric(length(x))
      probs[x >= LOD] <- pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=FALSE,log.p=log_prob) + (1-additional_probs)*pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      probs[x < LOD] <- dgumbel(x[x < LOD], mu=viral_loads, sigma=viral_load_sd, log=log_prob)*additional_probs
    } else {
      if(x >= LOD) {
        probs <- pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=FALSE,log.p=log_prob) + (1-additional_probs)*pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      } else {
        probs <- dgumbel(x, mu=viral_loads, sigma=viral_load_sd, log=log_prob)*additional_probs
      }
    }
  }
  probs
}
## Function to give probability of observing x given age a and the viral kinetics curve
p_a1 <- function(x, a, pars,viral_loads,obs_model="normal",log_prob=FALSE) {
  viral_load_sd <- pars["obs_sd"]
  LOD <- pars["intercept"]
  
  ## If Cts are normally or Gumbelly distributed
  if(obs_model=="normal"){
    ## If a vector or single observation
    if(length(x) > 1){
      probs <- numeric(length(x))
      probs[x >= LOD] <- pnorm(LOD, mean=viral_loads, sd=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      probs[x < LOD] <- dnorm(x[x < LOD], mean=viral_loads, sd=viral_load_sd, log=log_prob)
    } else {
      if(x >= LOD) {
        probs <- pnorm(LOD, mean=viral_loads, sd=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      } else {
        probs <- dnorm(x, mean=viral_loads, sd=viral_load_sd, log=log_prob)
      }
    }
  } else {
    ## If a vector or single observation
    if(length(x) > 1){
      probs <- numeric(length(x))
      
      probs[x >= LOD] <- pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      probs[x < LOD] <- dgumbel(x[x < LOD], mu=viral_loads, sigma=viral_load_sd, log=log_prob)
      
    } else {
      if(x >= LOD) {
        probs <- pgumbel(LOD, mu=viral_loads, sigma=viral_load_sd,lower.tail=TRUE,log.p=log_prob)
      } else {
        probs <- dgumbel(x, mu=viral_loads, sigma=viral_load_sd, log=log_prob)
      }
    }
  }
  probs
}

prop_detectable <- function(a, pars,viral_loads, obs_model="normal", additional_detect_process=FALSE){
  viral_load_sd <- pars["obs_sd"]
  LOD <- pars["intercept"]
  additional_probs <- rep(1, length(a))
  ## Additional detection process
  if(additional_detect_process){
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    ts <- a[a > t_switch] - t_switch
    additional_probs[which(a > t_switch)] <-  (1-pars["prob_detect"]*pars["t_unit"])^ts
    #additional_probs[which(a < pars["tshift"])] <-  0
  }
  
  ## Are Cts normally or Gumbelly distributed?
  if(obs_model == "normal"){
    main_probs <- pnorm(LOD, mean=viral_loads, sd=viral_load_sd,lower.tail=TRUE,log.p=FALSE) 
  } else {
    main_probs <- extraDistr::pgumbel(LOD,mu=viral_loads,sigma=viral_load_sd, lower.tail=TRUE, log.p=FALSE)
  }
  main_probs * additional_probs
}

## Negative log likelihood function
ParamNegLogLik <- function(beta, ## Growth rate
                           viral_loads, ## Predicted viral loads for the age vector a.vec
                           p_a_use, ## Function for distribution of positive viral loads as function of age
                           prop_detectable_use, ## Function for proportion of individuals detectable as a function of age
                           a.vec, ## Age vector to solve model over
                           x.vec, ## Observed Cts
                           pars, ## Viral kinetics model parameters
                           obs_model="normal", ## How are Cts distributed about the mean?
                           additional_detect_process=FALSE ## Do we take into account the secondary loss process?
                           ) {
  PerBeta <- function(b) {
    ## Optional, but re-scaling e^-ba to not deal with stupidly large numbers. Cancels out anyway
    exp_tmp <- exp(-a.vec*b)
    exp_tmp <- exp_tmp/max(exp_tmp)
    
    log.phi <- log(sum(prop_detectable_use(a.vec, pars, viral_loads, obs_model, additional_detect_process)*exp_tmp))
    log.x <- sapply(X=x.vec, FUN=function(y) log(sum(p_a_use(y, a.vec, pars, viral_loads, obs_model, FALSE,additional_detect_process)*exp_tmp)))
    return(length(x.vec)*log.phi - sum(log.x))
  }
  return(sapply(X=beta, FUN=function(beta.val) PerBeta(beta.val)))
}


## Get predicted distribution of positive Cts
pred_dist <- function(beta,viral_loads, p_a_use, 
                      a.vec,x.vec, pars, obs_model="normal",additional_detect_process) {
  log.x <- sapply(X=x.vec, FUN=function(y) sum(p_a(y, a.vec, pars, viral_loads, obs_model,FALSE,additional_detect_process)*exp(-a.vec*beta)))
  return(log.x)
}

create_lik_func <- function(parTab, data, ages, PRIOR_FUNC=NULL,solve_ver="likelihood", test_cts=seq(0,39,by=1),
                           viral_load_func_use=viral_load_func_hinge, obs_model="normal", 
                           additional_detect_process=FALSE,solve_likelihood=TRUE,
                            ...) {
  par_names <- parTab$names
  f <- function(pars){
    names(pars) <- par_names
    #browser()
    vl <- viral_load_func_use(pars, ages)
    if(solve_ver == "likelihood"){
      if(solve_likelihood){
        lik <- ParamNegLogLik(pars["beta"], 
                              viral_loads=vl,
                              p_a_use=p_a, 
                              prop_detectable_use=prop_detectable, 
                              ages,
                              data,
                              pars,
                              obs_model, 
                              additional_detect_process)
      } else {
        likelihood <- 10000
      }
      prior <- 0
      if(!is.null(PRIOR_FUNC)){
        prior <- PRIOR_FUNC(pars)
      }
      return(-lik + prior)
    } else {
      res <- pred_dist(pars["beta"], vl, 
                       p_a=p_a,
                       ages, 
                       test_cts,
                       pars, obs_model,additional_detect_process)
      return(res)
    }
  }
  f
}


to.png <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  png(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}
to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

real_fitting_func <- function(ages, real_dat, parTab, filename, nchains,
                             chainwd, savewd, mcmcPars,mcmcPars_multivariate=NULL,
                             PRIOR_FUNC=NULL,nsamp=100,
                             all_priors=NULL, plot_par_lines=FALSE,
                             viral_load_func_use=viral_load_func_hinge,
                             obs_model="normal", additional_detect_process=FALSE,
                             traceplot=TRUE,solve_likelihood=TRUE) {
  pars_fitted <- parTab$values
  names(pars_fitted) <- parTab$names
  
  p_growth_hist <- ggplot() + 
    geom_histogram(data=data.frame(real_dat=real_dat),aes(x=real_dat),binwidth=1) +
    export_theme +
    xlab("Ct") +
    ylab("Count") +
    scale_x_continuous(trans="reverse") +
    scale_y_continuous(expand=c(0,0))
  
  to.pdf(print(p_growth_hist), paste0(savewd,"/",filename,"_hist.pdf"),height=3,width=5)
  
  ## RUN MCMC
  filenames <- paste0(filename,"_",1:nchains)
  chainwd1 <- paste0(chainwd, "/",filename)
  ##res <- foreach(j=seq_along(filenames),.packages = c("lazymcmc","rethinking","extraDistr","ggthemes")) %dopar% {

  for(j in seq_along(filenames)){
    #source("~/Google Drive/nCoV/viral_loads_time/code_20200616/model_funcs.R")
    
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
    #startTab <- parTab
    output <- run_MCMC(parTab=startTab, data=real_dat, ages=ages, mcmcPars=mcmcPars, filename=filenames[j], 
                       CREATE_POSTERIOR_FUNC=create_lik_func, mvrPars=NULL, 
                       PRIOR_FUNC = PRIOR_FUNC, OPT_TUNING=0.2,
                       viral_load_func_use=viral_load_func_use,
                       obs_model=obs_model, 
                       additional_detect_process=additional_detect_process, solve_likelihood=solve_likelihood)
    
    if(!is.null(mcmcPars_multivariate)){
      chain <- read.csv(output$file)
      bestPars <- get_best_pars(chain)
      chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
      covMat <- cov(chain)
      mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
      
      startTab <- parTab
      startTab$values <- bestPars
      output <- run_MCMC(parTab=startTab, data=real_dat, ages=ages, mcmcPars=mcmcPars_multivariate, filename=filenames[j], 
                         CREATE_POSTERIOR_FUNC=create_lik_func, mvrPars=mvrPars, 
                         PRIOR_FUNC = PRIOR_FUNC, OPT_TUNING=0.2,
                         viral_load_func_use=viral_load_func_use,
                         obs_model=obs_model, 
                         additional_detect_process=additional_detect_process, solve_likelihood=solve_likelihood)
    }
  }
  chains <- load_mcmc_chains(chainwd1,parTab = parTab,burnin=mcmcPars["adaptive_period"],multi=!is.null(mcmcPars_multivariate))
  if (traceplot) {
  par(ask=FALSE,mfrow=c(1,1))
  to.pdf(plot(chains[[1]]), paste0(savewd,"/",filename,"_traceplot.pdf"))
  par(ask=FALSE,mfrow=c(1,1))
  }
  
  chains_tmp <- load_mcmc_chains(chainwd1, parTab, FALSE, 1, mcmcPars["adaptive_period"], !is.null(mcmcPars_multivariate))
  chain <- as.data.frame(chains_tmp[[2]])
  p_left <- plot_dist_fit(chain, parTab, real_dat, ages,nsamp,viral_load_func_use,
                          obs_model=obs_model, additional_detect_process=additional_detect_process)
  
  estimated_pars <- parTab[parTab$fixed == 0,"names"]
  melted_chain <- reshape2::melt(chain,id.vars=c("sampno"))
  melted_chain$ver <- "Posterior"
  melted_chain <- melted_chain[,c("variable","value","ver")]
  melted_chain <- melted_chain[melted_chain$variable %in% estimated_pars,]
  
  all_plot <- melted_chain
  if(!is.null(all_priors)){
    all_plot <- bind_rows(melted_chain, all_priors)
  }
  
  parTab_plot <- parTab
  colnames(parTab_plot)[2] <- "variable"
  
  par_key <- c("obs_sd"="Variance", 
               "viral_peak"="Peak viral load",
               "t_wane"="Time to undetectable",
               "beta"="Growth rate",
               "true_0"="Y-axis shift")
  all_plot$variable <- as.character(all_plot$variable)
  #all_plot$var <- par_key[all_plot$variable]
  #parTab_plot$var <- par_key[parTab_plot$variable]
  
  all_plot$var <- all_plot$variable
  parTab_plot$var <- parTab_plot$variable
  
  parTab_plot[parTab_plot$variable == "beta","values"] <- NA
  
  p_right <- ggplot(all_plot) + 
    geom_density(aes(x=value,fill=ver),alpha=0.25) + 
    facet_wrap(~var,scales="free",ncol=2) +
    scale_fill_lancet() +
    export_theme +
    xlab("Value") +
    ylab("Density") +
    theme(legend.position=c(0.8,0.2))
  
  if(plot_par_lines){
    p_right <- p_right + geom_vline(data=parTab_plot[parTab_plot$fixed == 0,],aes(xintercept=values),linetype="dashed")
  }
  #to.pdf(plot(p_left), paste0(savewd,"/",filename,"_fit_dist.pdf"),width=7,height=5)
  #to.pdf(plot(p_right), paste0(savewd,"/",filename,"_posteriors.pdf"),width=7,height=5)
  
  to.pdf(plot(p_left), paste0(savewd,"/",filename,"_fit_dist.pdf"),width=7,height=5)
  to.pdf(plot(p_right), paste0(savewd,"/",filename,"_posteriors.pdf"),width=7,height=5)
  
  return(list(chain=chains[[2]],plot_dist=p_left,plot_posteriors=p_right))
}



plot_dist_fit <- function(chain, parTab, dat, ages, nsamp=100,
                         viral_load_func_use=viral_load_func_hinge,
                         obs_model="normal", additional_detect_process=FALSE){
  test_cts <- seq(0,parTab[parTab$names == "intercept","values"]-1,by=1)
  f <- create_lik_func(parTab, dat, ages, PRIOR_FUNC=NULL, solve_ver="model",test_cts=test_cts,
                      viral_load_func_use=viral_load_func_use,
                      obs_model=obs_model, additional_detect_process=additional_detect_process)
  samps <- sample(1:nrow(chain),nsamp)
  
  vls <- matrix(0, nrow=length(samps),ncol=length(test_cts))
  for(i in seq_along(samps)){
    pars <- get_index_pars(chain, samps[i])
    tmp <- f(pars)
    tmp <- tmp/sum(tmp)
    vls[i,] <- tmp
  }
  
  quants <- as.data.frame(t(apply(vls, 2, function(x) quantile(x, c(0.025,0.5,0.975)))))
  colnames(quants) <- c("lower","median","upper")
  quants$t <- test_cts
  
  p <- ggplot(data=tibble(x=dat)) + 
    geom_histogram(aes(x=x, y=..density..),fill="grey70",col="black",binwidth=1) +
    geom_ribbon(data=quants,aes(x=t+0.5,ymin=lower,ymax=upper),fill="blue",alpha=0.25) +
    geom_line(data=quants,aes(y=median,x=t+0.5)) +
    scale_y_continuous(expand=c(0,0)) +
    xlab("Ct value") +
    ylab("Density") +
    scale_x_continuous(breaks=seq(0,40,by=5)) +
    export_theme
  p
}

