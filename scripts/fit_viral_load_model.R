## Main script 1: fit kinetics model based on expected parameters
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(lazymcmc)
library(patchwork)
library(lhs)
library(optimx)

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop must move to the correct working directory to source the model functions
#main_wd <- "~/Documents/Harvard/Research/Infectious Diseases/COVID-19/PCR Cts/Code & Results/ct_inference"
main_wd <- "~/Documents/GitHub/ct_dynamics_preprint/"
setwd(main_wd)

source("code/model_funcs.R")

######################################################
## 1. FIT GAMMA MODEL SUBJECTIVELY
######################################################
parTab_gamma <- read.csv("pars/partab_gamma fit.csv")

## Takes 2 days for viral loads to peak
parTab_gamma[parTab_gamma$names == "tshift","values"] <- 2
## Viral loads peak 4 days post infection
parTab_gamma[parTab_gamma$names == "desired_mode","values"] <- 4
## Dealing with log10 per ul
parTab_gamma[parTab_gamma$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
ages_observed <- c(-5, 0, 5, 10, 15, 20, 25, 30, 35,40,45,50) + assumed_incu_period
desired_probs <- c(0,1,0.9,0.72,0.49,0.29,0.14,0.04,0,0,0,0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_gamma <- parTab_gamma$values
names(pars_gamma) <- parTab_gamma$names

ages <- ages_observed

######################################################
## COST FUNCTION
cost_function_gamma <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1){
  pars1 <- pars_gamma
  ## Change true 0, gamma variance, observation variance and peak viral load
  pars1["true_0"] <- explore_pars[1]
  pars1["gamma_sd"] <- explore_pars[2]
  pars1["obs_sd"] <- explore_pars[3]
  pars1["viral_peak"] <- explore_pars[4]
  peak <- viral_load_func(pars1,pars1["desired_mode"],convert_ct=TRUE)
  cost1 <- 0
  if(use_ct_cost){
    cost1 <- (11.4 - qnorm(0.01,peak,pars1["obs_sd"]))^2
  }
  vls <- viral_load_func(pars1, ages, convert_ct=TRUE)
  obs <- prop_detectable(ages,pars1, vls)
  costs2 <- (obs-desired_probs)^2
  
  sum(costs2 + cost1*ct_cost_weight)
}

######################################################
## A) RUN OPTIM
######################################################
fit_gamma <- optim(c(-5,20,7,7), cost_function_gamma,control=list(maxit=10000,abstol=1e-8,reltol=1e-8),
             use_ct_cost=TRUE, ct_cost_weight=1)

pars1 <- pars_gamma
pars1["true_0"] <- fit_gamma$par[1]
pars1["gamma_sd"] <- fit_gamma$par[2]
pars1["obs_sd"] <- fit_gamma$par[3]
pars1["viral_peak"] <- fit_gamma$par[4]

parTab_gamma$values <- pars1
write_csv(parTab_gamma,"pars/partab_gamma_fitted.csv")

######################################################
## B) PLOT AGAINST PROPORTION DETECTABLE
######################################################
test_ages <- seq(0,50,by=0.1)
vls <- viral_load_func(pars1, test_ages, TRUE)
probs <- prop_detectable(test_ages,pars1,vls)

fitted_dat <- tibble(age=test_ages,prob=probs)

p_detect_gamma <- ggplot()+
  geom_line(data=fitted_dat, aes(x=age,y=prob)) +
  geom_point(data=targeted_dat,aes(x=age,y=prob),col="red") +
  geom_vline(xintercept=assumed_incu_period, linetype="dashed") +
  export_theme +
  ylab("Proportion detectable") +
  xlab("Days since infection") 





######################################################
## C) PLOT DISTRIBUTION OF VIRAL LOADS
######################################################
vls <- viral_load_func(pars1, test_ages, TRUE)
line_dat <- tibble(mean_load=vls,t=test_ages)
violin_times <- seq(0,50,by=5)
hists <- matrix(ncol=10000,nrow=length(violin_times))
for(i in 1:length(violin_times)){
  omg <- viral_load_func(pars1, violin_times[i],TRUE)
  hists[i,] <- rnorm(10000,omg, pars1["obs_sd"])
}
hists[hists > pars1["intercept"]] <- pars1["intercept"]
hists <- reshape2::melt(hists)
colnames(hists) <- c("t","i","value")
hists$t <- violin_times[hists$t]
hists <- as_tibble(hists)

ribbon_times <- seq(0,50,by=0.1)
hists_ribbon <- matrix(ncol=10000,nrow=length(ribbon_times))
for(i in 1:length(ribbon_times)){
  omg <- viral_load_func(pars1, ribbon_times[i],TRUE)
  hists_ribbon[i,] <- rnorm(10000,omg, pars1["obs_sd"])
}
hists_ribbon[hists_ribbon > pars1["intercept"]] <- pars1["intercept"]
hists_ribbon <- reshape2::melt(hists_ribbon)
colnames(hists_ribbon) <- c("t","i","value")
hists_ribbon$t <- ribbon_times[hists_ribbon$t]
hists_ribbon <- as_tibble(hists_ribbon)

quants <- hists_ribbon %>% group_by(t) %>% summarize(lower_quant=quantile(value, 0.25),
                                                     upper_quant = quantile(value, 0.75))
p_ct_gamma <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=value,group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  geom_line(data=line_dat,aes(x=t,y=mean_load)) +
  geom_line(data=quants, aes(x=t,y=lower_quant),linetype="dotted",col="blue") +
  geom_line(data=quants, aes(x=t,y=upper_quant),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["intercept"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,60,by=5),trans="reverse") +
  ylab("Ct value") +
  xlab("Days since infection") +
  export_theme

p_vl_gamma <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=(40-value)/log2(10) + pars1["LOD"],group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  geom_line(data=line_dat,aes(x=t,y=(40-mean_load)/log2(10) + pars1["LOD"])) +
  geom_line(data=quants, aes(x=t,y=(40-lower_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  geom_line(data=quants, aes(x=t,y=(40-upper_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["LOD"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,15,by=1)) +
  ylab("log10 RNA copies / ml") +
  xlab("Days since infection") +
  export_theme

p_final_gamma <- (p_ct_gamma | p_vl_gamma) / p_detect_gamma


png("figs/gamma_eyeball_fit.png",width=8,height=7,units="in",res=300)
p_final_gamma
dev.off()


######################################################
## 2. FIT HINGE MODEL SUBJECTIVELY
######################################################
parTab_hinge <- read.csv("pars/partab_hinge.csv")

## Takes 2 days for viral loads to peak
parTab_hinge[parTab_hinge$names == "tshift","values"] <- 2
## Viral loads peak 4 days post infection
parTab_hinge[parTab_hinge$names == "desired_mode","values"] <- 4
## Dealing with log10 per ul
parTab_hinge[parTab_hinge$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
ages_observed <- c(-5, 0, 5, 10, 15, 20, 25, 30, 35,40,45,50) + assumed_incu_period
desired_probs <- c(0,1,0.9,0.72,0.49,0.29,0.14,0.04,0,0,0,0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_hinge <- parTab_hinge$values
names(pars_hinge) <- parTab_hinge$names

ages <- ages_observed

######################################################
## COST FUNCTION
cost_function_hinge <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1){
  pars1 <- pars_hinge
  ## Change true 0, gamma variance, observation variance and peak viral load
  pars1["true_0"] <- explore_pars[1]
  pars1["t_wane"] <- explore_pars[2]
  pars1["obs_sd"] <- explore_pars[3]
  pars1["viral_peak"] <- explore_pars[4]
  peak <- viral_load_func_hinge(pars1,pars1["desired_mode"],convert_ct=TRUE)
  cost1 <- 0
  if(use_ct_cost){
    cost1 <- (11.4 - qnorm(0.01,peak,pars1["obs_sd"]))^2
  }
  vls <- viral_load_func_hinge(pars1, ages, convert_ct=TRUE)
  obs <- prop_detectable(ages,pars1, vls)
  costs2 <- (obs-desired_probs)^2
  
  sum(costs2 + cost1*ct_cost_weight)
}

######################################################
## A) RUN OPTIM
######################################################
fit_hinge <- optim(c(-5,25,7,7), cost_function_hinge,control=list(maxit=10000,abstol=1e-8,reltol=1e-8),
                   use_ct_cost=TRUE, ct_cost_weight=1)

pars1 <- pars_hinge
pars1["true_0"] <- fit_hinge$par[1]
pars1["t_wane"] <- fit_hinge$par[2]
pars1["obs_sd"] <- fit_hinge$par[3]
pars1["viral_peak"] <- fit_hinge$par[4]


parTab_hinge$values <- pars1
write_csv(parTab_hinge,"pars/partab_hinge_fitted.csv")
######################################################
## B) PLOT AGAINST PROPORTION DETECTABLE
######################################################
test_ages <- seq(0,50,by=0.1)
vls <- viral_load_func_hinge(pars1, test_ages, TRUE)
probs <- prop_detectable(test_ages,pars1,vls)

fitted_dat <- tibble(age=test_ages,prob=probs)

p_detect_hinge <- ggplot()+
  geom_line(data=fitted_dat, aes(x=age,y=prob)) +
  geom_point(data=targeted_dat,aes(x=age,y=prob),col="red") +
  geom_vline(xintercept=assumed_incu_period, linetype="dashed") +
  export_theme +
  ylab("Proportion detectable") +
  xlab("Days since infection") 

######################################################
## C) PLOT DISTRIBUTION OF VIRAL LOADS
######################################################
vls <- viral_load_func_hinge(pars1, test_ages, TRUE)
line_dat <- tibble(mean_load=vls,t=test_ages)
violin_times <- seq(0,50,by=5)
hists <- matrix(ncol=10000,nrow=length(violin_times))
for(i in 1:length(violin_times)){
  omg <- viral_load_func_hinge(pars1, violin_times[i],TRUE)
  hists[i,] <- rnorm(10000,omg, pars1["obs_sd"])
}
hists[hists > pars1["intercept"]] <- pars1["intercept"]
hists <- reshape2::melt(hists)
colnames(hists) <- c("t","i","value")
hists$t <- violin_times[hists$t]
hists <- as_tibble(hists)

ribbon_times <- seq(0,50,by=0.1)
hists_ribbon <- matrix(ncol=10000,nrow=length(ribbon_times))
for(i in 1:length(ribbon_times)){
  omg <- viral_load_func_hinge(pars1, ribbon_times[i],TRUE)
  hists_ribbon[i,] <- rnorm(10000,omg, pars1["obs_sd"])
}
hists_ribbon[hists_ribbon > pars1["intercept"]] <- pars1["intercept"]
hists_ribbon <- reshape2::melt(hists_ribbon)
colnames(hists_ribbon) <- c("t","i","value")
hists_ribbon$t <- ribbon_times[hists_ribbon$t]
hists_ribbon <- as_tibble(hists_ribbon)

quants <- hists_ribbon %>% group_by(t) %>% summarize(lower_quant=quantile(value, 0.25),
                                              upper_quant = quantile(value, 0.75))

p_ct_hinge <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=value,group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  geom_line(data=line_dat,aes(x=t,y=mean_load)) +
  geom_line(data=quants, aes(x=t,y=lower_quant),linetype="dotted",col="blue") +
  geom_line(data=quants, aes(x=t,y=upper_quant),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["intercept"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,60,by=5),trans="reverse") +
  ylab("Ct value") +
  xlab("Days since infection") +
  export_theme

p_vl_hinge <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=(40-value)/log2(10) + pars1["LOD"],group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  geom_line(data=line_dat,aes(x=t,y=(40-mean_load)/log2(10) + pars1["LOD"])) +
  geom_line(data=quants, aes(x=t,y=(40-lower_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  geom_line(data=quants, aes(x=t,y=(40-upper_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["LOD"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,15,by=1)) +
  ylab("log10 RNA copies / ml") +
  xlab("Days since infection") +
  export_theme

p_final_hinge<- (p_ct_hinge | p_vl_hinge) / p_detect_hinge



png("figs/hinge_eyeball_fit.png",width=8,height=7,units="in",res=300)
p_final_hinge
dev.off()




######################################################
## 3. FIT 2-STAGE HINGE MODEL SUBJECTIVELY
######################################################
parTab_hinge2 <- read.csv("pars/partab_2hinge.csv")

## Dealing with log10 per ul
parTab_hinge2[parTab_hinge2$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
ages_observed <- c(-5, 0, 5, 10, 15, 20, 25, 30, 35,40,45,50) + assumed_incu_period
desired_probs <- c(0,1, 0.9,0.72,0.49,0.29,0.14,0.04,0,0,0,0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_hinge2 <- parTab_hinge2$values
names(pars_hinge2) <- parTab_hinge2$names

ages <- ages_observed

######################################################
## COST FUNCTION
cost_function_hinge2 <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1){
  pars1 <- pars_hinge2
  ## Change true 0, gamma variance, observation variance and peak viral load
  #pars1["true_0"] <- explore_pars[1]
  pars1["t_switch"] <- explore_pars[1]
  pars1["obs_sd"] <- explore_pars[2]
  pars1["viral_peak"] <- explore_pars[3]
  
  pars1["level_switch"] <- explore_pars[4]
  pars1["wane_rate2"] <- explore_pars[5]
  pars1["prob_detect"] <- explore_pars[6]
  
  peak <- viral_load_func_asymp(pars1,pars1["desired_mode"]+pars1["tshift"],convert_ct=TRUE)
  cost1 <- 0
  if(use_ct_cost){
    #cost1 <- (11.4 - qnorm(0.01,peak,pars1["obs_sd"]))^2
    cost1 <- (11.4 - qgumbel(0.01,peak,pars1["obs_sd"]))^2
  }
  vls <- viral_load_func_asymp(pars1, ages, convert_ct=TRUE)
  #obs <- prop_detectable_gumbell(ages,pars1, vls)
  obs <- prop_detectable(ages,pars1, vls, obs_model="gumbel", additional_detect_process = TRUE)
  costs2 <- (obs-desired_probs)^2
  sum(costs2 + cost1*ct_cost_weight)
}
upper_bounds <- c(30,10,25,5,5,1)
lhs_pars <- randomLHS(10000,6)
test_pars <- t(apply(lhs_pars,1, function(x) x*upper_bounds))

res <- apply(test_pars, 1, function(x) cost_function_hinge2(x))
start_par <- test_pars[which.min(res),]
######################################################
## A) RUN OPTIM
######################################################
fit_hinge2 <- optimx(c(12,5,9,3.5,30,0.05), cost_function_hinge2,
                     lower=c(0,0,0,3,0,0),
                     upper=c(25,10,25,5,50,1),
                   use_ct_cost=TRUE, ct_cost_weight=1)

pars1 <- pars_hinge2

pars1["t_switch"] <- as.numeric(fit_hinge2[1,1])
pars1["obs_sd"] <- as.numeric(fit_hinge2[1,2])
pars1["viral_peak"] <- as.numeric(fit_hinge2[1,3])
pars1["level_switch"] <- as.numeric(fit_hinge2[1,4])
pars1["wane_rate2"] <- as.numeric(fit_hinge2[1,5])
pars1["prob_detect"] <- as.numeric(fit_hinge2[1,6])
pars1["t_unit"] <- 1

######################################################
## B) PLOT AGAINST PROPORTION DETECTABLE
######################################################
parTab_hinge2$values <- pars1
write_csv(parTab_hinge2,"pars/partab_hinge2_fitted.csv")

test_ages <- seq(0,55,by=1)
vls <- viral_load_func_asymp(pars1, test_ages, TRUE)
probs <- prop_detectable(test_ages,pars1,vls,obs_model="gumbel",additional_detect_process = TRUE)

fitted_dat <- tibble(age=test_ages,prob=probs)

p_detect_hinge2 <- ggplot()+
  geom_line(data=fitted_dat, aes(x=age,y=prob)) +
  geom_point(data=targeted_dat,aes(x=age,y=prob),col="red") +
  geom_vline(xintercept=assumed_incu_period, linetype="dashed") +
  export_theme +
  ylab("Proportion detectable") +
  xlab("Days since infection") 
p_detect_hinge2
######################################################
## C) PLOT DISTRIBUTION OF VIRAL LOADS
######################################################
vls <- viral_load_func_asymp(pars1, test_ages, TRUE)
line_dat <- tibble(mean_load=vls,t=test_ages)
violin_times <- seq(0,50,by=5)
hists <- matrix(ncol=10000,nrow=length(violin_times))
for(i in 1:length(violin_times)){
  omg <- viral_load_func_asymp(pars1, violin_times[i],TRUE)
  #hists[i,] <- rnorm(10000,omg, pars1["obs_sd"])
  hists[i,] <- rgumbel(10000,omg, pars1["obs_sd"])
  
}
hists[hists > pars1["intercept"]] <- NA# pars1["intercept"]
hists <- reshape2::melt(hists)
colnames(hists) <- c("t","i","value")
hists$t <- violin_times[hists$t]
hists <- as_tibble(hists)

ribbon_times <- seq(0,50,by=0.1)
hists_ribbon <- matrix(ncol=10000,nrow=length(ribbon_times))
for(i in 1:length(ribbon_times)){
  omg <- viral_load_func_asymp(pars1, ribbon_times[i],TRUE)
  hists_ribbon[i,] <- rnorm(10000,omg, pars1["obs_sd"])
}
hists_ribbon[hists_ribbon > pars1["intercept"]] <- NA# pars1["intercept"]
hists_ribbon <- reshape2::melt(hists_ribbon)
colnames(hists_ribbon) <- c("t","i","value")
hists_ribbon$t <- ribbon_times[hists_ribbon$t]
hists_ribbon <- as_tibble(hists_ribbon)

quants <- hists_ribbon %>% group_by(t) %>% summarize(lower_quant=quantile(value, 0.25,na.rm=TRUE),
                                                     upper_quant = quantile(value, 0.75,na.rm=TRUE))

p_ct_hinge2 <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=value,group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  #geom_line(data=line_dat,aes(x=t,y=mean_load)) +
  #geom_line(data=quants, aes(x=t,y=lower_quant),linetype="dotted",col="blue") +
  #geom_line(data=quants, aes(x=t,y=upper_quant),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["intercept"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,60,by=5),trans="reverse") +
  ylab("Ct value") +
  xlab("Days since infection") +
  export_theme

p_vl_hinge2 <- ggplot() +
  geom_violin(data=hists,aes(x=t,y=(40-value)/log2(10) + pars1["LOD"],group=t),scale="width",trim=TRUE,draw_quantiles = c(0.025,0.5,0.975),fill="grey70")+
  #geom_line(data=line_dat,aes(x=t,y=(40-mean_load)/log2(10) + pars1["LOD"])) +
  #geom_line(data=quants, aes(x=t,y=(40-lower_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  #geom_line(data=quants, aes(x=t,y=(40-upper_quant)/log2(10) + pars1["LOD"]),linetype="dotted",col="blue") +
  geom_hline(yintercept=pars1["LOD"],linetype="dashed") +
  scale_y_continuous(breaks=seq(0,15,by=1)) +
  ylab("log10 RNA copies / ml") +
  xlab("Days since infection") +
  export_theme

p_final_hinge2 <- (p_ct_hinge2 | p_vl_hinge2) / p_detect_hinge2
p_final_hinge2


png("figs/hinge2_eyeball_fit.png",width=8,height=7,units="in",res=300)
p_final_hinge2
dev.off()




