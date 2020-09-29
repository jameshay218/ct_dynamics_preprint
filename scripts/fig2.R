library(tidyverse)
library(ggplot2)
library(doParallel)
library(patchwork)
library(coda)
library(ggthemes)
library(ggsci)
library(lazymcmc)
library(data.table)
library(deSolve)
library(odin)
library(extraDistr)
library(dde)

set.seed(1234)

main_wd <- "~/Documents/GitHub/ct_dynamics_preprint/"
setwd(main_wd)

## File management
savewd <- "results/figure2"
chainwd <- "chains/figure2"

source("code/model_funcs.R")
source("code/simulation_functions.R")
source("code/odin_funcs.R")
source("code/priors.R")
source("code/sim_data_funcs.R")

set.seed(1234)

population_n <- 6900000
times <- 0:365
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 
reporting_frac <- 0.2
reporting_prob <- 0.1
rerun_mcmc <- FALSE

##############################
## Model parameters
best_pars_fitting <- read.csv("pars/seir_fits_best_pars.csv")
best_pars <- as.numeric(best_pars_fitting[best_pars_fitting$institution == "NH 1",])
names(best_pars) <- colnames(best_pars_fitting)
best_pars["R0"] <- 2
best_pars["I0"] <- 100/population_n
best_pars["t0"] <- 0

kinetics_pars_dat <- read_csv("pars/partab_hinge2_fitted.csv")
kinetics_pars <- kinetics_pars_dat$values
names(kinetics_pars) <- kinetics_pars_dat$names
##############################3

## Simulate SEIR dynamics, incidence, growth rates, Rt etc
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,pars=best_pars, ver="odin",tstep=1)
incidence <- seir_dynamics$incidence/population_n


## Simulate onset times, confirmation delays etc
observed_individuals <- simulate_observations_wrapper(seir_dynamics$incidence,times=times,population_n=population_n)

## Simulate time-varying reporting **fraction**
## This represents a fraction of the population being tested each day
frac_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.00001,end_prob=0.001,
                                                            growth_rate=0.1,switch_point=120), ver="increase")
frac_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.0001,end_prob=0.00001,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Simulate time-varying reporting **probability**
## This gives the probability of an individual getting reported if they have symptoms
prob_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.01,end_prob=0.1,
                                                                   growth_rate=0.05,switch_point=250), ver="increase")
prob_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.1,end_prob=0.01,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Probability of getting reported
## NOTE - need to be careful this does not surpass 1, as each individual can only be observed once
1-prod(1-prob_report_increase$prob)

bind_rows(prob_report_decrease, prob_report_increase) %>%
  ggplot() +
  geom_line(aes(x=t,y=prob,col=ver)) +
  theme_bw() +
  scale_y_continuous(limits=c(0,max(prob_report_increase$prob)*1.1)) +
  ylab("Probability of report") +
  xlab("Time") +
  theme(legend.position=c(0.8,0.5))

## Random surveillance
observed_indivs_flat <- simulate_reporting(observed_individuals, frac_report=reporting_frac,timevarying_prob=NULL,
                                           solve_times=times, symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_flat$sampled_individuals,
                                                      kinetics_pars=kinetics_pars,viral_load_func_ver="hinge2",obs_ver="gumbel",
                                                      additional_detect_process=TRUE)

 
coeff <- 30000
p1A <- ggplot(seir_dynamics$seir_outputs)+ 
  geom_line(aes(x=step,y=inc),col="red")+
  geom_line(aes(x=step,y=Rt*coeff),col="forestgreen") +
scale_y_continuous(expand=c(0,0),limits=c(0,80000),
                   sec.axis=sec_axis(~.*1/coeff, name="Rt")) +
  scale_x_continuous(limits=c(50,325),breaks=seq(0,365,by=25)) +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("New infections") 
p1A



plot_dat <- simulated_viral_loads %>% mutate(week = floor(sampled_time/7)) %>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < kinetics_pars["intercept"]) 
use_days <- plot_dat %>% 
  group_by(first_day) %>%
  tally() %>% 
  filter(n >= 50) %>% 
  pull(first_day)

subset_dat <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  group_by(first_day) %>% sample_n(pmin(50,n()))

p2 <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  ggplot() + 
  geom_violin(aes(x=first_day, y=ct_obs,group=week),scale="width",
              fill="grey70",alpha=0.75,color=NA) +
  #geom_jitter(aes(x=first_day, y=panther_Ct,group=week),width=2,height=0,
  #            fill="grey70",size=0.25,alpha=0.25) +
  geom_dotplot(data=subset_dat,aes(x=first_day, y=ct_obs,group=first_day),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data= . %>% group_by(first_day) %>% summarize(median_ct=median(ct_obs)),
              aes(x=first_day,y=median_ct),col="blue",se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(45, 5),expand=c(0,0)) +
  scale_x_continuous(limits=c(50,325),breaks=seq(0,365,by=25)) +
  geom_hline(yintercept=40,linetype="dashed") +
  theme_bw() +
  theme_classic()+
  xlab("Day of sample") +
  ylab("Ct value") +
  theme(legend.position="none",
        plot.title=element_blank(),
        panel.grid.minor =element_blank()) +
  labs(tag="B")
p2

p1A/p2

seir_weekly <- seir_dynamics$seir_outputs %>% 
  as_tibble() %>% 
  rename(sampled_time=step) %>%
  mutate(week=floor(sampled_time/7))%>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  group_by(week, first_day) %>%
  summarize(mean_rt=mean(Rt))
  

combined_dat <- simulated_viral_loads %>% 
  mutate(week=floor(sampled_time/7))%>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < kinetics_pars["intercept"]) %>%
  #filter(sampled_time >= min(use_days) & sampled_time <= max(use_days)) %>%
  #group_by(sampled_time) %>% 
  #filter(ct_obs < kinetics_pars["intercept"]) %>% 
  group_by(week, first_day) %>%
  summarize(median_ct=median(ct_obs),
            skew_ct=moments::skewness(ct_obs),
            n=n()) %>%
  #rename(date=first_day) %>%
  full_join(seir_weekly) %>%
  filter(n >= 50)


pA <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=skew_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=skew_ct),size=2) +
  geom_vline(xintercept=1,linetype="dashed") +
  scale_y_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  theme_classic() +
  ylab("Skewness of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="C")

pB <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=median_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=median_ct),size=2)+
  geom_vline(xintercept=1,linetype="dashed") +
  theme_classic() +
  scale_y_continuous(trans="reverse",limits=c(35,28),breaks=seq(28,35,by=1)) +
  ylab("Median of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="D")

pC <- ggplot(combined_dat) +
  geom_point(aes(x=skew_ct,y=median_ct,col=mean_rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1,
                        limits=c(0,2),
                        guide=guide_colorbar(title="Rt",
                                            barwidth=1,ticks=FALSE,barheight=4))+
  scale_x_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  scale_y_continuous(trans="reverse",limits=c(35,28),breaks=seq(28,35,by=1)) +
  theme_classic() +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  theme(legend.position=c(0.2,0.8),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6)) +
  labs(tag="E")

main_p <- p1A/p2/(pA|pB|pC)

ggsave("figs/Fig2.pdf",main_p,
       width=8,height=8,units="in")
ggsave("figs/Fig2.png",main_p,
       width=8,height=8,units="in")

#############################################
## MCMC - Simulation recovery
#############################################
ages <- seq(0,35,by=1)
prior_func_use <- prior_func_hinge2
prior_table1 <- data.frame(par=c("beta","obs_sd","viral_peak","wane_rate2","t_switch","level_switch","prob_detect"),
                           prior_ver=c("normal","normal","normal","normal","normal","normal","beta"),
                           mean=c(0,6,9,30,16,3.7,0.1),
                           sd=c(0.4,1,1,3,2,0.25,0.1)/2)
all_priors <- generate_all_priors(prior_table1)

ts <- seq(50,350,by=5)

## MCMC management
library(doParallel)
nchains <- 2
n_clusters <- 12
cl <- makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

mcmcPars <- c("iterations"=20000,"popt"=0.44,"opt_freq"=100,
              "thin"=10,"adaptive_period"=10000,"save_block"=1000)
mcmcPars2 <- c("iterations"=50000,"popt"=0.234,"opt_freq"=100,
              "thin"=10,"adaptive_period"=20000,"save_block"=1000)
if(!file.exists(savewd)) dir.create(savewd,recursive=TRUE)
if(!file.exists(chainwd)) dir.create(chainwd,recursive=TRUE)

ns_hinge <- numeric(length(ts))
parTab <- read.csv("pars/partab_hinge2_fitted.csv",stringsAsFactors=FALSE)
parTab[parTab$names != "beta", "fixed"] <- 1
#for(i in seq_along(ts)){
ns_hinge <- foreach(i=seq_along(ts),.packages = c("lazymcmc","extraDistr","ggthemes","tidyverse","ggsci")) %dopar% {
  source(paste0(main_wd,"/code/model_funcs.R"))
#i <- 5
  time <- ts[i]
  filename <- paste0("sim_","hinge2","_","gumbel","_",time)
  sim_dat1 <- simulated_viral_loads %>% filter(sampled_time ==time & ct_obs < 40) %>% pull(ct_obs)
  n_used <- 0
  if(length(sim_dat1) > 0) { 
    n_used <- length(sim_dat1)
    if(rerun_mcmc){
      real_fitting_func(ages, sim_dat1,parTab, filename,nchains,paste0(main_wd,"/",chainwd), paste0(main_wd,"/",savewd), 
                        mcmcPars, 
                        mcmcPars_multivariate=mcmcPars2,PRIOR_FUNC = prior_func_use,
                        all_priors=all_priors, plot_par_lines=TRUE,
                        viral_load_func_use=viral_load_func_asymp,
                        obs_model="gumbel", additional_detect_process = TRUE)
    }
  }
  n_used
}

n_hinge_dat_hinge <- simulated_viral_loads %>% 
  filter(ct_obs < 40) %>% 
  group_by(sampled_time) %>% 
  tally() %>%
  rename(ts=sampled_time)

## plot inferred growth rates against SEIR and Cts
all_betas <- NULL
prob_growth <- rep(NA,length(ts))
setwd(main_wd)
for(time in ts){
  print(time)
  sim_dat1 <- simulated_viral_loads %>% filter(sampled_time ==time & ct_obs < 40) %>% pull(ct_obs)
  if(length(sim_dat1) > 0){
    filename <- paste0("sim_","hinge2","_","gumbel","_",time)
    chainwd1 <- paste0(chainwd, "/",filename)
    chains_tmp <- load_mcmc_chains(chainwd1, parTab, TRUE, 10, mcmcPars["adaptive_period"], TRUE)
    chain <- chains_tmp[[2]]
    betas <- tibble(beta=chain[,"beta"],t=time)
    n_growing <- chain[chain[,"beta"] > 0,]
    n_sampd <- 0
    if(is.null(nrow(n_growing))) {
      if(!is.null(n_growing)) {
        n_sampd <- 1
      } else {
        n_sampd <- 0
      }
    } else {
      n_sampd <- nrow(n_growing)
    }
    prob_growth[match(time, ts)] <- n_sampd/nrow(chain)
    all_betas <- rbind(all_betas, betas)
  }
}
growth_probs <- tibble(ts=ts, p=prob_growth)

beta_estimates <- all_betas %>%# mutate(beta = beta*beta_scale + max(incidence)) %>%
  group_by(t) %>% summarise(lower=quantile(beta,0.025),
                            upper=quantile(beta,0.975),
                            median=quantile(beta,0.5)) %>% 
  rename(ts=t) %>%
  left_join(n_hinge_dat_hinge) %>%
  rename(t=ts)

seir_dat <- tibble(t=times+1,y=incidence)
max_beta <- max(beta_estimates$upper)
max_inc <- max(incidence)
inc_scale <- max_beta/max_inc

GR_daily <- log(incidence[2:350]/incidence[1:349])
GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
GR_full <- NULL
lastday <- 35
for (i in (lastday+1):length(GR_daily)) {
  #for (i in 1:length(GR_daily)) {
  end_index <- i-1
  start_index <- max(1, (i-lastday))
  GR_full <- c(GR_full,mean(GR_daily[start_index:end_index], na.rm=TRUE))
}
GR <- data.frame(Days=(lastday+1):length(GR_daily), GR=GR_full, DailyGR=GR_daily[(lastday+1):length(GR_daily)])

inc_scale <- 75
y_max <- 0.5
p1 <- ggplot() +
  geom_line(data=seir_dat,aes(x=t,y=y*inc_scale + -y_max),size=1,col="#000000")+
  geom_line(data=GR,aes(x=Days,y=GR),col="#E69F00",size=1) +
  geom_point(data=beta_estimates%>% filter(t %in% seq(50, 350, by=10)),aes(x=t,y=median,col=n),size=0.5) +
  geom_errorbar(data=beta_estimates %>% filter(t %in% seq(50, 350, by=10)),aes(x=t,ymin=lower,ymax=upper,col=n),
                width=0.5) +
  #geom_line(data=GR,aes(x=Days,y=DailyGR)) +
  scale_y_continuous(breaks=seq(-y_max,y_max,by=0.2), 
                     limits=c(-y_max,y_max),
                     sec.axis = sec_axis(~ . *(1/inc_scale) - -y_max/inc_scale,name="Incidence per capita")) +
  #coord_cartesian(ylim=c(-0.3, 0.3)) +
  scale_x_continuous(limits=range(ts)) +
  scale_color_gradient(limits=c(0,800),low="blue",high="red")+#,breaks=seq(0,3000,by=500)) +
  guides(color=guide_colorbar("Detectable samples", direction="horizontal",title.position="top", 
                              barwidth=6,barheight=1,ticks=FALSE)) +
  geom_hline(yintercept=0.0,linetype="dashed") +
  export_theme +
  ylab("Growth rate") +
  xlab("Days since start of epidemic") +
  theme(
    axis.title.x=element_text(size=8),
    axis.text.x=element_text(size=6),
    axis.title.y=element_text(size=8),
    axis.text.y=element_text(size=6),

    legend.position="none",
        legend.background = element_rect(fill="white",color="white"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.size= unit(0.5, "cm"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        plot.tag = element_text(vjust=-1.5, size=12)) +
  labs(tag="A")

merged_dat <- GR %>% as_tibble() %>% rename(t=Days) %>% left_join(beta_estimates)

p2 <- ggplot(merged_dat %>% filter(n > 25)) + 
  geom_rect(data=data.frame(xmin=0.08,xmax=0,ymin=0, ymax=0.6),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="green",alpha=0.05) +
  geom_rect(data=data.frame(xmin=-0.08,xmax=0,ymin=0, ymax=-0.6),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="green",alpha=0.05) +
  geom_rect(data=data.frame(xmin=0.08,xmax=0,ymin=0, ymax=-0.6),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="orange",alpha=0.05) +
  geom_rect(data=data.frame(xmin=-0.08,xmax=0,ymin=0, ymax=0.6),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="orange",alpha=0.05) +
  geom_abline(slope=-1,col="black",linetype="dashed") +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  #geom_point(aes(x=GR,y=median,col=n),size=1) +
  geom_point(aes(x=GR,y=median,col=n),size=1) +
  geom_errorbar(aes(x=GR,ymin=lower,ymax=upper,col=n)) +
  #geom_pointrange(aes(x=GR,y=median,ymin=lower,ymax=upper,col=n),size=0.5,alpha=0.75) +
  scale_color_gradient(limits=c(0,800),low="blue",high="red")+#,breaks=seq(0,3000,by=500)) +
  guides(color=guide_colorbar("Detectable samples", direction="vertical",
                              title.position = "left",
                              barwidth=1,barheight=6,ticks=FALSE)) +
  scale_x_continuous(trans="reverse") +
  coord_cartesian(ylim=c(-0.5,0.5),xlim=c(0.06,-0.06)) +
  export_theme +
  theme(legend.position="none",
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        legend.title = element_text(angle=90),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        plot.tag = element_text(vjust=-1.5,size=12)) +
  xlab("True growth rate") +
  ylab("Estimated growth rate") +
  labs(tag="B")

all_betas %>%
  filter(t %in% seq(50,350,by=10)) %>%
  mutate(is_growing=ifelse(beta > 0, 1,0)) %>%
  group_by(t) %>% 
  summarise(prob_growing=sum(is_growing)/n()) %>%
  ggplot() +
  geom_line(aes(x=t,y=prob_growing)) +
  geom_hline(yintercept=0.5,linetype="dashed") +
  export_theme

t_test <- c(150,200,250)

all_preds_exp <- NULL
for(index in seq_along(t_test)){
  ti <- t_test[index]
  tmp_beta_tmp <- all_betas %>% filter(t == ti)
  tmp_exp <- matrix(ncol=length(ages),nrow=nrow(tmp_beta_tmp))
  for(i in 1:nrow(tmp_beta_tmp)){
    beta1 <- tmp_beta_tmp[i,] %>% pull(beta)
    y <- exp(-ages*beta1)
    tmp_exp[i,] <- y/max(y)
  }
  tmp_exp <- as.data.frame(tmp_exp)
  colnames(tmp_exp) <- ages
  tmp_exp$t <- ti
  all_preds_exp[[index]] <- tmp_exp
}
all_preds_exp <- do.call("bind_rows",all_preds_exp)

exp_dat <- all_preds_exp %>% pivot_longer(-t) %>%
  mutate(name=as.numeric(name)) %>% 
  group_by(t, name) %>%
  summarize(median=mean(value),
            lower=quantile(value,0.025),
            lower_mid=quantile(value,0.25),
            upper=quantile(value,0.975),
            upper_mid=quantile(value,0.75)) %>%
  mutate(t_name=ifelse(t==t_test[1],paste0("Day ",t_test[1], " (growth)"), as.character(t)),
         t_name=ifelse(t==t_test[2],paste0("Day ",t_test[2], " (peak)"), t_name),
         t_name=ifelse(t==t_test[3],paste0("Day ",t_test[3], " (decline)"), t_name))
        

p4 <- ggplot(exp_dat) +
  geom_ribbon(aes(x=name,ymin=lower,ymax=upper),alpha=0.25,fill="#D55E00") +
  geom_ribbon(aes(x=name,ymin=lower_mid,ymax=upper_mid),alpha=0.5,fill="#D55E00") +
  geom_line(aes(x=name,y=median),col="#D55E00") +
  scale_x_continuous(trans="reverse") +
  facet_wrap(~t_name,scales="free_y",ncol=1) +
  ylab("Relative probability of infection") +
  xlab("Days since infection") +
  scale_y_continuous(expand=c(0,0)) +
  export_theme +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        plot.tag = element_text(vjust=-1.5,size=12),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6)) +
  labs(tag="D")


quants_all <- NULL
for(index in seq_along(t_test)) {
  filename <- paste0("sim_","hinge2","_","gumbel","_",t_test[index])
  chainwd1 <- paste0(chainwd, "/",filename)
  chains_tmp <- load_mcmc_chains(chainwd1, parTab, FALSE, 10, mcmcPars["adaptive_period"], TRUE)
  chain <- chains_tmp[[2]]
  real_dat <- simulated_viral_loads %>% filter(sampled_time ==t_test[index] & ct_obs < 40) %>% pull(ct_obs)
  test_cts <- seq(0,parTab[parTab$names == "intercept","values"]-1,by=1)
  f <- create_lik_func(parTab, real_dat, ages, PRIOR_FUNC=NULL, solve_ver="model",test_cts=test_cts,
                       viral_load_func_use=viral_load_func_asymp,
                       obs_model="gumbel", additional_detect_process=TRUE)
  samps <- sample(1:nrow(chain),200)
  
  vls <- matrix(0, nrow=length(samps),ncol=length(test_cts))
  for(i in seq_along(samps)){
    pars <- get_index_pars(chain, samps[i])
    tmp <- f(pars)
    tmp <- tmp/sum(tmp)
    vls[i,] <- tmp
  }
  
  quants <- as.data.frame(t(apply(vls, 2, function(x) quantile(x, c(0.025,0.5,0.975)))))
  colnames(quants) <- c("lower","median","upper")
  quants$age <- test_cts
  quants$time_test <- t_test[index]
  quants_all[[index]] <- quants
}
quants_all_comb <- do.call("bind_rows", quants_all)

dat_use_last_p <- simulated_viral_loads %>% 
  filter(sampled_time %in% t_test & ct_obs < 40) %>% 
  rename(t=sampled_time) %>%
  mutate(t_name=ifelse(t==t_test[1],paste0("Day ",t_test[1], " (growth)"), as.character(t)),
         t_name=ifelse(t==t_test[2],paste0("Day ",t_test[2], " (peak)"), t_name),
         t_name=ifelse(t==t_test[3],paste0("Day ",t_test[3], " (decline)"), t_name))

quants_all_comb <- quants_all_comb %>%
  rename(t=time_test) %>%
  mutate(t_name=ifelse(t==t_test[1],paste0("Day ",t_test[1], " (growth)"), as.character(t)),
         t_name=ifelse(t==t_test[2],paste0("Day ",t_test[2], " (peak)"), t_name),
         t_name=ifelse(t==t_test[3],paste0("Day ",t_test[3], " (decline)"), t_name))

p3 <- ggplot(data=dat_use_last_p) + 
  geom_histogram(aes(x=ct_obs, y=..density..),fill="grey70",col="black",binwidth=1) +
  geom_ribbon(data=quants_all_comb,aes(x=age+0.5,ymin=lower,ymax=upper),fill="blue",alpha=0.25) +
  geom_line(data=quants_all_comb,aes(y=median,x=age+0.5)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Ct value") +
  ylab("Density") +
  scale_x_continuous(breaks=seq(0,40,by=5)) +
  facet_wrap(~t_name,ncol=1) +
  export_theme +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      plot.tag = element_text(vjust=-1.5,size=12),
      axis.title.x=element_text(size=8),
      axis.text.x=element_text(size=6),
      axis.title.y=element_text(size=8),
      axis.text.y=element_text(size=6)) +
  labs(tag="C")
p3

fig4 <- ((p1/p2) | p3 | p4) + plot_layout(widths=c(2,1,1))

p2_b <- p2 +  
  theme(legend.position="right")

ggsave("figs/fig4.png",fig4, width=8,height=5,units="in")
ggsave("figs/fig4.pdf",fig4, width=8,height=5,units="in")
ggsave("figs/fig4_legend.pdf",p2_b, width=8,height=5,units="in")
``