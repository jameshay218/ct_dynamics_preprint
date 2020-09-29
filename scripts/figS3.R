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
reporting_frac <- 0.1
reporting_prob <- 0.1
rerun_mcmc <- TRUE

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
observed_indivs_flat <- simulate_reporting(observed_individuals, frac_report=frac_report,timevarying_prob=NULL,
                                           solve_times=times, symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_flat$sampled_individuals,
                                                      kinetics_pars=kinetics_pars,viral_load_func_ver="hinge2",obs_ver="gumbel",
                                                      additional_detect_process=TRUE)

## And add some symptomatic individuals
## Symptomatic surveillance
observed_indivs_symptom <- simulate_reporting(observed_individuals, frac_report=reporting_prob, 
                                              solve_times=times, symptomatic=TRUE)

simulated_viral_loads_symptom <- simulate_viral_loads_wrapper(observed_indivs_symptom$sampled_individuals,
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



plot_dat_symptom <- simulated_viral_loads_symptom %>% filter(confirmed_time <= 125) %>% mutate(week = floor(sampled_time/7)) %>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < kinetics_pars["intercept"]) 
plot_dat_normal <- simulated_viral_loads %>% mutate(week = floor(sampled_time/7)) %>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < kinetics_pars["intercept"]) 
plot_dat <- bind_rows(plot_dat_normal, plot_dat_symptom)

use_days <- plot_dat %>% 
  group_by(first_day) %>%
  tally() %>% 
  filter(n >= 50) %>% 
  pull(first_day)

subset_dat <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  group_by(first_day) %>% sample_n(pmin(50,n()))

plot_dat_normal_use <- plot_dat_normal %>%
  filter(first_day %in% use_days) %>%
  group_by(first_day) %>%
  summarize(median_ct=median(ct_obs))

p2 <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  ggplot() + 
  geom_violin(aes(x=first_day, y=ct_obs,group=week),scale="width",
              fill="grey70",alpha=0.75,color=NA) +
  #geom_jitter(aes(x=first_day, y=panther_Ct,group=week),width=2,height=0,
  #            fill="grey70",size=0.25,alpha=0.25) +
  #geom_dotplot(data=subset_dat,aes(x=first_day, y=ct_obs,group=first_day),binaxis="y",
  #             binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data= . %>% group_by(first_day) %>% summarize(median_ct=median(ct_obs)),
              aes(x=first_day,y=median_ct),col="purple",se=FALSE) +
  geom_smooth(data= plot_dat_normal_use,
              aes(x=first_day,y=median_ct),col="blue",se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(45, 5),expand=c(0,0)) +
  scale_x_continuous(limits=c(50,325),breaks=seq(0,365,by=25)) +
  geom_hline(yintercept=40,linetype="dashed") +
  geom_vline(xintercept=125,linetype="dashed",col="red") +
  theme_bw() +
  theme_classic()+
  xlab("Day of sample") +
  ylab("Ct value") +
  theme(legend.position="none",
        plot.title=element_blank(),
        panel.grid.minor =element_blank())
p2


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

main_p <- p2/(pA|pB|pC)

ggsave("figs/FigS3.pdf",p2,
       width=5,height=3,units="in")
ggsave("figs/FigS3.png",p2,
       width=5,height=3,units="in")

