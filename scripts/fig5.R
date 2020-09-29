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

set.seed(7)

main_wd <- "~/Documents/GitHub/ct_inference_preprint/"
setwd(main_wd)

source("code/model_funcs.R")
source("code/simulation_functions.R")
source("code/odin_funcs.R")
source("code/sim_data_funcs.R")

set.seed(1234)

population_n <- 6900000
times <- 0:365
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 
reporting_frac <- 0.1
reporting_prob <- 0.1


##############################
## Model parameters
best_pars_fitting <- read.csv("pars/seir_fits_best_pars.csv")
best_pars <- as.numeric(best_pars_fitting[best_pars_fitting$institution == "NH 1",2:(ncol(best_pars_fitting)-1)])
names(best_pars) <- colnames(best_pars_fitting)[2:(ncol(best_pars_fitting)-1)]
best_pars["R0"] <- 2
best_pars["I0"] <- 100/population_n
best_pars["t0"] <- 0

kinetics_pars_dat <- read_csv("pars/partab_hinge2_fitted.csv")
kinetics_pars <- kinetics_pars_dat$values
names(kinetics_pars) <- kinetics_pars_dat$names
##############################3

## Simulate SEIR dynamics, incidence, growth rates, Rt etc
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,pars=best_pars, ver="odin",tstep=1)

## Simulate onset times, confirmation delays etc
observed_individuals <- simulate_observations_wrapper(seir_dynamics$incidence,times=times,population_n=population_n,symp_frac = 1)

## Simulate time-varying reporting **fraction**
## This represents a fraction of the population being tested each day
frac_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.00001,end_prob=0.001,
                                                            growth_rate=0.05,switch_point=250), ver="increase")

observed_indivs_flat <- simulate_reporting(observed_individuals, frac_report=0.2,timevarying_prob=NULL,
                                           solve_times=times, symptomatic=FALSE)
observed_indivs_increasing <- simulate_reporting(observed_individuals, timevarying_prob=frac_report_increase, 
                                                 solve_times=times, symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_increasing$sampled_individuals,
                                                      kinetics_pars=kinetics_pars,viral_load_func_ver="hinge2",obs_ver="gumbel",
                                                      additional_detect_process=TRUE)

simulated_viral_loads_flat <- simulate_viral_loads_wrapper(observed_indivs_flat$sampled_individuals,
                                                      kinetics_pars=kinetics_pars,viral_load_func_ver="hinge2",obs_ver="gumbel",
                                                      additional_detect_process=TRUE)

dat_flat <- simulated_viral_loads_flat %>% filter(ct_obs < 40) %>% ungroup() %>% 
  sample_n(simulated_viral_loads %>% filter(ct_obs < 40) %>% nrow) %>%
  group_by(sampled_time) %>% tally()
dat_increasing <- simulated_viral_loads %>% filter(ct_obs < 40) %>% 
  group_by(sampled_time) %>% tally()

p1 <- ggplot() +
  geom_line(data=dat_flat,aes(x=sampled_time,y=n)) +
  geom_line(data=dat_increasing,aes(x=sampled_time,y=n),col="red")

p2 <- ggplot(seir_dynamics$seir_outputs) +
  geom_line(aes(x=step,y=Rt))
p1/p2
