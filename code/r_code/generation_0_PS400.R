########## 
##########
# This code contains some data preparation and model fitting for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-03-29
##########
##########

# set-up =======================================================================

library(here)
library(tidyverse)
library(rstan)
library(gdata)
library(bayesplot)
library(brms)
library(beepr)
library(parallel)
library(R2WinBUGS)
library(deSolve)
detectCores()


dir = 'C:/Users/brookson/Documents/Github/Schur-etal-Data/' #private data repo

growth_data = read_csv(paste0(dir, 'growth_data_gen0.csv'))
reproduction_data = read_csv(paste0(dir, 'offspring_data.csv'), 
                             guess_max = 15000)
unique(reproduction_data$treatment)

# get the data required ========================================================

# keep only the LFC and generation 0
reproduction_data = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'PS_400') 

# make the length into mm
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment == 'MP400') %>% #keep only control treatment
  filter(test == 1) %>% 
  filter(day < 23) %>%
  mutate(length_mm = length*0.001)  

# prep data ====================================================================

# create data matrix
l_y_obs_data = as.matrix(growth_data %>% 
                           filter(day < 22) %>% 
                           select(day, length_mm, replicate) %>% 
                           pivot_wider(id_cols = c(replicate, day, length_mm),
                                       names_from = replicate,
                                       values_from = length_mm) %>% 
                           select_if(~ !any(is.na(.))) %>% 
                           select(-day))

# identify replicates that didn't survive
reps_to_remove = reproduction_data %>% 
  select(day, replicate, offspring) %>% 
  group_by(day, replicate) %>% 
  filter(offspring == 'X')

# remove reps that didn't surviv
reproduction_data_x_removed = reproduction_data %>% 
  filter(!replicate %in% reps_to_remove$replicate)

# make numeric
reproduction_data_x_removed$offspring = 
  as.numeric(reproduction_data_x_removed$offspring)

# turn NA's into zeros
reproduction_data_x_removed$offspring[which(
  is.na(reproduction_data_x_removed$offspring))] = 0

# create data matrix
r_y_data = as.matrix(reproduction_data_x_removed %>% 
                       filter(day < 22) %>%  
                       pivot_wider(id_cols = c(replicate, day, offspring),
                                   names_from = replicate,
                                   values_from = offspring, 
                                   values_fn = max))
r_y_data = r_y_data[,2:ncol(r_y_data)]

#declare all variables
N_obs = 10
N_mis = 12
ii_obs = c(1, 3, 6, 8, 10, 13, 15, 17, 20, 22)
ii_mis = c(2, 4, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21) 
ll_init = 0.2
cq_init = 400
l_y_obs = l_y_obs_data[,1]#this is taking first replicate (column)
ts = 1:nrow(r_y_data)
r_y = r_y_data[,1] #this is taking first replciate (column )

gen_0_ps400_data = list(
  N_obs = N_obs, 
  N_mis = N_mis,
  ll_init = array(ll_init),
  cq_init = array(cq_init),
  ii_obs = ii_obs,
  ii_mis = ii_mis, 
  l_y_obs = l_y_obs,
  ts = ts,
  r_y = r_y)

# fit model ====================================================================
gen_0_ps400_onerep_fit = stan(file = 
                here('./code/stan_files/organism_costs_model_gen0_ps400.stan'),
                                data = gen_0_ps400_data,
                                chains = 4,
                                cores = 8,
                                warmup = 5000,
                                iter = 20000,
                                seed = 12,
                                verbose = TRUE,
                                #open_progress = TRUE,
                                control = list(adapt_delta = 0.9999)); beep(3)
saveRDS(gen_0_ps400_onerep_fit, 
        here('/output/intermediate-objects/gen_0_ps400_onerep_fit.RDS'))

# diagnose model fit
gen_0_ps400_summ = print(gen_0_ps400_onerep_fit, 
                       pars=c("theta_ll[1]", "theta_cq[1]", "cstar", "NEC",
                              "Lp", "Rm", "Lm", "tau_l", "tau_r"),
                       probs=c(0.1, 0.5, 0.9), digits = 3)

parms = c("theta_ll[1]", "theta_cq[1]", "cstar", "NEC",
           "Lp", "Rm", "Lm", "tau_l", "tau_r")
gen_0_fit_output1 = rstan::extract(gen_0_ps400_summ,
                                    permuted=TRUE,include=TRUE)

gen_0_control_onerep_fit_Trace = stan_trace(gen_0_ps400_onerep_fit,parms)
gen_0_control_onerep_fit_Dens = mcmc_dens(gen_0_ps400_onerep_fit,parms)
gen_0_control_onerep_fit_Overlay = mcmc_dens_overlay(gen_0_ps400_onerep_fit,parms)
gen_0_control_onerep_fit_Violin = mcmc_violin(gen_0_ps400_onerep_fit,parms,probs = c(0.1, 0.5, 0.9))
gen_0_control_onerep_fit_Pairs = mcmc_pairs(gen_0_ps400_onerep_fit,parms)





