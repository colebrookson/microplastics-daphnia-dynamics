########## 
##########
# This code contains some data preparation and model fitting for an analysis
# investigating the effects of microplastic pollution on Daphnia magna
# presented in Brookson et al. ()
##########
##########
# AUTHOR: Cole B. Brookson
# DATE OF CREATION: 2021-03-22
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

# keep only the LFC and generation 0
reproduction_data = reproduction_data %>% 
  dplyr::select(-`X1`) %>% 
  filter(treatment == 'LFC') 

# make the length into mm
growth_data = growth_data %>% 
  select(-`x1`) %>% 
  filter(treatment == 'LFC') %>% #keep only control treatment
  filter(test == 1) %>% 
  filter(day < 23) %>%
  mutate(length_mm = length*0.001)  

# keep only the first replicate for each and only 21 for reproduction
reproduction_data_meanreps = reproduction_data %>% 
  group_by(day) %>% 
  summarize(mean_offspring = mean(offspring, na.rm = TRUE))
  
# growth_data = growth_data %>% 
#   filter(day < 23) %>% 
#   group_by(day) %>% 
#   summarize(mean_length = mean(length_mm, na.rm = TRUE))

growth_repro_data = left_join(reproduction_data_meanreps, 
                              growth_data, 
                              by = 'day')
growth_repro_data[is.na(growth_repro_data['mean_offspring']),
                  'mean_offspring']=0

# add concentration (0)
growth_data$conc = 0
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
x_obs = 22
N_obs = 10
N_mis = 12
rep_r = ncol(r_y_data)
rep_l = ncol(l_y_obs_data)
ii_obs = c(1, 3, 6, 8, 10, 13, 15, 17, 20, 22)
ii_mis = c(2, 4, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21) 
ll_init = rep(0.2, rep_l)
l_y_obs = l_y_obs_data
ts = 1:nrow(r_y_data)
r_y = r_y_data

gen_0_control_data = list(
                          N_obs = N_obs, 
                          N_mis = N_mis,
                          rep_r = rep_r,
                          rep_l = rep_l,
                          ll_init = ll_init,
                          ii_obs = ii_obs,
                          ii_mis = ii_mis, 
                          l_y_obs = l_y_obs,
                          ts = ts,
                          r_y = r_y)

N_obs = 10
N_mis = 12
ii_obs = c(1, 3, 6, 8, 10, 13, 15, 17, 20, 22)
ii_mis = c(2, 4, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21) 
ll_init = 0.2
l_y_obs = l_y_obs_data[,1]
ts = 1:nrow(r_y_data)
r_y = r_y_data[,1]

gen_0_control_onerep_data = list(
  N_obs = N_obs, 
  N_mis = N_mis,
  ll_init = array(ll_init),
  ii_obs = ii_obs,
  ii_mis = ii_mis, 
  l_y_obs = l_y_obs,
  ts = ts,
  r_y = r_y)

# fit model ====================================================================
gen_0_control_fit = stan(file = 
                           here('./code/stan_files/organism_costs_model_control.stan'),
               data = gen_0_control_data,
               chains = 1,
               cores = 8,
               warmup = 5000,
               iter = 10000,
               seed = 12,
               verbose = TRUE,
               #open_progress = TRUE,
               control = list(adapt_delta = 0.9999)); beep(3)
gen_0_control_onerep_fit = stan(file = 
                         here('./code/stan_files/organism_costs_model_control_onerep.stan'),
                         data = gen_0_control_onerep_data,
                         chains = 4,
                         cores = 8,
                         warmup = 5000,
                         iter = 10000,
                         seed = 12,
                         verbose = TRUE,
                         #open_progress = TRUE,
                         control = list(adapt_delta = 0.9999)); beep(3)
saveRDS(gen_0_control_onerep_fit, here('/output/intermediate-objects/gen_0_fit_onerep.RDS'))
### diagnose model fit
gen_0_fit_summ = print(gen_0_control_onerep_fit, 
                     pars=c("theta_ll", "Lp", "Rm", "Lm", "tau_l", "tau_r"),
                     probs=c(0.1, 0.5, 0.9), digits = 3)

parms <- c("theta_ll[1]", "Lp", "Rm", "Lm", "tau_l", "tau_r")
gen_0_fit_output1 <- rstan::extract(gen_0_control_onerep_fit,permuted=TRUE,include=TRUE)

# Parms
gen_0_control_onerep_fit_Trace <- stan_trace(gen_0_control_onerep_fit,parms)
gen_0_control_onerep_fit_Dens <- mcmc_dens(gen_0_control_onerep_fit,parms)
gen_0_control_onerep_fit_Overlay <- mcmc_dens_overlay(gen_0_control_onerep_fit,parms)
gen_0_control_onerep_fit_Violin <- mcmc_violin(gen_0_control_onerep_fit,parms,probs = c(0.1, 0.5, 0.9))
gen_0_control_onerep_fit_Pairs <- mcmc_pairs(gen_0_control_onerep_fit,parms)
gen_0_control_onerep_fit_Trace
gen_0_control_onerep_fit_Dens
gen_0_control_onerep_fit_Overlay
gen_0_control_onerep_fit_Violin
gen_0_control_onerep_fit_Pairs

png(here('./output/model-fit-figs/gen_0_control_onerep_fit_Trace.png'))
gen_0_control_onerep_fit_Trace
dev.off()
ggsave(here('./output/model-fit-figs/gen_0_control_onerep_fit_Dens.png'),
       gen_0_control_onerep_fit_Dens, dpi = 200,
       width = 6, height = 6)
ggsave(here('./output/model-fit-figs/gen_0_control_onerep_fit_Overlay.png'),
       gen_0_control_onerep_fit_Overlay, dpi = 200, 
       width = 9, height = 6)
ggsave(here('./output/model-fit-figs/gen_0_control_onerep_fit_Violin.png'),
       gen_0_control_onerep_fit_Violin, dpi = 200, 
       width = 9, height = 6)
ggsave(here('./output/model-fit-figs/gen_0_control_onerep_fit_Pairs.png'),
       gen_0_control_onerep_fit_Pairs, dpi = 200, 
       width = 6, height = 6)


